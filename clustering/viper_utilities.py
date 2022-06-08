import os, shutil, gzip
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

# CheckV Copyright (c) 2020, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy). All rights reserved.
# ===>
def parse_blast(handle):
    for line in handle:
        r = line.split()
        yield {
            "qname": r[0],
            "tname": r[1],
            "pid": float(r[2]),
            "len": float(r[3]),
            "qcoords": sorted([int(r[6]), int(r[7])]),
            "tcoords": sorted([int(r[8]), int(r[9])]),
            "qlen": float(r[-2]),
            "tlen": float(r[-1]),
            "evalue": float(r[-4]),
        }


def yield_alignment_blocks(handle):
    # init block with 1st record
    key, alns = None, None
    for aln in parse_blast(handle):
        if aln["qname"] == aln["tname"]:
            continue
        key = (aln["qname"], aln["tname"])
        alns = [aln]
        break
    # loop over remaining records
    for aln in parse_blast(handle):
        # skip self hits
        if aln["qname"] == aln["tname"]:
            continue
        # extend block
        elif (aln["qname"], aln["tname"]) == key:
            alns.append(aln)
        # yield block and start new one
        else:
            yield alns
            key = (aln["qname"], aln["tname"])
            alns = [aln]
    yield alns


def prune_alns(alns, min_len=0, min_evalue=1e-3, min_pid=90):
    # remove short aligns
    alns = [
        aln
        for aln in alns
        if aln["len"] >= min_len
        and aln["evalue"] <= min_evalue
        and aln["pid"] >= min_pid
    ]
    return alns


def compute_ani(alns):
    return round(
        sum(a["len"] * a["pid"] for a in alns) / sum(a["len"] for a in alns), 2
    )


def compute_cov(alns):

    # merge qcoords
    coords = sorted([a["qcoords"] for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        # overlapping, update start coord
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        # non-overlapping, append to list
        else:
            nr_coords.append([start, stop])
    # compute qry_cov
    alen = sum([stop - start + 1 for start, stop in nr_coords])
    qcov = round(100.0 * alen / alns[0]["qlen"], 2)

    # merge tcoords
    coords = sorted([a["tcoords"] for a in alns])
    nr_coords = [coords[0]]
    for start, stop in coords[1:]:
        # overlapping, update start coord
        if start <= (nr_coords[-1][1] + 1):
            nr_coords[-1][1] = max(nr_coords[-1][1], stop)
        # non-overlapping, append to list
        else:
            nr_coords.append([start, stop])
    # compute qry_cov
    alen = sum([stop - start + 1 for start, stop in nr_coords])
    tcov = round(100.0 * alen / alns[0]["tlen"], 2)

    return qcov, tcov


# <===


def anicalc(blast_input, pid=90, length=0):
    with open(blast_input, "r") as input:
        anicalc_dict = {}
        index = 0
        for alns in yield_alignment_blocks(input):
            alns = prune_alns(alns, min_pid=pid, min_len=length)
            if len(alns) == 0:
                continue
            index += 1
            qname, tname = alns[0]["qname"], alns[0]["tname"]
            ani = compute_ani(alns)
            qcov, tcov = compute_cov(alns)
            anicalc_dict[index] = {
                "qname": qname,
                "tname": tname,
                "num_alns": len(alns),
                "pid": ani,
                "qcov": qcov,
                "tcov": tcov,
            }
        anicalc_df = pd.DataFrame.from_dict(anicalc_dict, orient="index")
        return anicalc_df


def parse_seqs(path):
    handle = gzip.open(path) if path.split(".")[-1] == "gz" else open(path)
    id = next(handle).split()[0][1:]
    seq = ""
    for line in handle:
        if line[0] == ">":
            yield id, seq
            id = line.split()[0][1:]
            seq = ""
        else:
            seq += line.rstrip()
    yield id, seq
    handle.close()


def aniclust(
    fasta, anicalc_df, minlength=0, min_qcov=0, min_tcov=85, min_ani=95  # outname,
):
    # list seqs, sorted by length
    print("\nreading sequences...")
    seqs = {}

    for index, r in enumerate(parse_seqs(fasta)):
        id, seq = r
        if len(seq) < minlength:
            continue
        else:
            seqs[id] = len(seq)

    seqs = [x[0] for x in sorted(seqs.items(), key=lambda x: x[1], reverse=True)]
    print("%s sequences retained from fna" % len(seqs))

    # store edges
    print("\nstoring edges...")
    num_edges = 0
    edges = {x: [] for x in seqs}

    df2 = anicalc_df.loc[
        (anicalc_df["qname"] != anicalc_df["tname"])
        & (anicalc_df["qname"].isin(edges.keys()))
        & (anicalc_df["tname"].isin(edges.keys()))
        & (anicalc_df["pid"] >= min_ani)
        & (anicalc_df["qcov"] >= min_qcov)
        & (anicalc_df["tcov"] >= min_tcov),
        ["qname", "tname"],
    ]
    for row in df2.to_dict("records"):
        edges[row["qname"]].append(row["tname"])
        num_edges += 1

    print("%s edges retained from blastani" % num_edges)
    print("%s edges currently stored" % sum([len(_) for _ in edges.values()]))

    # cluster
    print("\nclustering...")

    clust_to_seqs = {}
    seq_to_clust = {}
    # loop over seqs in sorted order
    for seq_id in seqs:
        # seq already been assigned; cant be centroid
        if seq_id in seq_to_clust:
            continue
        # seq is centroid for new cluster
        else:
            # add self to cluster
            clust_to_seqs[seq_id] = [seq_id]
            seq_to_clust[seq_id] = seq_id
            # update with cluster members
            for mem_id in edges[seq_id]:
                if mem_id not in seq_to_clust:
                    clust_to_seqs[seq_id].append(mem_id)
                    seq_to_clust[mem_id] = seq_id

    print("%s total clusters" % len(clust_to_seqs))

    # write
    print("\nwriting clusters...")
    # keep = set()
    # for seq_id, mem_ids in clust_to_seqs.items():
    #    keep.add(seq_id)
    # with open(outname, 'w') as out:
    #    for seq_id, mem_ids in clust_to_seqs.items():
    #        keep.add(seq_id)
    #        out.write(seq_id + '\t' + ','.join(mem_ids)+'\n')

    return clust_to_seqs


def clustering(fasta, output, threads, pid=95, cov=85, returnDict=False):
    if os.path.exists("blastdb"):
        shutil.rmtree("blastdb")
    os.mkdir("blastdb")
    makedb = NcbimakeblastdbCommandline(
        dbtype="nucl", input_file=fasta, out="blastdb/" + output + "_db"
    )
    blastn = NcbiblastnCommandline(
        query=fasta,
        db="blastdb/" + output + "_db",
        outfmt="6 std qlen slen",
        max_target_seqs=10000,
        perc_identity=90,
        num_threads=threads,
        out=output + ".out",
    )
    makedb()
    blastn()

    anicalc_df = anicalc(output + ".out")

    aniclust_dict = aniclust(fasta, anicalc_df, min_ani=pid, min_tcov=cov)

    with open(output + ".fasta", "w") as f:
        for seq in SeqIO.parse(fasta, "fasta"):
            if seq.id in aniclust_dict.keys():
                SeqIO.write(seq, f, "fasta")
    shutil.rmtree("blastdb")
    os.remove(output + ".out")

    if returnDict:
        return aniclust_dict
