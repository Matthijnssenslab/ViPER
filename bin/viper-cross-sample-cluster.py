#!/usr/bin/env python3

import os, argparse, textwrap, shutil, sys
import pandas as pd
from clustering import viper_utilities as vu
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description=f"""
        Script to cluster contigs across a study and reinclude contigs that had contamination issues.
    """,
    )
    parser.add_argument(
        "-c",
        "--cluster-fasta",
        dest="cluster_fasta",
        required=True,
        type=str,
        metavar="PATH",
        help="Fasta file.",
    )
    parser.add_argument(
        "-r",
        "--reinclude-fasta",
        dest="reinclude_fasta",
        required=True,
        type=str,
        metavar="PATH",
        help="Fasta file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        metavar="STR",
        help="Output name for files.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        metavar="INT",
        help="Number of threads to use.",
        default=1,
    )
    parser.add_argument(
        "--min-identity",
        dest="pid",
        type=int,
        metavar="INT",
        help="Minimum average nucleotide identity (ANI) for sequences to be clustered.",
        default=95,
    )
    parser.add_argument(
        "--min-coverage",
        dest="cov",
        type=int,
        metavar="INT",
        help="Minimum coverage %% of the shortest sequence that should be covered before clustering.",
        default=85,
    )
    return vars(parser.parse_args())


def main():
    args = parse_arguments()
    cluster_fasta = args["cluster_fasta"]
    reinclude_fasta = args["reinclude_fasta"]
    output = args["output"]
    threads = args["threads"]

    clust_seqs = vu.clustering(cluster_fasta, output, threads, returnDict=True)

    makedb = NcbimakeblastdbCommandline(
        dbtype="nucl", input_file=output + ".fasta", out="blastdb/" + output + "_db"
    )
    makedb()

    blastn = NcbiblastnCommandline(
        query=reinclude_fasta,
        db="blastdb/" + output + "_db",
        outfmt="6 std qlen slen",
        max_target_seqs=10000,
        perc_identity=90,
        num_threads=threads,
        out=output + "_reinclude.out",
    )
    blastn()

    anicalc_df = vu.anicalc(output + "_reinclude.out")
    anicalc_dict = anicalc_df.to_dict("records")

    qcov85 = set()
    tcov85 = set()
    singletons = set()

    for row in anicalc_dict:
        if (
            row["pid"] >= args["pid"]
            and row["qcov"] >= args["cov"]
            and row["qname"] not in qcov85
        ):
            qcov85.add(row["qname"])

    for row in anicalc_dict:
        if (
            row["pid"] >= args["pid"]
            and row["tcov"] >= args["cov"]
            and row["qname"] not in qcov85.union(tcov85)
        ):
            tcov85.add(row["qname"])

    reinclude_sequences = SeqIO.parse(reinclude_fasta, "fasta")

    for fasta in reinclude_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if name not in qcov85.union(tcov85, singletons):
            singletons.add(name)

    fasta_length_dictionary = {}
    fasta_sequences = SeqIO.parse(output + ".fasta", "fasta")
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_length_dictionary[name] = len(sequence)

    qcov_dict = {}
    for contig in qcov85:
        data = []
        for row in anicalc_dict:
            if (
                row["qname"] == contig
                and row["pid"] >= args["pid"]
                and row["qcov"] >= args["cov"]
            ):
                data.append(row["tname"])
        subset = {key: fasta_length_dictionary[key] for key in data}
        qcov_dict[contig] = max(subset, key=subset.get)

    tcov_dict = {}
    for contig in tcov85:
        data = []
        for row in anicalc_dict:
            if (
                row["qname"] == contig
                and row["pid"] >= args["pid"]
                and row["tcov"] >= args["cov"]
            ):
                data.append(row)
        df = pd.DataFrame(data)
        df.sort_values(by=["pid"], inplace=True, ascending=False)
        qcov_sum = df["qcov"].sum()
        if qcov_sum < 30:
            tcov_dict[contig] = contig
            continue
        while qcov_sum > 110:
            df = df[:-1]
            qcov_sum = df["qcov"].sum()
        tcov_dict[contig] = df["tname"].tolist()

    single_dict = {}
    for contig in singletons:
        single_dict[contig] = contig

    reinclude = {**qcov_dict, **tcov_dict, **single_dict}

    for k, v in reinclude.items():
        if isinstance(v, list):
            for i in v:
                clust_seqs[i].append(k)
        elif v not in clust_seqs:
            clust_seqs[v] = list(k)
        else:
            clust_seqs[v].append(k)

    with open(output + "_" + args["pid"] + "-" + args["cov"] + ".fasta", "w") as f:
        reinclude_sequences = SeqIO.parse(reinclude_fasta, "fasta")  # generator
        for fasta in reinclude_sequences:
            if fasta.id in clust_seqs.keys():
                SeqIO.write(fasta, f, "fasta")

        fasta_sequences = SeqIO.parse(output + ".fasta", "fasta")
        for fasta in fasta_sequences:
            if fasta.id in clust_seqs.keys():
                SeqIO.write(fasta, f, "fasta")

    with open(output + "_cluster_representatives.txt", "w") as out:
        for seq_id, mem_ids in clust_seqs.items():
            out.write(seq_id + "\t" + ",".join(mem_ids) + "\n")

    shutil.rmtree("blastdb")
    os.remove(output + "_reinclude.out")


if __name__ == "__main__":
    main()
