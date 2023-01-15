#!/usr/bin/env python3

import argparse
import logging
import os
import shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline

from clustering import viper_utilities as vu


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description=f"""
        Script to cluster contigs across a study and reinclude contigs that had contamination issues as assessed by the ViPER per sample clustering.
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
        help="Number of threads to use. (default: %(default)s)",
        default=1,
    )
    parser.add_argument(
        "--min-identity",
        dest="pid",
        type=int,
        metavar="INT",
        help="Minimum average nucleotide identity (ANI) for sequences to be clustered. (default: %(default)s)",
        default=95,
    )
    parser.add_argument(
        "--min-coverage",
        dest="cov",
        type=int,
        metavar="INT",
        help="Minimum coverage %% of the shortest sequence that should be covered before clustering. (default: %(default)s)",
        default=85,
    )
    return vars(parser.parse_args())


def main():
    args = parse_arguments()
    cluster_fasta = args["cluster_fasta"]
    reinclude_fasta = args["reinclude_fasta"]
    threads = args["threads"]
    output = args["output"]
    output_name = Path(args["output"]).name
    output_dir = Path(args["output"]).parent

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    logger = vu.get_logger()

    clust_seqs = vu.clustering(cluster_fasta, output, threads, returnDict=True)

    logger.newline()
    logger.info(
        f"Calculating ANI and coverage of sequences to be reincluded against clustered sequences..."
    )
    makedb = NcbimakeblastdbCommandline(
        dbtype="nucl",
        input_file=output + ".fasta",
        out=os.path.join("blastdb", output_name + "_db"),
    )
    makedb()

    blastn = NcbiblastnCommandline(
        query=reinclude_fasta,
        db=os.path.join("blastdb", output_name + "_db"),
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

    logger.newline()
    logger.info(f"Add sequences with possible duplications...")
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

    logger.info(f"Add possible chimeric sequences...")
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

    logger.info(f"Add singleton sequences...")
    single_dict = {}
    for contig in singletons:
        single_dict[contig] = contig

    reinclude = {**qcov_dict, **tcov_dict, **single_dict}

    for k, v in reinclude.items():
        if isinstance(v, list):
            for i in v:
                clust_seqs[i].append(k)
        elif v not in clust_seqs:
            clust_seqs[v] = [k]
        else:
            clust_seqs[v].append(k)

    logger.newline()
    logger.info(f"Write fasta file with clustered sequences...")
    with open(
        output + "_" + str(args["pid"]) + "-" + str(args["cov"]) + ".fasta", "w"
    ) as f:
        reinclude_sequences = SeqIO.parse(reinclude_fasta, "fasta")  # generator
        for fasta in reinclude_sequences:
            if fasta.id in clust_seqs.keys():
                SeqIO.write(fasta, f, "fasta")

        fasta_sequences = SeqIO.parse(output + ".fasta", "fasta")
        for fasta in fasta_sequences:
            if fasta.id in clust_seqs.keys():
                SeqIO.write(fasta, f, "fasta")

    logger.info(
        f"Write file with clusters and their respective cluster representatives..."
    )
    with open(output + "_cluster_representatives.txt", "w") as out:
        for seq_id, mem_ids in clust_seqs.items():
            out.write(seq_id + "\t" + ",".join(mem_ids) + "\n")

    shutil.rmtree("blastdb")
    os.remove(output + "_reinclude.out")


if __name__ == "__main__":
    main()
