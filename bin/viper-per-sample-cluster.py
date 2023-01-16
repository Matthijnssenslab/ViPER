#!/usr/bin/env python3

import argparse
import logging
import os
import shutil
import sys
import tempfile
import textwrap
from pathlib import Path

import checkv
import pandas as pd
import pybedtools
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from clustering import viper_utilities as vu

logger = vu.get_logger()


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description=f"""This script takes the unclustered contigs from ViPER's triple assembly as input, 
    runs CheckV (https://bitbucket.org/berkeleylab/checkv) to split proviruses in 'host' and 'viral' contigs, and
    subsequently clusters the contigs. This script also ensures that contigs with contamination, duplications or contigs that
    are longer than expected, as assessed by CheckV, do not end up as (wrong) cluster representatives.
    """,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="fasta",
        required=True,
        type=str,
        metavar="PATH",
        help="Fasta file.",
    )
    parser.add_argument(
        "-d",
        "--database",
        dest="db",
        required=True,
        type=str,
        metavar="PATH",
        help="CheckV database.",
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
        "-l",
        "--length",
        dest="length",
        type=int,
        metavar="INT",
        help="Minimum length of the host or viral region. (default: %(default)s)",
        default=1000,
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
    parser.add_argument(
        "--keep-bed",
        dest="bed",
        action="store_true",
        help="Keep BED files with viral and host regions. (default: %(default)s)",
        default=False,
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        help="Keep temporary files for debugging. (default: %(default)s)",
        default=False,
    )
    return vars(parser.parse_args())


def make_bed(contamination, output, minlength, keepBed=False):
    """Function to generate BED files/dataframes from CheckV contamination output."""
    # Read in contamination file from CheckV
    df = pd.read_csv(contamination, sep="\t", dtype=object)

    # Define columns to keep from the file
    col_list = ["contig_id", "region_types", "region_lengths", "region_coords_bp"]

    # Subset dataframe with only proviruses and the desired columns
    df1 = df.loc[(df["provirus"] == "Yes"), col_list]

    # Define possible contamination combo's
    contamination_list = ["viral,host", "host,viral", "host,viral,host"]

    # Check if there is another contamination possibility than defined
    if not df1.loc[~df1["region_types"].isin(contamination_list)].empty:
        logger.warninf("Other types of contamination:")
        print(*df1["contig_id"].values, sep="\n")

    # Subset the dataframe with only contaminations in 'contamination_list' and continue
    df1 = df1.loc[df1["region_types"].isin(contamination_list)]

    # Check if there are proviruses found, otherwise return none
    if df1.empty:
        logger.info(f"No proviruses found.")
        return None, None

    # Create new df with viral/host regions split up
    df2 = pd.concat(
        [
            df1["contig_id"],
            df1[col_list[1]].str.split(",", expand=True),
            df1[col_list[2]].str.split(",", expand=True),
            df1[col_list[3]].str.split(",", expand=True),
        ],
        axis=1,
    ).reset_index(drop=True)

    # create list of dataframes with all contamination regions
    data = []
    for i in range(3):
        try:
            df2[i]
        except:
            if i == 0:
                logger.info("No contamination in contigs.")
        else:
            tmpdf = pd.concat([df2["contig_id"], df2[i]], axis=1)
            tmpdf.columns = col_list
            data.append(tmpdf)
    df3 = pd.concat(data).sort_index().reset_index(drop=True)

    # Split the region coordinates in two columns and name them 'start'/'end'
    df4 = pd.concat(
        [
            df3.drop("region_coords_bp", axis=1),
            df3["region_coords_bp"].str.split("-", expand=True),
        ],
        axis=1,
    )
    df4.rename(columns={0: "start", 1: "end"}, inplace=True)

    # Convert to numeric
    df4[["region_lengths", "start", "end"]] = df4[
        ["region_lengths", "start", "end"]
    ].apply(pd.to_numeric)
    # convert from CheckV to BED coords
    df4[["start", "end"]] -= 1

    # Extract all viral regions from dataframe that are larger than the given minimum length
    viraldf = df4.loc[
        (df4["region_types"] == "viral") & (df4["region_lengths"] >= minlength),
        ["contig_id", "region_lengths", "start", "end"],
    ]
    # Rename the contig id with info on the viral region length
    viraldf.loc[:, "bed_name"] = [
        x.replace("_cov", "v" + str(y) + "_cov").replace("_length", "v_length")
        for x, y in viraldf[["contig_id", "region_lengths"]].to_numpy()
    ]
    viraldf.drop("region_lengths", inplace=True, axis=1)

    # Extract all host regions from dataframe that are larger than the given minimum length
    hostdf = df4.loc[
        (df4["region_types"] == "host") & (df4["region_lengths"] >= minlength),
        ["contig_id", "region_lengths", "start", "end"],
    ]
    # Rename the contig id with info on the host region length
    hostdf.loc[:, "bed_name"] = [
        x.replace("_cov", "h" + str(y) + "_cov").replace("_length", "h_length")
        for x, y in hostdf[["contig_id", "region_lengths"]].to_numpy()
    ]
    hostdf.drop("region_lengths", inplace=True, axis=1)

    # Keep BED files for debugging purposes
    if keepBed:
        viraldf.to_csv(output + "_viral.bed", sep="\t", header=False, index=False)
        hostdf.to_csv(output + "_host.bed", sep="\t", header=False, index=False)

    return viraldf, hostdf


def bedtools(bed, fasta):
    """Function to run bedtools with BED file from string on fasta sequences."""
    a = pybedtools.BedTool(bed, from_string=True)
    a = a.sequence(fi=fasta, name=True)
    return open(a.seqfn).read()


def bedtools_fasta_dict(fasta_text):
    """Function to return a dictionary with fasta sequences from bedtools function."""
    fasta_dict = {}  # make dictionary
    for i in fasta_text.split("\n"):
        if i.startswith(">"):
            key = i.split(":")[0]
            fasta_dict[key] = ""
        elif i == "":
            continue
        else:
            fasta_dict[key] = i
    return fasta_dict


def quality_summary_selection(checkv_summary):
    """Function to select sequences with duplication/high k-mer warnings as assessed by CheckV."""
    # Read in CheckV's quality_summary file
    df = pd.read_csv(checkv_summary, sep="\t", dtype={"warnings": "str"})
    # Fill NA in warning column with empty string
    df["warnings"] = df["warnings"].fillna(
        ""
    )  # TO DO: find faster way to get rid of Nan

    include = set()
    exclude = set()
    # Define warnings
    warnings = [
        "high kmer_freq may indicate large duplication",
        "contig >1.5x longer than expected genome length",
    ]
    # Change quality_summary df in dictionary
    df_dict = df.to_dict("records")
    for row in df_dict:
        # Select any provirus contig with duplication/longer than expected warning and add it to the 'exclude' set
        if row["provirus"] == "Yes" or any(x in row["warnings"] for x in warnings):
            exclude.add(row["contig_id"])

        # Select all contigs with warnings and add it to the 'include' set
        if any(x in row["warnings"] for x in warnings):
            include.add(row["contig_id"])

    return include, exclude


def write_fasta(dictionary, name):
    """Function to write fasta file from dictionary."""
    wrapper = textwrap.TextWrapper(width=60)
    with open(name, "w") as fh:
        for k, v in dictionary.items():
            if k.startswith(">"):
                fh.write(k + "\n")
                text_list = wrapper.wrap(text=v)
                for i in text_list:
                    fh.write(i + "\n")
            else:
                logger.error(
                    "Dictionary key does not start with >, your fasta file might be misformatted."
                )
                sys.exit(1)


def biopython_fasta(dictionary):
    """Function to write fasta file from biopython dictionary."""
    biopython_dict = {}
    for k, v in dictionary.items():
        k = k.lstrip(">")
        biopython_dict[k] = SeqRecord(Seq(v), id=k, name=k, description="")
    return biopython_dict


def main():
    # Parse script arguments
    args = parse_arguments()
    fasta = args["fasta"]
    minlength = args["length"]
    bed = args["bed"]
    debug = args["debug"]
    output = args["output"]
    output_name = Path(args["output"]).name
    output_dir = Path(args["output"]).parent

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Make tmp file for CheckV and other calculations
    tmpdir = tempfile.mkdtemp(dir=output_dir)
    # if os.path.exists("tmp_clustering"):
    #    shutil.rmtree("tmp_clustering")
    #
    # os.mkdir("tmp_clustering")

    # Define CheckV arguments
    checkv_arguments = {
        "input": fasta,
        "db": args["db"],
        "output": tmpdir,
        "threads": args["threads"],
        "restart": True,
        "quiet": False,
        "remove_tmp": True,
    }

    logger.info(f"Running CheckV...")

    # Run CheckV
    checkv.end_to_end.main(checkv_arguments)

    # Define CheckV output paths
    contamination = os.path.join(tmpdir, "contamination.tsv")
    qsummary = os.path.join(tmpdir, "quality_summary.tsv")

    logger.info(
        f"Excluding contigs with contamination, longer than expected and duplication issues..."
    )

    # Return sets of contigs to include/exclude
    include, exclude = quality_summary_selection(qsummary)

    logger.newline()
    logger.info(f"Making BED for host and viral regions...")

    # Make viral and host BED dfs for provirus contigs
    viralbed, hostbed = make_bed(contamination, output, minlength, keepBed=bed)

    # Define dictionaries
    hdict = {}
    vdict = {}
    clean_v_dict = {}

    # When there are proviruses in the data
    if viralbed is not None and hostbed is not None:
        logger.newline()
        logger.info(f"Splitting host sequence from viral contigs...")
        pysam.faidx(fasta)
        # Make fasta for host/viral regions
        host = bedtools(hostbed.to_csv(header=None, index=False, sep="\t"), fasta)
        viral = bedtools(viralbed.to_csv(header=None, index=False, sep="\t"), fasta)

        # Add them to a dictionary (key: contig_id, value: sequence)
        hdict = bedtools_fasta_dict(host)
        vdict = bedtools_fasta_dict(viral)

        # Remove provirus full contig names from include set and add their new names (with viral region info)
        # Also add the new names to a 'viral exclude' set
        viral_exclude = set()
        for row in viralbed.to_dict("records"):
            if row["contig_id"] in include:
                include.remove(row["contig_id"])
                include.add(row["bed_name"])
                viral_exclude.add(row["bed_name"])

        # Remove provirus contigs with duplication/longer than expected issues from provirus fasta dictionary
        clean_v_dict = {
            k: v for k, v in vdict.items() if k.lstrip(">") not in viral_exclude
        }

        # Remove fasta index file
        os.remove(fasta + ".fai")

    # Open unclustered fasta file
    fasta_seqs = {}
    with open(fasta, "r") as fh:
        lines = fh.readlines()
        for line in lines:
            if line.startswith(">"):
                key = line.strip()
                fasta_seqs[key] = ""
            else:
                fasta_seqs[key] += line.strip()

    # Remove sequences that are in the exclude set
    clean_fasta_dict = {
        k: v for k, v in fasta_seqs.items() if k.lstrip(">") not in exclude
    }

    # Combine the host, clean provirus and clean fasta dictionaries
    cluster_dict = {**clean_v_dict, **hdict, **clean_fasta_dict}
    write_fasta(cluster_dict, os.path.join(tmpdir, output_name + "_to_cluster.fasta"))

    # Generate reinclude fasta
    if not include:
        logger.newline()
        logger.info(
            f"No duplication issues in sequences, reinclude fasta is not generated."
        )
    else:
        logger.newline()
        logger.info(
            f"Writing fasta file with contigs to re-include after cross-sample clustering..."
        )
        # Combine general fasta with the viral regions fasta dictionary
        inclv_dict = {**vdict, **fasta_seqs}
        # Extract all sequences (with issues) to reinclude
        reinclude_dict = {
            k: v for k, v in inclv_dict.items() if k.lstrip(">") in include
        }
        write_fasta(reinclude_dict, output + "_re-include.fasta")

    logger.newline()

    # Cluster remaining sequences
    vu.clustering(
        os.path.join(tmpdir, output_name + "_to_cluster.fasta"),
        output + "_clustered",
        args["threads"],
        args["pid"],
        args["cov"],
    )

    if not debug:
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    main()
