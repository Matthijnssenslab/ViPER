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
import genomad
import pandas as pd
import pybedtools
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from viper_clustering import viper_utilities as vu

logger = vu.get_logger()


def parse_arguments():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.MetavarTypeHelpFormatter,
        description=f"""This script takes the unclustered contigs from ViPER's triple assembly as input,
    runs geNomad to find proviruses (https://github.com/apcamargo/genomad), splits the proviruses in
    'host' and 'viral' contigs and subsequently clusters the contigs. This script also runs CheckV
    (https://bitbucket.org/berkeleylab/checkv) to ensure that contigs with contamination, duplications
    or contigs that are longer than expected do not end up as (wrong) cluster representatives.
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
        "-o",
        "--output",
        dest="output",
        required=True,
        type=str,
        metavar="STR",
        help="Output name for files.",
    )
    parser.add_argument(
        "--min-identity",
        dest="pid",
        type=float,
        metavar="FLOAT",
        help="Minimum average nucleotide identity (ANI) for sequences to be clustered. [0-1] (default: %(default)s).",
        default=0.95,  # change from 95 to 99 to prevent over-clustering during per-sample clustering
    )
    parser.add_argument(
        "--min-coverage",
        dest="cov",
        type=float,
        metavar="FLOAT",
        help="Minimum coverage %% of the shortest sequence that should be covered before clustering. [0-1] (default: %(default)s).",
        default=0.85,  # change from 85 to 99 to prevent over-clustering during per-sample clustering
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
        "--checkv-db",
        dest="checkv_db",
        type=str,
        metavar="PATH",
        help="CheckV database.",
    )
    parser.add_argument(
        "--genomad-db",
        dest="genomad_db",
        type=str,
        metavar="PATH",
        help="geNomad database.",
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
        "--sensitivity-marker-search",
        dest="sens_marker",
        type=float,
        metavar="FLOAT",
        help="Sensitivity of MMseqs2 marker search against geNomad markers DB. Higher values will annotate more proteins, but the search will be slower and consume more memory. x>=0.0 (default: %(default)s)",
        default=4.2,
    )
    parser.add_argument(
        "--evalue-marker-search",
        dest="eval_marker",
        type=float,
        metavar="FLOAT",
        help="Maximum accepted E-value in the MMseqs2 marker search (default: %(default)s)",
        default=0.001,
    )
    parser.add_argument(
        "--splits",
        dest="splits",
        type=int,
        metavar="INT",
        help="Split the data for the MMseqs2 search. Higher values will reduce memory usage, but will make the search slower. If the MMseqs2 search is failing, try to increase the number of splits. x>=0 (default: %(default)s)",
        default=0,
    )
    parser.add_argument(
        "--crf-threshold",
        dest="ct",
        type=float,
        metavar="FLOAT",
        help="Minimum gene-level score to flag a provirus gene using the conditional random field model. Lower values will result in longer proviruses but will increase the probability of host genes being flagged as part of proviruses. 0.0<=x<=1.0 (default: %(default)s)",
        default=0.4,
    )
    parser.add_argument(
        "--marker-threshold",
        dest="mt",
        type=float,
        metavar="FLOAT",
        help="Minimum total virus marker score allowed for proviruses that do not encode integrases or are not located at scaffold edges. Lower values will increase the sensitivity but reduce the precision of the provirus identification procedure (default: %(default)s)",
        default=12.0,
    )
    parser.add_argument(
        "--marker-threshold-integrase",
        dest="mti",
        type=float,
        metavar="FLOAT",
        help="Minimum total virus marker score allowed for proviruses that encode integrases (default: %(default)s)",
        default=8,
    )
    parser.add_argument(
        "--marker-threshold-edge",
        dest="mte",
        type=float,
        metavar="FLOAT",
        help="Minimum total virus marker score allowed for proviruses that are located at scaffold edges (default: %(default)s)",
        default=8,
    )
    parser.add_argument(
        "--max-integrase-distance",
        dest="mid",
        type=int,
        metavar="INT",
        help="Maximum allowed distance between provirus boundaries and the integrases used for boundary extension. x>=0 (default: %(default)s)",
        default=10000,
    )
    parser.add_argument(
        "--max-trna-distance",
        dest="mtd",
        type=int,
        metavar="INT",
        help="Maximum allowed distance between provirus boundaries and the tRNAs used for boundary extension. x>=0 (default: %(default)s)",
        default=5000,
    )
    parser.add_argument(
        "--sensitivity-integrase-search",
        dest="sens_integrase",
        type=float,
        metavar="FLOAT",
        help="Sensitivity of MMseqs2 integrase search during geNomad proviruses identification. Higher values will annotate more proteins, but the search will be slower and consume more memory. x>=0.0 (default: %(default)s)",
        default=8.2,
    )
    parser.add_argument(
        "--evalue-integrase-search",
        dest="eval_integrase",
        type=float,
        metavar="FLOAT",
        help="Maximum accepted E-value in the MMseqs2 integrase search (default: %(default)s)",
        default=0.001,
    )

    parser.add_argument(
        "--keep-bed",
        dest="keep_bed",
        action="store_true",
        help="Keep BED files with viral and host regions. (default: %(default)s)",
        default=False,
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        help="Keep temporary files for debugging purposes. (default: %(default)s)",
        default=False,
    )
    return vars(parser.parse_args())


def proviruses_bed_coordinates(proviruses, minlength, output):
    """Function to generate the viral BED file from geNomad find_proviruses output."""
    # Read in find_proviruses file from geNomad
    viraldf = pd.read_csv(proviruses, sep="\t", dtype=object)

    # Check if there are proviruses found, otherwise return none
    if viraldf.empty:
        logger.info(f"No proviruses found.")
        return None

    viraldf[["start", "end", "length"]] = viraldf[["start", "end", "length"]].apply(
        pd.to_numeric, downcast="integer"
    )

    # Convert from Genomad to BED coords (only the start coordinate needs to be adjusted)
    viraldf["start"] -= 1

    # Keep only the proviral regions that are longer than the given minimum length
    viraldf = viraldf.loc[
        (viraldf["length"] >= minlength),
        ["source_seq", "start", "end", "length"],
    ]

    # Rename the contig id with info on the viral region length
    viraldf.loc[:, "bed_name"] = [
        x.replace("_cov", "V" + str(y) + "_cov").replace("_length", "V_length")
        for x, y in viraldf[["source_seq", "length"]].to_numpy()
    ]
    # Check for duplicates in the "bed_name" column (viral regions of same length)
    if viraldf["bed_name"].duplicated().any():
        logger.warning(
            f"Two (or more) viral regions from the same contig have identical length, resulting in identical contig names."
        )

    viraldf.drop("length", inplace=True, axis=1)

    # Write the viraldf to a BED file
    viraldf.to_csv(output + "_viral.bed", sep="\t", header=False, index=False)

    return viraldf


def source_seq_lengths_function(viraldf, output):
    """Function to generate a txt file with the lengths of the source sequences. It is needed to calculate the BED coordinates of the host regions"""
    df = pd.DataFrame(viraldf["source_seq"])

    # Create a "length" column and extract the length from the contig name
    df["length"] = df["source_seq"].str.split("_").apply(lambda x: x[3])

    # Remove duplicate source sequences (in case of multiple proviral regions in the same contig)
    df = df.drop_duplicates(subset=["source_seq"])[["source_seq", "length"]]

    df.to_csv(output + "_source_seq_lengths.txt", sep="\t", header=False, index=False)


def host_bed_coordinates(
    viralbed, source_seq_lengths, fasta, minlength, output, keepBed=False
):
    """Function to generate the host BED file using bedtools complement."""
    # Generate a panda dataframe with the BED coordinates of the host regions
    hostbed = pybedtools.BedTool().complement(i=viralbed, g=source_seq_lengths)
    hostdf = pd.read_table(hostbed.fn, names=["source_seq", "start", "end"])

    # Generate a sub-fasta file using the BED coordinates and the original fasta file
    hostfasta = hostbed.sequence(fi=fasta)
    fasta_file = open(hostfasta.seqfn).read()

    # Calculate the lengths of sequences in the output_fasta file and add to the hostdf
    lengths = []
    for i in fasta_file.split("\n"):
        if not i.startswith(">") and len(i) > 0:
            lengths.append(len(i))
    hostdf["length"] = lengths

    # Keep only the proviral regions that are longer than the given minimum length
    hostdf = hostdf[hostdf["length"] >= minlength]

    # Rename the contig id with info on the viral region length
    hostdf.loc[:, "bed_name"] = [
        x.replace("_cov", "H" + str(y) + "_cov").replace("_length", "H_length")
        for x, y in hostdf[["source_seq", "length"]].to_numpy()
    ]
    # Check for duplicates in the "bed_name" column (host regions of same length)
    if hostdf["bed_name"].duplicated().any():
        logger.warning(
            f"Two (or more) host regions from the same contig have identical length, resulting in identical contig names."
        )

    hostdf.drop("length", inplace=True, axis=1)

    # Keep BED file for debugging purposes
    if keepBed:
        hostdf.to_csv(output + "_host.bed", sep="\t", header=False, index=False)
    else:
        os.remove(output + "_viral.bed")
    return hostdf


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


def lists_proviruses_warnings(checkv_summary, viraldf):
    """Function to select sequences that are proviruses and/or have duplication/high k-mer warnings as assessed by CheckV."""
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
    if viraldf is not None:
        # Create a list of all the unique contigs names in the "source_seq" column of viraldf
        proviruses_source_seq_names = set(viraldf["source_seq"])

        for row in df_dict:
            # create exclude list: containing contigs that are proviruses or have CheckV warnings
            if row["contig_id"] in proviruses_source_seq_names or any(
                x in row["warnings"] for x in warnings
            ):
                exclude.add(row["contig_id"])

            # create include list: containing contigs with CheckV warnings
            if any(x in row["warnings"] for x in warnings):
                include.add(row["contig_id"])
    else:
        for row in df_dict:
            if any(x in row["warnings"] for x in warnings):
                exclude.add(row["contig_id"])
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


def run_checkv(fasta, checkv_db, threads, tmpdir):
    # Define CheckV arguments
    checkv_arguments = {
        "input": fasta,
        "db": checkv_db,
        "output": tmpdir,
        "threads": threads,
        "restart": True,
        "quiet": False,
        "remove_tmp": True,
    }

    logger.info(f"Running CheckV...")

    # Run CheckV
    checkv.end_to_end.main(checkv_arguments)


def run_genomad(
    fasta,
    tmpdir,
    genomad_db,
    threads,
    sens_marker,
    eval_marker,
    splits,
    ct,
    mt,
    mti,
    mte,
    mid,
    mtd,
    sens_integrase,
    eval_integrase,
):
    logger.info(f"Running geNomad...")

    # Run geNomad annotate
    genomad.annotate.main(
        input_path=Path(fasta),
        output_path=Path(tmpdir),
        database_path=Path(genomad_db),
        use_minimal_db=False,
        restart=True,
        threads=threads,
        verbose=True,
        conservative_taxonomy=False,
        sensitivity=sens_marker,
        evalue=eval_marker,
        splits=splits,
        cleanup=True,
    )

    # Run geNomad find-proviruses
    genomad.find_proviruses.main(
        input_path=Path(fasta),
        output_path=Path(tmpdir),
        database_path=Path(genomad_db),
        cleanup=False,
        restart=True,
        skip_integrase_identification=False,
        skip_trna_identification=False,
        threads=threads,
        verbose=True,
        crf_threshold=ct,
        marker_threshold=mt,
        marker_threshold_integrase=mti,
        marker_threshold_edge=mte,
        max_integrase_distance=mid,
        max_trna_distance=mtd,
        sensitivity=sens_integrase,
        evalue=eval_integrase,
    )


def exclude_sequences_for_clustering(
    fasta, output, output_name, tmpdir, minlength, keep_bed, debug
):
    # If debugging option is set also keep the bed files, regardless of the original value
    if debug:
        keep_bed = True

    # Define CheckV output paths
    qsummary = os.path.join(tmpdir, "quality_summary.tsv")

    # Define geNomad output paths
    genomad_find_proviruses_output_folder = str(Path(fasta).stem) + "_find_proviruses"
    genomad_tsv_output_file = str(Path(fasta).stem) + "_provirus.tsv"
    proviruses = os.path.join(
        tmpdir, genomad_find_proviruses_output_folder, genomad_tsv_output_file
    )

    logger.newline()
    logger.info(f"Making BED for host and viral regions...")

    # Make viral and host BED dfs for provirus contigs
    viralbed = proviruses_bed_coordinates(proviruses, minlength, output)
    if viralbed is not None:
        source_seq_lengths_function(viralbed, output)
        hostbed = host_bed_coordinates(
            output + "_viral.bed",
            output + "_source_seq_lengths.txt",
            fasta,
            minlength,
            output,
            keepBed=keep_bed,
        )
        os.remove(output + "_source_seq_lengths.txt")
    else:
        hostbed = None

    logger.newline()
    logger.info(
        f"Excluding contigs with contamination, longer than expected and duplication issues..."
    )

    # Return sets of contigs to include/exclude for further clustering
    include, exclude = lists_proviruses_warnings(qsummary, viralbed)

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
            if row["source_seq"] in include:
                include.remove(row["source_seq"])
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


def main():
    # Parse script arguments
    args = parse_arguments()
    fasta = args["fasta"]
    minlength = args["length"]
    keep_bed = args["keep_bed"]
    debug = args["debug"]
    output = args["output"]
    output_name = Path(args["output"]).name
    output_dir = Path(args["output"]).parent
    threads = args["threads"]

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Make tmp file for CheckV and geNomad outputs
    tmpdir = tempfile.mkdtemp(dir=output_dir)

    if args["checkv_db"] is not None and args["genomad_db"] is not None:
        # Run CheckV
        checkv_dict = {
            "fasta": fasta,
            "tmpdir": tmpdir,
            "checkv_db": Path(args["checkv_db"]),
            "threads": threads,
        }

        run_checkv(**checkv_dict)

        # Run genomad
        genomad_dict = {
            "fasta": fasta,
            "tmpdir": tmpdir,
            "genomad_db": Path(args["genomad_db"]),
            "threads": threads,
            "sens_marker": args["sens_marker"],
            "eval_marker": args["eval_marker"],
            "splits": args["splits"],
            "ct": args["ct"],
            "mt": args["mt"],
            "mti": args["mti"],
            "mte": args["mte"],
            "mid": args["mid"],
            "mtd": args["mtd"],
            "sens_integrase": args["sens_integrase"],
            "eval_integrase": args["eval_integrase"],
        }

        run_genomad(**genomad_dict)

        # Make fastas with sequences to cluster and sequences to reinclude in a later clustering stage (across your study)
        exclude_sequences_for_clustering(
            fasta, output, output_name, tmpdir, minlength, keep_bed, debug
        )

        logger.newline()

        # Cluster remaining sequences
        vu.clustering(
            os.path.join(tmpdir, output_name + "_to_cluster.fasta"),
            output,
            threads,
            args["pid"],
            args["cov"],
        )
    else:
        logger.info(f"Only clustering, no provirus identification.")
        vu.clustering(
            fasta,
            output,
            threads,
            args["pid"],
            args["cov"],
        )

    if not debug:
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    main()
