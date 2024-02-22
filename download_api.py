#!/usr/bin/env python3
import argparse
import json
import os
import subprocess as sp
import shutil
import sys
import tempfile
import zipfile

import pandas as pd
import requests


def create_taxon_csv(genome_dir, ncbi_jsonl, taxons=["genus"], ranked_lineage_path="rankedlineage.dmp"):
    """Creates a csv with the filename and the taxon"""
    ranked_lineage = pd.read_table(ranked_lineage_path, sep=r"\t\|\t")
    tsv_ncbi = "ncbi_tsv.tsv"
    with open(tsv_ncbi, "w") as tsv_out:
        sp.run([
            "dataformat",
            "tsv",
            "genome",
            "--fields",
            "organism-tax-id,current-accession",
            "--inputfile",
            ncbi_jsonl], stdout=tsv_out)

    ncbi_df = pd.read_csv(tsv_ncbi, sep="\t")
    ncbi_df.columns = ["tax_id", "assembly_accession"]
    ranked_lineage.columns = [ "tax_id",
                               "tax_name",
                               "species",
                               "genus",
                               "family",
                               "order",
                               "class",
                               "phylum",
                               "kingdom",
                               "superkingdom"
                              ]
    assemblies_taxon = pd.merge(ncbi_df, ranked_lineage[taxons + ["tax_id"]], how="left", on="tax_id")
    genome_files = [g for g in os.listdir(genome_dir) if g.endswith((".fna", ".fasta", ".fa"))]
    genomes_df = pd.DataFrame(genome_files, columns=["genome"])
    genomes_df["assembly_accession"] = genomes_df["genome"].apply(lambda x: x[:15])
    res_df = pd.merge(genomes_df, assemblies_taxon, how="left", on="assembly_accession")
    return res_df[taxons + ["genome"]]


def extract_zip(path_zip, out_genome_dir):
    """Extracts and reorganises ncbi's api's zip output"""
    os.makedirs(out_genome_dir, exist_ok=True)

    zf = zipfile.ZipFile(path_zip)
    with tempfile.TemporaryDirectory() as tmpdirname:
        zf.extractall(tmpdirname)
        data_dir = os.path.join(tmpdirname, "ncbi_dataset/data/")
        for dirpath, _, filename in os.walk(data_dir):
            if filename[0].endswith((".fna", ".fasta", ".fa")):
                # copy genomes to genome directory
                shutil.copy(os.path.join(dirpath, filename[0]), out_genome_dir)
            else:
                # copy jsons to current wd
                for fn in filename:
                    shutil.copy(os.path.join(dirpath, fn), os.getcwd())


def SimpleFastaParser(handle):
    """Iterate over Fasta records as string tuples.

    Arguments:
     - handle - input stream opened in text mode

    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.

    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')

    adapted from biopython : https://github.com/biopython/biopython/blob/master/Bio/SeqIO/FastaIO.py
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")


def fasta_writer(tiselist, fasta_out):
    os.makedirs(os.path.dirname(fasta_out), exist_ok=True)
    with open(fasta_out, "w") as fileout:
        for title, seq in tiselist:
            lines = [f">{title}\n"]
            for i in range(0, len(seq), 60):
                lines.append(seq[i : i + 60] + "\n")
            fileout.write("".join(lines))


def write_only_chromosomes(refseq_dir_in, refseq_dir_out):
    """Filters out plasmids from a fasta file, write a new fasta to a new location."""
    fasta_list = [fa for fa in os.listdir(refseq_dir_in) if fa.endswith((".fasta", ".fa", ".fna"))]
    for fasta in fasta_list:
        tiselist = []
        with open(os.path.join(refseq_dir_in, fasta), "r") as filein:
            for title, seq in SimpleFastaParser(filein):
                if not "plasmid" in title:
                    tiselist.append((title, seq))
        os.makedirs(refseq_dir_out, exist_ok=True)
        fasta_writer(tiselist, os.path.join(refseq_dir_out, fasta))


def main():
    parser = argparse.ArgumentParser(
        description="""
        Download genomes from NCBI based on a csv containing assembly accessions.""")
    parser.add_argument(
        "--api_server",
        type=str,
        help="The base url of the ncbi api server",
        default= "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
    )
    parser.add_argument(
        "accession_csv",
        type=str,
        help="The csv with a column assembly_accession containing assemblies to download"
    )
    parser.add_argument(
        "rankedlineage",
        type=str,
        help="the rankedlineage.dmp file from ncbi taxonomy newtaxdump.gz"
    )
    parser.add_argument(
        "zip_outfile",
        type=str,
        help="zip output file path"
    )
    parser.add_argument(
        "genome_dir",
        type=str,
        help="output genomes_directory"
    )
    parser.add_argument(
        "taxon_csv",
        type=str,
        help="The output csv linkning genomes and taxon levels"
    )

    args = parser.parse_args()
    # expects a csv with a column assembly_accession
    accession_file = args.accession_csv
    zip_outfile = args.zip_outfile
    # base address of ncbi api server
    api_server = args.api_server
    genome_dir = args.genome_dir
    out_csv = args.taxon_csv

    accession_df = pd.read_csv(accession_file)
    list_accessions = list(accession_df["assembly_accession"])

    full_url = api_server + "/genome/download"
    payload = {"accessions": list_accessions, "include_annotation_type" : ["GENOME_FASTA"]}
    params = {"api-key": os.getenv("NCBI_API_KEY")}
    res = requests.post(full_url, params=params, json=payload)
    if res:
        with open(zip_outfile, "wb") as fileout:
            fileout.write(res.content)
        extract_zip(zip_outfile, genome_dir)
        taxon_csv = create_taxon_csv(genome_dir, "./assembly_data_report.jsonl", ["genus", "family", "species"], args.rankedlineage)
        taxon_csv.to_csv(out_csv)
    else:
        print(f"Request failed : {res.status_code}")
        sys.exit()

    genome_dir_no_plasmid = genome_dir + "_no_plasmid"
    write_only_chromosomes(genome_dir, genome_dir_no_plasmid)

if __name__ == "__main__":
    main()
