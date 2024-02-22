#!/usr/bin/env python3
import argparse
import os
import subprocess as sp
import shutil
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

    extr_ass = pd.read_csv(accession_file)
    list_accessions = list(extr_ass["assembly_accession"])
    all_accessions = "%2C".join(list_accessions)

    full_url = api_server + "/genome/accession/" + all_accessions + "/download"
    params = {"api-key": os.getenv("NCBI_API_KEY"),
              "include_annotation_type": "GENOME_FASTA"}
    res = requests.get(full_url, params=params)
    with open(zip_outfile, "wb") as fileout:
        for chunk in res.iter_content(chunk_size=255):
              if chunk:
                fileout.write(chunk)

    extract_zip(zip_outfile, genome_dir)
    taxon_csv = create_taxon_csv(genome_dir, "./assembly_data_report.jsonl", ["genus", "family", "species"], args.rankedlineage)
    taxon_csv.to_csv(out_csv)


if __name__ == "__main__":
    main()
