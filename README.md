Small guide of how to use this to download genomes from NCBI.
1. Get a list of assembly accessions (for instance using the metadata of GTDB). They start with GCF (RefSeq) or GCA (GenBank). The file should be a single column csv with column name "assembly_accession". You can use `bac120_metadata_r220.tsv.gz` to filter the genomes according to what you need (release 220 is here : https://data.gtdb.ecogenomic.org/releases/release220/).
2. Download the necessary files :
    1. `ranked_lineage.dmp`obtained from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/ (decompress `new_taxdump.tar.gz`) (it is the only way to get the taxonmy of NCBI)
    2. `bac120_taxonomy_r220.tsv.gz` the taxonomy of GTDB
3. Get an API key from NCBI : https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/rest-api/authentication/
4. Export the api key as an environment variable called NCBI_API_KEY `export NCBI_API_KEY="yourapikey"`
5. Create the necessary conda environment `conda env create -f ncbi_datasets_env.yaml` and activate it.
6. Run the program : `download_api.py [-h] accession_csv rankedlineage gtdb_taxonomy genome_dir taxon_csv`. genome dir and taxon csv are (necessary) outputs.
⚠ : the download is processed by chunks of 1000 genomes (~ 5 Gb), you need at least this much free RAM (or edit lines 229 and 230).
