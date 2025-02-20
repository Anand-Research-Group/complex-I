import os
import numpy as np
import pandas as pd
import pytaxonkit
import requests
from tqdm import tqdm
from pathlib import Path
from config import PROKARYOTES_FILE, GENOME_DATASET_FILE, NCBI_GENOME_RECORDS_DIR

# NCBI FTP link for latest prokaryotes.txt
NCBI_PROKARYOTES_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"

def download_latest_prokaryotes():
    """Download the latest prokaryotes.txt file from NCBI if the local copy is missing or outdated."""
    prokaryotes_path = Path(PROKARYOTES_FILE)
    user_choice = input("Would you like to download the latest prokaryotes.txt from NCBI? (y/n): ").strip().lower()
    if user_choice == 'y':
        print(f"Downloading latest prokaryotes.txt from {NCBI_PROKARYOTES_URL}...")
        try:
            response = requests.get(NCBI_PROKARYOTES_URL, stream=True, timeout=60)
            response.raise_for_status()
            with open(prokaryotes_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"✅ Download complete: {prokaryotes_path}")
        except requests.RequestException as e:
            print(f"❌ Failed to download latest prokaryotes.txt: {e}")

def fetch_taxonomy(data, cpu=4):
    data["TaxID"] = data["TaxID"].astype(str)
    data_taxonkit = pytaxonkit.lineage(data["TaxID"].unique(), formatstr="{s}", threads=cpu)
    data_taxonkit["TaxID"] = data_taxonkit["TaxID"].astype(str)

    taxonomic_info = []
    for _, row in tqdm(data_taxonkit.iterrows(), total=data_taxonkit.shape[0], desc="Processing taxonomic information"):
        lineage_ranks = row['FullLineageRanks'].split(';')
        full_lineage = row['FullLineage'].split(';')
        taxonomy = dict(zip(lineage_ranks, full_lineage))
        taxonomic_info.append({
            'TaxID': row['TaxID'],
            'Organism': full_lineage[1] if len(full_lineage) > 1 else None,
            'Superkingdom': taxonomy.get('superkingdom'),
            'Phylum': taxonomy.get('phylum'),
            'Class': taxonomy.get('class'),
            'Order': taxonomy.get('order'),
            'Family': taxonomy.get('family'),
            'Genus': taxonomy.get('genus'),
            'Species': taxonomy.get('species'),
            'Strain': taxonomy.get('strain'),
        })
    return pd.DataFrame(taxonomic_info).drop_duplicates()

def load_prokaryotes():
    """Loads prokaryotes data and prepares it for processing."""
    prokaryotes_path = Path(PROKARYOTES_FILE)
    if not prokaryotes_path.exists():
        print(f"⚠️ {prokaryotes_path} not found.")
        download_latest_prokaryotes()
    prokaryotes = pd.read_csv(prokaryotes_path, sep='\t', low_memory=False)
    prokaryotes['FTP Path'] = prokaryotes['FTP Path'].replace('-', np.nan)
    prokaryotes.dropna(subset=['FTP Path'], inplace=True)
    ftp_to_path = lambda x: f"{x}/{x.split('/')[-1]}"
    prokaryotes["GenomePath"] = prokaryotes["FTP Path"].apply(lambda x: ftp_to_path(x) + "_genomic.fna.gz")
    prokaryotes["CDSPath"] = prokaryotes["FTP Path"].apply(lambda x: ftp_to_path(x) + "_cds_from_genomic.fna.gz")
    return prokaryotes

def save_output(genomes_records, dir_path):
    """Saves processed genome records to files."""
    dir_path = Path(dir_path)
    dir_path.mkdir(parents=True, exist_ok=True)
    output_files = {
        'CDSPath': dir_path / 'ncbi_cds_ftp.txt',
        'GenomePath': dir_path / 'ncbi_genomes_ftp.txt'
    }

    for column, filepath in output_files.items():
        if column in genomes_records.columns:
            genomes_records[[column]].to_csv(filepath, index=False, header=False)
        else:
            print(f"⚠️ Warning: '{column}' column is missing. Skipping file: {filepath}")

download_latest_prokaryotes()
data = load_prokaryotes()
data = data[data['Status'].isin(['Complete Genome', 'Chromosome', 'Complete', 'Chromosome(s)'])]
save_output(data, NCBI_GENOME_RECORDS_DIR)
data = data[['TaxID', 'Group', 'SubGroup', 'Size (Mb)', 'GC%', 'Genes', 'Proteins', 'Assembly Accession', 'Reference', 'FTP Path']]
taxonomy_data = fetch_taxonomy(data)
data['FTP Path'] = data['FTP Path'].replace('-', np.nan)
data.dropna(subset=['FTP Path'], inplace=True)
data = taxonomy_data[['Organism', 'Species', 'Strain', 'TaxID']].merge(data, on='TaxID', how='inner')
taxonomy_data.to_csv(NCBI_GENOME_RECORDS_DIR / 'taxonomy.csv', index=False)
data.to_csv(GENOME_DATASET_FILE, index=False)
