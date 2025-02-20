import os
import sys
import json
import ssl
from time import sleep
from urllib import request
from urllib.error import HTTPError
from tqdm import tqdm
import pandas as pd
from config import HMM_INTERPRO_SEQS_DIR, INTERPRO_CSV  # Use updated config paths

class InterProFetcher:
    HEADER_SEPARATOR = "|"
    LINE_LENGTH = 80

    def __init__(self, interpro_accession, ncbi_taxid, subunit, organism):
        self.interpro_accession = interpro_accession
        self.ncbi_taxid = ncbi_taxid
        self.subunit = subunit
        self.organism = organism
        self.base_url = f"https://www.ebi.ac.uk:443/interpro/api/protein/reviewed/entry/InterPro/{interpro_accession}/taxonomy/uniprot/{ncbi_taxid}/?page_size=200&extra_fields=sequence"
        self.output_file = HMM_INTERPRO_SEQS_DIR / f"{subunit.lower()}_{organism}_{interpro_accession.lower()}_interpro.faa"

    def fetch_and_save_sequences(self):
        """Fetches InterPro sequences and saves them in FASTA format."""
        next_url = self.base_url

        with open(self.output_file, "w") as file:
            while next_url:
                data = self.fetch_data(next_url)
                if not data:
                    break

                for item in data["results"]:
                    self.write_fasta(item, file)

                next_url = data.get("next")
                sleep(1)  

# Ensure output directory exists
HMM_INTERPRO_SEQS_DIR.mkdir(parents=True, exist_ok=True)

# Check if the InterPro file exists
if not INTERPRO_CSV.exists():
    sys.stderr.write(f"⚠️ Missing InterPro classification file: {INTERPRO_CSV}\n")
    sys.exit(1)

# Load InterPro accessions
interpro_data = pd.read_csv(INTERPRO_CSV)
print(f"✅ Loaded {len(interpro_data)} InterPro accessions from {INTERPRO_CSV}")

# Define taxonomic groups
prok_ncbi_taxid = [('archaea', 2157), ('bacteria', 2)]

# Fetch sequences
for organism, ncbi_taxid in prok_ncbi_taxid:
    for _, row in tqdm(interpro_data.iterrows(), total=len(interpro_data), desc=f"Fetching sequences for {organism}"):
        subunit_protein = row['Protein']
        interpro_accession = row['InterPro Accession']
        fetcher = InterProFetcher(interpro_accession, ncbi_taxid, subunit_protein, organism)
        fetcher.fetch_and_save_sequences()

print("✅ InterPro sequence fetching complete.")
