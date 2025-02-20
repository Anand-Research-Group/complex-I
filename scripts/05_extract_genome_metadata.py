import os
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from config import GENOMES_DIR, GENOME_METADATA_DIR  # Import standardized paths

# Define output file
GENOME_METADATA_FILE = GENOME_METADATA_DIR / "genomes_metadata.csv"

def extract_genome_metadata():
    """Extracts genome metadata from FASTA headers."""
    metadata = []

    for genome in tqdm(os.listdir(GENOMES_DIR), desc="Extracting genome metadata"):
        fasta_path = GENOMES_DIR / genome

        for record in SeqIO.parse(fasta_path, "fasta"):
            accession = record.id
            description = record.description
            seq_length = len(record.seq) / 1_000_000  # Convert to megabases

            if 'plasmid' in description.lower():
                replicon_type = 'Plasmid'
            elif 'chromosome' in description.lower() or 'genome' in description.lower():
                replicon_type = 'Chromosome'
            else:
                replicon_type = 'Undefined'

            metadata.append((accession, replicon_type, genome, seq_length))

    return pd.DataFrame(metadata, columns=['Accession', 'Replicon', 'GenomeFile', 'SequenceLength(Mb)'])

# Process and save genome metadata
genome_df = extract_genome_metadata()
if not genome_df.empty:
    genome_df.to_csv(GENOME_METADATA_FILE, index=False)
    print(f"✅ Saved genome metadata to {GENOME_METADATA_FILE}")
else:
    print("⚠️ No genome metadata extracted.")

