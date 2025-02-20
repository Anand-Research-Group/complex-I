import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from config import CDS_DIR, CDS_METADATA_DIR, HMM_CDS_SEQS_DIR  # Correct output directory

# Ensure output directory exists
HMM_CDS_SEQS_DIR.mkdir(parents=True, exist_ok=True)

def process_sequences(cds_df):
    """
    Process and save protein sequences from CDS data.

    :param cds_df: DataFrame containing CDS information
    """
    # Filter for chromosome-associated CDS entries
    chr_cds = cds_df[cds_df['Replicon'] == 'Chromosome']

    # Group CDS entries by file, mapping them to (Subunit, Header)
    cds_to_info = chr_cds.groupby('CDSFile')[['Subunit', 'Header']].apply(lambda x: x.values.tolist()).to_dict()
    
    # Initialize a dictionary to store sequences by subunit
    sequence_data = {subunit: [] for subunit in cds_df['Subunit'].unique()}

    # Process each CDS file
    for fasta in tqdm(chr_cds['CDSFile'].unique(), desc="Extracting sequences"):
        fasta_path = CDS_DIR / fasta

        if not fasta_path.exists():
            print(f"⚠️ Skipping missing file: {fasta_path}")
            continue

        for record in SeqIO.parse(fasta_path, "fasta"):
            for subunit, header in cds_to_info[fasta]:
                if header in record.description:
                    sequence_data[subunit].append(record.seq.translate(table=11, to_stop=True))

    # Write sequences to the correct HMM directory
    for subunit, sequences in sequence_data.items():
        records = [SeqRecord(seq, id=f"{subunit}_{i+1}", description=f"{subunit} Subunit") for i, seq in enumerate(sequences)]
        if records:
            output_file = HMM_CDS_SEQS_DIR / f"{subunit.lower()}_cds.faa"
            SeqIO.write(records, output_file, "fasta")
            print(f"✅ Saved sequences to: {output_file}")

# Load preprocessed CDS metadata
NUO_CDS_FILE = CDS_METADATA_DIR / "nuo_cds_prescreened.csv"
NDU_CDS_FILE = CDS_METADATA_DIR / "ndu_cds_prescreened.csv"

# Read input data
if NUO_CDS_FILE.exists():
    nuo_data = pd.read_csv(NUO_CDS_FILE)
    process_sequences(nuo_data)
else:
    print(f"⚠️ Missing file: {NUO_CDS_FILE}")

if NDU_CDS_FILE.exists():
    ndu_data = pd.read_csv(NDU_CDS_FILE)
    process_sequences(ndu_data)
else:
    print(f"⚠️ Missing file: {NDU_CDS_FILE}")
