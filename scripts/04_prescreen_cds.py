import os
import re
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import warnings
from config import CDS_DIR, CDS_METADATA_DIR  # Import standardized paths

warnings.filterwarnings("ignore")

# Define output file paths
NUO_CDS_FILE = CDS_METADATA_DIR / "nuo_cds_prescreened.csv"
NDU_CDS_FILE = CDS_METADATA_DIR / "ndu_cds_prescreened.csv"

def extract_from_header(header, patterns):
    """Extracts values from a FASTA header using regex patterns."""
    return {key: (match.group(1) if (match := re.search(pattern, header)) else None) for key, pattern in patterns.items()}

def clean_gene_symbol(gene_symbol):
    """Formats gene symbols by removing special characters and standardizing casing."""
    return gene_symbol.lower().replace('[h', '').replace('[c', '').strip().replace('nuo', '').upper().translate(str.maketrans('', '', '-_/'))

def parse_cds_files(gene_initial="nuo"):
    """Parses CDS FASTA files to extract gene-specific records."""
    data = []
    patterns = {
        "accession": r"lcl\|(.*?)_cds",
        "gene": r"\[gene=(.*?)\]",
        "protein": r"\[protein=(.*?)\]"
    }

    for fasta in tqdm(os.listdir(CDS_DIR), desc=f"Processing {gene_initial.upper()} CDS files"):
        fasta_path = CDS_DIR / fasta
        for record in SeqIO.parse(fasta_path, "fasta"):
            if f"gene={gene_initial}" in record.description.lower():
                prot_seq = record.seq.translate(table=11, to_stop=True)
                record_info = extract_from_header(record.description, patterns)
                data.append([fasta, record.description] + list(record_info.values()) + [len(prot_seq)])

    df = pd.DataFrame(data, columns=['CDSFile', 'Header', 'Accession', 'GeneName', 'ProteinName', 'ProteinLength'])

    if not df.empty:
        df['GeneName'] = df['GeneName'].fillna('NONE').apply(clean_gene_symbol)
        df.drop_duplicates(subset=['CDSFile', 'GeneName'], inplace=True)
        df['Subunit'] = "Nuo" + df['GeneName'].str.replace(r'\d', '', regex=True)

        # Remove unwanted subunits
        unwanted_subunits = {"NuoBC", "Nuo", "NuoII", "NuoLM", "NuoP"}
        df = df[~df['Subunit'].isin(unwanted_subunits)]

        # Handle fused subunits
        fused_patterns = r"subunit C,D|subunit C/D|chain C,D|chain C/D|chain CD|NADH.*fused CD subunit"
        df.loc[df['Subunit'].eq('NuoC') & df['ProteinName'].str.contains(fused_patterns, case=False, na=False), 'Subunit'] = "NuoCD"

    return df

# Process both 'nuo' and 'nduf' genes
for gene, output_file in [("nuo", NUO_CDS_FILE), ("nduf", NDU_CDS_FILE)]:
    result_df = parse_cds_files(gene)
    if not result_df.empty:
        result_df.to_csv(output_file, index=False)
        print(f"✅ Saved {output_file}")
    else:
        print(f"⚠️ No data found for {gene}, skipping...")