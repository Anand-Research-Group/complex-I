import os
from pathlib import Path

# Prompt user for the base directory
base_dir = Path(input("Enter the base directory for the project: ")).resolve()

# Define required subdirectories
subdirs = {
    "sequence_data": base_dir / "data/sequence_data",
    "genomes": base_dir / "data/sequence_data/genomes",
    "cds": base_dir / "data/sequence_data/cds",
    "proteomes": base_dir / "data/sequence_data/proteomes",
    
    # Metadata directories
    "genomic_metadata": base_dir / "data/genomic_metadata",
    "genome_metadata": base_dir / "data/genomic_metadata/genome_metadata",
    "ncbi_genome_records": base_dir / "data/genomic_metadata/ncbi_genome_records",
    "cds_metadata": base_dir / "data/cds_metadata",
    
    # External metadata
    "external_metadata": base_dir / "data/external_metadata",

    # HMM directories
    "hmm_analysis": base_dir / "data/hmm_data",
    "hmm_profiles": base_dir / "data/hmm_data/profiles",
    "hmm_interpro_seqs": base_dir / "data/hmm_data/interpro_prot_seqs",
    "hmm_cds_seqs": base_dir / "data/hmm_data/cds_prot_seqs",
    "hmm_clust_seqs": base_dir / "data/hmm_data/clustered_prot_seqs",
    "hmm_msa_seqs": base_dir / "data/hmm_data/clustered_msa_seqs",
    "hmm_combined_seqs": base_dir / "data/hmm_data/combined_interpro_cds_seqs",
}

# Create directories if they don't exist
for key, path in subdirs.items():
    path.mkdir(parents=True, exist_ok=True)
    print(f"✅ Created: {path}")

# Generate `config.py`
config_content = f"""from pathlib import Path

BASE_DIR = Path("{base_dir}")

# Data directories
SEQUENCE_DATA_DIR = BASE_DIR / "data/sequence_data"
GENOMES_DIR = SEQUENCE_DATA_DIR / "genomes"
CDS_DIR = SEQUENCE_DATA_DIR / "cds"
PROTEOMES_DIR = SEQUENCE_DATA_DIR / "proteomes"

# Metadata directories
GENOMIC_METADATA_DIR = BASE_DIR / "data/genomic_metadata"
GENOME_METADATA_DIR = GENOMIC_METADATA_DIR / "genome_metadata"
NCBI_GENOME_RECORDS_DIR = GENOMIC_METADATA_DIR / "ncbi_genome_records"
CDS_METADATA_DIR = BASE_DIR / "data/cds_metadata"

# External metadata directory (NEW)
EXTERNAL_METADATA_DIR = BASE_DIR / "data/external_metadata"

# HMM analysis
HMM_ANALYSIS_DIR = BASE_DIR / "data/hmm_data"
HMM_PROFILES_DIR = HMM_ANALYSIS_DIR / "profiles"
HMM_INTERPRO_SEQS_DIR = HMM_ANALYSIS_DIR / "interpro_prot_seqs"
HMM_CDS_SEQS_DIR = HMM_ANALYSIS_DIR / "cds_prot_seqs"
HMM_CLUST_SEQS_DIR = HMM_ANALYSIS_DIR / "clustered_prot_seqs"
HMM_MSA_SEQS_DIR = HMM_ANALYSIS_DIR / "clustered_msa_seqs"
HMM_COMBINED_SEQS_DIR = HMM_ANALYSIS_DIR / "combined_interpro_cds_seqs"

# Output directories
OUTPUT_DIR = SEQUENCE_DATA_DIR / "clustered_protein_sequences"

# Specific file paths
PROKARYOTES_FILE = NCBI_GENOME_RECORDS_DIR / "prokaryotes.txt"
GENOME_DATASET_FILE = GENOME_METADATA_DIR / "genomes_dataset.csv"
NUO_CDS_FILE = CDS_METADATA_DIR / "nuo_cds_prescreened.csv"
NDU_CDS_FILE = CDS_METADATA_DIR / "ndu_cds_prescreened.csv"
GENOME_METADATA_FILE = GENOME_METADATA_DIR / "genomes_metadata.csv"

# InterPro classification file (NEW)
INTERPRO_CSV = EXTERNAL_METADATA_DIR / "nuo_interpro_classification_accessions.csv"
"""

# Write config.py
config_path = base_dir / "config.py"
with open(config_path, "w", encoding="utf-8") as f:
    f.write(config_content)

print(f"\n✅ Config file created: {config_path}")
