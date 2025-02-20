from pathlib import Path

BASE_DIR = Path("/Users/akshayonly/Work/Submission-Data/demo")

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
NUO_CDS_FILE = CDS_METADATA_DIR / "cds_subunits_metadata/nuo_cds_prescreened.csv"
