# NuoHMMER Pipeline

## Repository Overview
This repository contains a set of scripts and workflows designed to search for **Respiratory Complex I (NADH Ubiquinone Oxidoreductase) subunits** in prokaryotic genomes and proteomes. The pipeline automates:
- Downloading genome and coding sequence (CDS) data.
- Pre-screening for annotated subunits.
- Merging sequences with InterPro data.
- Clustering and performing multiple sequence alignments (MSA).
- Building and utilizing HMM profiles for protein searches.
- Filtering and analyzing results, including statistical assessments of oxygen tolerance relationships.

## Directory Structure
Upon execution of `01_setup_dirs.py`, the following directories will be created:

```
project_root/
│-- config.py  # Configuration file with project paths
│-- data/
│   ├── sequence_data/  # Stores genomic, CDS, and proteome data
│   │   ├── genomes/
│   │   ├── cds/
│   │   ├── proteomes/
│   ├── genomic_metadata/  # Metadata related to genomes
│   │   ├── genome_metadata/
│   │   ├── ncbi_genome_records/
│   ├── cds_metadata/  # Metadata specific to CDS
│   ├── external_metadata/  # Additional metadata from external sources
│   ├── hmm_data/  # Contains HMM profiles and related sequence data
│   │   ├── profiles/
│   │   ├── msa/
│   │   ├── interpro_prot_seqs/
│   │   ├── cds_prot_seqs/
│   │   ├── clustered_prot_seqs/
│   │   ├── clustered_msa_seqs/
│   │   ├── combined_interpro_cds_seqs/
```

## Scripts

### 1. **Setup and Configuration**
- **`01_setup_dirs.py`**: Creates the necessary directory structure and generates a `config.py` file with essential paths.

### 2. **Fetching Taxonomy and Preparing Downloads**
- **`02_fetch_taxonomy_prepare_downloads.py`**: Retrieves taxonomy information and generates FTP links for genome and CDS downloads.

### 3. **Downloading Genomes and CDS**
- **`03_download_genomes_cds.sh`**: Bash script that downloads genome and CDS files using generated FTP links.

### 4. **Pre-screening and Metadata Extraction**
- **`04_prescreen_cds.py`**: Screens CDS files for annotated Complex I subunits.
- **`05_extract_genome_metadata.py`**: Extracts genome metadata for further processing.
- **`06_extract_seqs_cds.py`**: Extracts sequences from CDS files for later HMM profiling.

### 5. **InterPro Data Integration**
- **`07_fetch_interpro_seqs.py`**: Retrieves InterPro sequences and integrates them with extracted CDS sequences.

### 6. **HMM-based Searches**
- **`08_hmm_pipeline.py`**: Constructs HMM profiles from MSA of clustered sequences.
- **`09_hmmer_search.py`**: Uses HMMER to search proteomes for Complex I subunits.
- **`10_process_hmmer_results.py`**: Process HMMER search results, combines into dataframe and saves them into csv format.

### 7. **Post-Processing and Analysis**
- **`post_search_01.ipynb`**: SAME as '10_process_hmmer_results.py'.
- **`post_search_02.ipynb`**: Intergenic distances and hits cluster analysis.
- **`post_search_03.ipynb`**: KDE evalues and hits cluster analysis.

## How to Run the Pipeline
1. **Set Up Project Directories**:
   ```bash
   python 01_setup_dirs.py
   ```
2. **Fetch Taxonomy and Prepare FTP Links**:
   ```bash
   python 02_fetch_taxonomy_prepare_downloads.py
   ```
3. **Download Genome and CDS Files**:
   ```bash
   bash 03_download_genomes_cds.sh
   ```
4. **Pre-screen CDS and Extract Metadata**:
   ```bash
   python 04_prescreen_cds.py
   python 05_extract_genome_metadata.py
   python 06_extract_seqs_cds.py
   ```
5. **Integrate InterPro Data and Build HMM Models**:
   ```bash
   python 07_fetch_interpro_seqs.py
   python 08_hmm_pipeline.py
   ```
6. **Run HMMER Search and Process Results**:
   ```bash
   python 09_hmmer_search.py
   python 10_process_hmmer_results.py
   ```
7. **Perform Post-Processing Analysis**:
   Open and run Jupyter notebooks `post_search_01.ipynb` and `post_search_02.ipynb` for visualization and statistical assessments.

## Requirements
- Python (≥3.8)
- Biopython
- HMMER
- MAFFT
- MMSeqs2
- Fasttree
- iqtree
- wget (for bash scripts)
- Jupyter Notebook (for post-processing analysis)

## Notes
- Users can modify scripts to search for different protein families by adjusting HMM profiles and filtering criteria.
- The `config.py` file centralizes all path references, ensuring easy modification if the project structure changes.

## Contact
For questions or issues, please reach out via GitHub Issues.

---
