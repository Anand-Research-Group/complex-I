Project-Repo/
│-- setup/                      # Scripts for initial project setup
│   ├── 01_setup_dirs.py        # Creates directory structure & Generates config file
│
│-- scripts/                    # Core bioinformatics workflow scripts
│   ├── 02_fetch_taxonomy_prepare_downloads.py
│   ├── 03_download_genomes_cds.sh  # Fetches required genome/CDS data
│   ├── 04_prescreen_cds.py
│   ├── 05_extract_genome_metadata.py
│   ├── 06_extract_seqs_cds.py
│   ├── 07_fetch_interpro_seqs.py
│   ├── 08_hmm_pipeline.py
│   ├── 09_hmmer_search.py
│   ├── 10_process_hmmer_results.py
│
│-- notebooks/                  # Jupyter Notebooks for interactive analysis
│   ├── 01_data_overview.ipynb  # Exploratory data analysis
│   ├── 02_hmm_results.ipynb    # Visualization of HMM search results
│   ├── 03_phylogenetics.ipynb  # Phylogenetic tree analysis
│
│-- config.py                   # Configuration file with all directory paths
│-- requirements.txt            # List of dependencies
│-- README.md                   # General project overview
│-- INSTALLATION.md             # Full setup guide
│-- WORKFLOW.md                 # Explanation of analysis steps
│
│-- results/                    # Output files (optional placeholder)
│   ├── figures/                # Barplots, scatterplots, phylogenetic trees
│   ├── tables/                 # Summary tables for manuscript
