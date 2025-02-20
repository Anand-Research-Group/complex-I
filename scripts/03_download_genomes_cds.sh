#!/bin/bash

# Define directories
BASE_DIR=$(pwd)  # Adjust if needed
GENOMES_DIR="$BASE_DIR/data/sequence_data/genomes"
CDS_DIR="$BASE_DIR/data/sequence_data/cds"
GENOMES_FTP_FILE="$BASE_DIR/data/genomic_metadata/ncbi_genomes_ftp.txt"
CDS_FTP_FILE="$BASE_DIR/data/genomic_metadata/ncbi_cds_ftp.txt"

# Create directories if they don't exist
mkdir -p "$GENOMES_DIR" "$CDS_DIR"

# Function to download files
download_files() {
    local ftp_file=$1
    local target_dir=$2
    
    if [ -f "$ftp_file" ]; then
        echo "Downloading files from $ftp_file to $target_dir"
        wget -i "$ftp_file" -P "$target_dir" --no-clobber --continue
    else
        echo "⚠️ File $ftp_file not found, skipping..."
    fi
}

# Download genome and CDS files
download_files "$GENOMES_FTP_FILE" "$GENOMES_DIR"
download_files "$CDS_FTP_FILE" "$CDS_DIR"

echo "🎯 All genome and CDS downloads are complete!"
