import os
import shutil
import subprocess
import logging
from pathlib import Path
from tqdm import tqdm
from config import (
    HMM_CDS_SEQS_DIR, HMM_INTERPRO_SEQS_DIR, HMM_COMBINED_SEQS_DIR, HMM_MSA_SEQS_DIR,
    HMM_CLUST_SEQS_DIR, HMM_PROFILES_DIR
)

# Setup logging
LOG_FILE = Path(__file__).parent / "hmm_pipeline.log"
logging.basicConfig(
    filename=LOG_FILE,
    filemode="w",
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logging.getLogger().addHandler(console_handler)

def setup_directories(directories):
    """Creates the required directories if they do not exist."""
    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)
        logging.info(f"✅ Directory ensured: {directory}")

def run_command(command):
    """Runs a shell command, logs output and errors."""
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info(f"✅ Command succeeded: {' '.join(command)}")
        return result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"⚠️ Error running command: {' '.join(command)}\n{e.stderr}")
        return None

def gather_sequences(directory):
    """Gathers sequence files and associates them with their respective subunits."""
    sequences = [(fasta.split('_')[0].lower(), fasta) for fasta in os.listdir(directory) if fasta.endswith(".faa")]
    logging.info(f"✅ Gathered {len(sequences)} sequences from {directory}")
    return sequences

def combining_sequences(sequences, combined_dir, source_directory):
    """Copies sequence files to a combined directory organized by subunits."""
    combined_dir = Path(combined_dir)
    combined_dir.mkdir(parents=True, exist_ok=True)
    for subunit, fasta in sequences:
        destination_directory = combined_dir / subunit.upper()
        destination_directory.mkdir(exist_ok=True)
        source = Path(source_directory) / fasta
        destination = destination_directory / fasta
        shutil.copy(source, destination)
    logging.info(f"✅ Sequences copied from {source_directory} to {combined_dir}")

def concatenate_sequences(combined_dir, seq_dir):
    """Concatenates sequence files from combined directories into single files per subunit."""
    combined_dir, seq_dir = Path(combined_dir), Path(seq_dir)
    for subunit_dir in tqdm(combined_dir.iterdir(), desc="Concatenating sequences"):
        if subunit_dir.is_dir():
            concat_file = seq_dir / f"combined_cds_interpro_{subunit_dir.name.lower()}.faa"
            with concat_file.open('wb') as outfile:
                for fasta in subunit_dir.glob("*.faa"):
                    with fasta.open('rb') as infile:
                        shutil.copyfileobj(infile, outfile)
            logging.info(f"✅ Concatenated sequences for subunit: {subunit_dir.name}")

def run_mmseqs_commands(fasta_file, basename, output_dir, threshold=0.85):
    """Runs MMSeqs2 clustering commands on the provided fasta file."""
    output_dir = Path(output_dir)
    temp_dir = output_dir / "temp"
    temp_dir.mkdir(exist_ok=True)
    db_name, cluster_db, subset_db = temp_dir / f"{basename}_db", f"{db_name}_clu", f"{cluster_db}_rep"
    output_fasta = output_dir / f"{basename}_clustered_mmseq_{int(100 * threshold)}.fasta"

    commands = [
        ["mmseqs", "createdb", str(fasta_file), str(db_name)],
        ["mmseqs", "cluster", str(db_name), str(cluster_db), str(temp_dir), "--min-seq-id", str(threshold)],
        ["mmseqs", "createsubdb", str(cluster_db), str(db_name), str(subset_db)],
        ["mmseqs", "convert2fasta", str(subset_db), str(output_fasta)]
    ]

    for cmd in commands:
        run_command(cmd)

def run_mafft(msa_input_dir, msa_output_dir):
    """Performs multiple sequence alignment using MAFFT."""
    msa_sequences = sorted(Path(msa_input_dir).glob("*.fasta"))
    setup_directories([msa_output_dir])

    for seq in tqdm(msa_sequences, desc="Running MAFFT"):
        input_file, output_file = seq, Path(msa_output_dir) / seq.name.replace(".fasta", "_mafft_msa.fasta")
        mafft_command = f"mafft --localpair --maxiterate 400 --quiet --thread 8 {input_file} > {output_file}"
        result = run_command(mafft_command.split())
        if result:
            logging.info(f"✅ MAFFT alignment saved: {output_file}")

def run_hmmbuild(msa_input_dir, profile_dir):
    """Generates HMM profiles using HMMER's hmmbuild with optimized parameters."""
    msa_files = sorted(Path(msa_input_dir).glob("*.fasta"))
    setup_directories([profile_dir])

    for msa_file in tqdm(msa_files, desc="Building HMM profiles"):
        profile_name, profile_file_path = msa_file.stem, Path(profile_dir) / f"{msa_file.stem}.hmm"

        hmmbuild_command = [
            "hmmbuild",
            "--amino",
            "--cpu", "8",
            "-n", profile_name,
            "--wnone",
            "--symfrac", "0.6",
            "--fragthresh", "0.3",
            "--plaplace",
            str(profile_file_path),
            str(msa_file)
        ]

        result = run_command(hmmbuild_command)
        if result:
            logging.info(f"✅ HMM profile generated: {profile_file_path}")

# **Workflow Execution**
if __name__ == "__main__":
    logging.info("🚀 HMM Pipeline Execution Started")

    setup_directories([HMM_COMBINED_SEQS_DIR, HMM_MSA_SEQS_DIR, HMM_CLUST_SEQS_DIR, HMM_PROFILES_DIR])

    interpro_sequences, cds_sequences = gather_sequences(HMM_INTERPRO_SEQS_DIR), gather_sequences(HMM_CDS_SEQS_DIR)
    combining_sequences(interpro_sequences, HMM_COMBINED_SEQS_DIR, HMM_INTERPRO_SEQS_DIR)
    combining_sequences(cds_sequences, HMM_COMBINED_SEQS_DIR, HMM_CDS_SEQS_DIR)

    concatenate_sequences(HMM_COMBINED_SEQS_DIR, HMM_CLUST_SEQS_DIR)

    for fasta_file in tqdm(sorted(Path(HMM_CLUST_SEQS_DIR).glob("*.faa")), desc="Running MMSeqs2 clustering"):
        run_mmseqs_commands(fasta_file, fasta_file.stem, HMM_CLUST_SEQS_DIR, threshold=0.85)

    run_mafft(HMM_CLUST_SEQS_DIR, HMM_MSA_SEQS_DIR)
    run_hmmbuild(HMM_MSA_SEQS_DIR, HMM_PROFILES_DIR)

    logging.info("✅ HMM Pipeline Execution Complete")
    print("✅ HMM pipeline execution complete! Logs saved to hmm_pipeline.log")
