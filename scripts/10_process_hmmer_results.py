import os
import re
import glob
import numpy as np
import pandas as pd
import logging
from tqdm import tqdm
from pathlib import Path
from config import HMM_RESULTS_DIR

# Setup logging
LOG_FILE = Path(__file__).parent / "hmmer_results.log"
logging.basicConfig(
    filename=LOG_FILE,
    filemode="w",
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logging.getLogger().addHandler(console_handler)

def parse_results_tblout_output(file_fullpath):
    """
    Parses an HMMER `.tblout` output file, extracting key information.

    Args:
        file_fullpath (str): Full path to the HMMER results file.

    Returns:
        pd.DataFrame: DataFrame containing parsed hit data.
    """
    accession_pattern = re.compile(r'_[0-9]+$')
    hits = []

    try:
        with open(file_fullpath, 'r') as handle:
            lines = handle.readlines()[3:]  # Skip header lines

        for line in lines:
            if not line.strip() or line.startswith('#'):
                continue  # Skip empty lines and comments

            cols = [x for x in line.strip().split() if x]

            hit = {
                'Accession': accession_pattern.sub('', cols[0]),
                'ProteinAccession': cols[0],
                'Profile': cols[2],
                'evalue': np.float64(cols[4]),
                'BitScore': np.float64(cols[5]),
                'Bias': np.float64(cols[6]),
                'SequenceDesc': " ".join(cols[18:])
            }

            hits.append(hit)

        if not hits:
            logging.warning(f"⚠️ No valid hits found in: {file_fullpath}")
            return pd.DataFrame()  # Return empty DataFrame if no hits

        return pd.DataFrame(hits)

    except Exception as e:
        logging.error(f"❌ Error processing file {file_fullpath}: {str(e)}")
        return pd.DataFrame()

def process_hmmer_results(result_dir, pattern_str=r'#\s*(\d+)\s*#\s*(\d+)\s*'):
    """
    Processes all HMMER `.tblout` results in a directory, cleaning and formatting them.

    Args:
        result_dir (str): Directory containing `.tblout` HMMER output files.
        pattern_str (str): Regex pattern for extracting 'Start' and 'End' from 'SequenceDesc'.

    Returns:
        pd.DataFrame: Processed results from all `.tblout` files.
    """
    result_dir = Path(result_dir)
    file_paths = list(result_dir.glob("**/*.txt"))  # Find all .txt results recursively

    if not file_paths:
        logging.error(f"❌ No result files found in {result_dir}")
        return pd.DataFrame()

    logging.info(f"📂 Found {len(file_paths)} result files in {result_dir}")

    all_dataframes = [parse_results_tblout_output(file) for file in tqdm(file_paths, desc="Processing HMMER results")]
    results = pd.concat([df for df in all_dataframes if not df.empty], ignore_index=True)

    if results.empty:
        logging.warning("⚠️ No valid results found after processing all files.")
        return results

    # Sort by 'evalue' and remove duplicates
    results.sort_values(by='evalue', inplace=True)
    results.drop_duplicates(subset=['Accession', 'ProteinAccession'], keep='first', inplace=True)

    # Split the 'Profile' column into three new columns: Subunit, SeqsClustThreshold, HMMParameter
    results[['Subunit', 'SeqsClustThreshold', 'HMMParameter']] = results['Profile'].str.split('_', expand=True)

    # Extract 'Start' and 'End' from 'SequenceDesc' using regex
    pattern = re.compile(pattern_str)
    matches = results['SequenceDesc'].apply(lambda x: pattern.search(x) if pd.notnull(x) else None)

    results['Start'] = matches.apply(lambda m: int(m.group(1)) if m else 0)
    results['End'] = matches.apply(lambda m: int(m.group(2)) if m else 0)

    # Apply log10 transformation to the 'evalue' column (avoid log(0) errors)
    results['log10evalue'] = results['evalue'].apply(lambda x: np.log10(x) if x > 0 else np.nan)

    # Remove the 'Profile' column
    results.drop(columns=['Profile', 'HMMParameter', 'SeqsClustThreshold'], inplace=True)

    # Reset index
    results.reset_index(drop=True, inplace=True)

    logging.info("✅ HMMER results processing complete.")
    return results

# **Execution**
if __name__ == "__main__":
    logging.info("🚀 Starting HMMER results processing...")
    processed_results = process_hmmer_results(HMM_RESULTS_DIR)

    if not processed_results.empty:
        output_file = HMM_RESULTS_DIR / "processed_hmmer_results.csv"
        processed_results.to_csv(output_file, index=False)
        logging.info(f"✅ Processed results saved to {output_file}")
    else:
        logging.warning("⚠️ No results were processed successfully.")

    print("✅ HMMER results processing complete! Logs saved to hmmer_results.log")
