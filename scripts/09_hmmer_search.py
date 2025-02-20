import os
import psutil
import time
import subprocess
import platform
import logging
import argparse
from pathlib import Path
from tqdm import tqdm
from config import HMM_PROFILES_DIR, HMM_PROTEOMES_DIR, HMM_RESULTS_DIR

# Setup logging
LOG_FILE = Path(__file__).parent / "hmmer_search.log"
logging.basicConfig(
    filename=LOG_FILE,
    filemode="w",
    format="%(asctime)s - %(levelname)s - %(message)s",
    level=logging.INFO
)
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logging.getLogger().addHandler(console_handler)

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Run HMMER searches with power-aware execution.")
parser.add_argument("--force-run", action="store_true", help="Bypass laptop shutdown on battery")
args = parser.parse_args()

def detect_system_type():
    """Determines if the system is a laptop or a desktop."""
    try:
        if psutil.sensors_battery():
            return "laptop"
    except AttributeError:
        pass  # Some systems may not support battery info
    return "desktop"

def is_plugged_in():
    """Checks if the laptop is plugged in."""
    battery = psutil.sensors_battery()
    return battery is not None and battery.power_plugged

def shutdown_mac():
    """Shuts down a Mac system."""
    logging.warning("⚠️ System running on battery. Initiating shutdown (Mac).")
    os.system("osascript -e 'tell app \"System Events\" to shut down'")

def shutdown_linux():
    """Shuts down a Linux system."""
    logging.warning("⚠️ System running on battery. Initiating shutdown (Linux).")
    os.system("shutdown -h now")

def run_command(command):
    """Runs a shell command, logs output and errors."""
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info(f"✅ Command succeeded: {' '.join(command)}")
        return result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Command failed: {' '.join(command)}\n{e.stderr}")
        return None

# Detect system type
system_type = detect_system_type()
num_cpus = psutil.cpu_count(logical=True)

# Adjust CPU allocation
if system_type == "laptop":
    cpu_allocation = min(4, num_cpus)  # Use at most 4 CPUs on laptops
else:
    cpu_allocation = num_cpus  # Use all CPUs on desktops

logging.info(f"🖥️  Detected system: {system_type.upper()}")
logging.info(f"🔢 Allocating {cpu_allocation} CPUs for HMMER")

# **Laptop Power Management: Shutdown unless `--force-run` is used**
if system_type == "laptop" and not is_plugged_in():
    if args.force_run:
        logging.warning("⚠️ System is on battery, but `--force-run` is enabled. Continuing execution.")
    else:
        logging.critical("⚠️ System is running on battery. Shutting down to prevent power loss.")
        if platform.system() == "Darwin":
            shutdown_mac()
        elif platform.system() == "Linux":
            shutdown_linux()
        exit(1)  # Stop further execution

# Ensure required directories exist
HMM_PROFILES_DIR.mkdir(parents=True, exist_ok=True)
HMM_RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# List all proteome files
proteome_files = sorted(HMM_PROTEOMES_DIR.glob("*.faa"))
logging.info(f"📂 Found {len(proteome_files)} proteome files to process.")

# Run HMMER search for each profile
if HMM_PROFILES_DIR.is_dir():
    for profile_file in HMM_PROFILES_DIR.glob("*.hmm"):
        profile_results_dir = HMM_RESULTS_DIR / HMM_PROFILES_DIR.name / profile_file.stem
        profile_results_dir.mkdir(parents=True, exist_ok=True)

        logging.info(f"🔍 Processing HMM profile: {profile_file.name}")

        # Run HMMER for each proteome file
        for proteome_file in tqdm(proteome_files, desc=f"{profile_file.name} Search"):
            result_filename = f"{proteome_file.stem}_results.txt"
            result_file_path = profile_results_dir / result_filename

            hmmer_command = [
                "hmmsearch",
                "--cpu", str(cpu_allocation),
                "--noali",
                "--tblout", str(result_file_path),
                str(profile_file),
                str(proteome_file)
            ]

            result = run_command(hmmer_command)
            if result:
                logging.info(f"✅ HMMER search completed: {profile_file.name} → {proteome_file.name}")

        # Cooling period after processing each profile
        logging.info(f"🛑 Cooling period after processing {profile_file.name}. Waiting for 5 minutes...")
        time.sleep(300)  # 5 minutes pause

logging.info("✅ HMMER search completed successfully!")
print("✅ HMMER search completed! Logs saved to hmmer_search.log")
