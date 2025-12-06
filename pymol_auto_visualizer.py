import json
import pandas as pd
import os
import subprocess
import re
import shutil
import sys
import argparse
import urllib.request

# 1. Data Structure (The PDB Dictionary)
GENE_PDB_MAP = {
    'gyrA': '1KZN',
    'parC': '1Z4U',
    'blaCTX-M': '1YLJ',
    'tetA': '4TQU',
    'emrE': '2I68',
    'emrK': '3B5D',
    'emrA': 'AF-P27303', # AlphaFold UniProt ID for E. coli EmrA
    'emrB': 'AF-P0AEJ0'  # AlphaFold UniProt ID for E. coli EmrB
}

PDB_CACHE_DIR = "PDB_Cache"

def get_pdb_file(pdb_id):
    """
    Checks for PDB file in cache, downloads if missing.
    Returns the absolute path to the PDB file.
    Supports standard PDB IDs (RCSB) and AlphaFold IDs (via API).
    """
    if not os.path.exists(PDB_CACHE_DIR):
        os.makedirs(PDB_CACHE_DIR)
    
    # Determine filename and URL
    if pdb_id.startswith("AF-"):
        # Format: AF-UniProtID
        uniprot_id = pdb_id.split("-")[1]
        pdb_filename = f"{pdb_id}.pdb"
        
        # Fetch URL from AlphaFold API to get latest version
        api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        try:
            with urllib.request.urlopen(api_url) as response:
                data = json.loads(response.read().decode())
                # data is a list, take first entry
                if data and 'pdbUrl' in data[0]:
                    url = data[0]['pdbUrl']
                else:
                    print(f"Error: No PDB URL found in AlphaFold API for {uniprot_id}")
                    return None
        except Exception as e:
            print(f"Error fetching AlphaFold API for {uniprot_id}: {e}")
            return None
            
    else:
        pdb_filename = f"{pdb_id}.pdb"
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    pdb_path = os.path.join(PDB_CACHE_DIR, pdb_filename)
    
    if not os.path.exists(pdb_path):
        print(f"Downloading structure {pdb_id} from {url}...")
        try:
            urllib.request.urlretrieve(url, pdb_path)
            print(f"Downloaded {pdb_id} to {pdb_path}")
        except Exception as e:
            print(f"Error downloading {pdb_id}: {e}")
            return None
    
    return os.path.abspath(pdb_path)

def parse_mutations(file_path):
    """
    Parses the input file (CSV or Excel) and yields mutation details.
    """
    print(f"Loading data from {file_path}...")
    try:
        if file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path)
        else:
            df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    # Filter for variants
    if 'Variant' not in df.columns or 'Gene Name' not in df.columns:
        print("Error: Required columns 'Variant' or 'Gene Name' not found.")
        return

    filtered_df = df[df['Variant'] != "No variants found"]

    for _, row in filtered_df.iterrows():
        gene_raw = str(row['Gene Name'])
        variant_raw = str(row['Variant'])

        # Extract Gene Name
        gene_name = None
        for key in GENE_PDB_MAP:
            if key in gene_raw:
                gene_name = key
                break
        
        if not gene_name:
            parts = gene_raw.split('_')
            if parts and parts[0] in GENE_PDB_MAP:
                gene_name = parts[0]
            else:
                gene_name = parts[0] if parts else gene_raw

        # Extract Residue Number
        residue_match = re.search(r'(\d+)', variant_raw)
        residue_number = residue_match.group(1) if residue_match else None

        if gene_name and residue_number:
            pdb_id = GENE_PDB_MAP.get(gene_name)
            yield (gene_name, pdb_id, residue_number, variant_raw)

def validate_residue(pdb_path, residue):
    """
    Parses the PDB file to check if the residue number exists in any chain.
    Returns True if found, False otherwise.
    """
    residue_str = str(residue)
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # PDB format: Residue sequence number is columns 22-26
                    # But often it's just space separated. Let's rely on fixed width for standard PDB.
                    # Columns 22-26 (1-based) -> 22:26 in Python slice? No, 22 is 23rd char.
                    # Python slice: line[22:26] gives 4 chars. PDB spec says resSeq is 23-26.
                    # Let's try splitting, which is safer for non-standard files, 
                    # but PDB ATOM lines are strict.
                    # Standard: ATOM   1    N   ASP A   1      ...
                    # Split: ['ATOM', '1', 'N', 'ASP', 'A', '1', ...] -> index 5 is res num.
                    parts = line.split()
                    if len(parts) > 5:
                        # Check if residue matches. 
                        # Note: Insertion codes might complicate this, but simple check is enough.
                        if parts[5] == residue_str:
                            return True
                        # Sometimes chain is merged: "ASP A1" -> need to be careful.
                        # Fixed width is better: line[22:26].strip()
                        res_seq = line[22:26].strip()
                        if res_seq == residue_str:
                            return True
    except Exception as e:
        print(f"Warning: Could not validate residue in {pdb_path}: {e}")
        return True # Assume true if we can't parse
    return False

def create_pml_script(gene, pdb_path, residue, base_filename, interactive=False):
    """
    Writes a temporary .pml file with PyMOL commands.
    """
    pml_filename = f"{base_filename}.pml"
    png_filename = f"{base_filename}.png"
    
    # Use 'load' for local file instead of 'fetch'
    # Escape backslashes for PyMOL (Windows paths)
    pdb_path_safe = pdb_path.replace('\\', '/')
    
    pml_content = f"""
load {pdb_path_safe}
hide everything
show cartoon
color white
bg_color white
select mutation_site, resi {residue}
show spheres, mutation_site
color firebrick, mutation_site
zoom mutation_site, 15
png {png_filename}, width=1200, height=1200, ray=1
"""
    if not interactive:
        pml_content += "quit\n"
        
    with open(pml_filename, 'w') as f:
        f.write(pml_content)
    return pml_filename

def main():
    parser = argparse.ArgumentParser(description="Automate PyMOL visualizations from mutation CSV.")
    parser.add_argument("input_file", nargs='?', help="Path to input CSV/Excel file.")
    parser.add_argument("-i", "--interactive", action="store_true", help="Keep PyMOL open after loading (do not quit).")
    args = parser.parse_args()

    # Define input path
    input_dir = "Input"
    input_file = args.input_file
    
    if not input_file:
        # Auto-detect file in Input directory
        if os.path.exists(input_dir):
            files = os.listdir(input_dir)
            for f in files:
                if f.endswith('.xlsx') or f.endswith('.csv'):
                    input_file = os.path.join(input_dir, f)
                    break
        
        if not input_file:
            input_file = "AMR_Full_Analysis.xlsx - Sheet1.csv"

    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return

    output_dir = "PyMOL_Visuals"
    os.makedirs(output_dir, exist_ok=True)

    # Check for PyMOL
    pymol_exe = shutil.which("pymol")
    if not pymol_exe:
        print("Warning: 'pymol' executable not found in PATH.")
        print("Scripts will be generated but not executed.")

    print(f"Processing {input_file}...")
    
    count = 0
    for gene, pdb_id, residue, variant in parse_mutations(input_file):
        if not pdb_id:
            print(f"Skipping {gene}: No PDB ID mapped.")
            continue

        # Get PDB file (cached)
        pdb_path = get_pdb_file(pdb_id)
        if not pdb_path:
            continue

        # Validate Residue
        if not validate_residue(pdb_path, residue):
            print(f"Skipping {gene} {variant}: Residue {residue} not found in PDB {pdb_id}.")
            continue

        base_name = os.path.join(output_dir, f"{gene}_{variant}")
        pml_file = create_pml_script(gene, pdb_path, residue, base_name, interactive=args.interactive)
        
        # Run PyMOL
        if pymol_exe:
            print(f"> Generating visual for {gene} {variant}...")
            try:
                # If interactive, we generally don't want to block the loop or open 50 windows.
                # But for now, let's assume the user knows what they are doing if they use -i.
                # Or, if interactive, maybe we just generate the scripts and don't run them?
                # The user said "files generated are not opening".
                # So they are likely running the scripts manually.
                # So we just need to ensure the scripts don't have 'quit'.
                
                # If interactive, we skip execution to avoid spamming windows, unless user wants it.
                # Let's execute only if NOT interactive.
                if not args.interactive:
                    subprocess.run(["pymol", "-c", "-q", pml_file], check=True)
                    if os.path.exists(pml_file):
                        os.remove(pml_file)
                else:
                    print(f"> [Interactive] Created {pml_file} (Manual open required)")
            except subprocess.CalledProcessError as e:
                print(f"Error running PyMOL for {gene} {variant}: {e}")
        else:
            print(f"> [Dry Run] Created {pml_file} for {gene} {variant}")
        
        count += 1

    print(f"Done. Processed {count} mutations.")

if __name__ == "__main__":
    main()
