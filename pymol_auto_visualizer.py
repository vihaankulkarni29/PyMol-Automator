import json
import pandas as pd
import os
import subprocess
import re
import shutil
import sys
import argparse
import urllib.request
import time
from Bio.PDB import PDBParser, PPBuilder, PPBuilder, Select, PDBIO
from Bio.SeqUtils import seq1

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
ALIGNMENT_DIR = "Input/WildTypeAligner Output"
FASTA_DIR = "Input/FastaAAExtractor Output"

# Amino acid single letter to three letter code mapping
AA_CODE_MAP = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

def aa_single_to_three(single_code):
    """Convert single letter amino acid code to three letter code."""
    return AA_CODE_MAP.get(single_code.upper(), 'UNK')

def parse_alignment_file(alignment_path):
    """
    Parse alignment file to extract mutations.
    Format: Query line, marker line (| for match, space for mismatch), Reference line
    Returns dict with: wt_sequence, mut_sequence, mutations (list of tuples: position, wt_aa, mut_aa)
    """
    if not os.path.exists(alignment_path):
        return None
    
    try:
        with open(alignment_path, 'r') as f:
            lines = f.readlines()
        
        wt_seq = ""
        mut_seq = ""
        mutations = []
        
        # Parse the alignment blocks
        i = 0
        while i < len(lines):
            line = lines[i].rstrip('\n')
            parts = line.split()
            
            # Skip header and empty lines
            if not parts or line.startswith('#') or line.startswith('Needleman') or \
               line.startswith('Query:') or line.startswith('Reference:') or \
               line.startswith('Score:') or line.startswith('Identities:') or \
               line.startswith('Gaps:') or line.startswith(' '):
                i += 1
                continue
            
            # Check if this is a query sequence line (CP107114.1_emrA with digit position)
            if len(parts) >= 3 and parts[1].isdigit() and 'CP' in parts[0] or 'NZ_' in parts[0]:
                query_start = int(parts[1])
                query_seq = parts[2]
                
                # Next line should be markers (all spaces at start)
                # Line after that should be reference
                if i + 2 < len(lines):
                    ref_line = lines[i + 2].rstrip('\n')
                    ref_parts = ref_line.split()
                    
                    # Reference line: "emrA MG1655.faa 1 SEQUENCE 60"
                    # Need to find the numeric position and sequence
                    ref_seq = None
                    for j, part in enumerate(ref_parts):
                        if part.isdigit() and int(part) == query_start:
                            # Sequence should be the next part
                            if j + 1 < len(ref_parts):
                                ref_seq = ref_parts[j + 1]
                                break
                    
                    if ref_seq:
                        # Add sequences
                        mut_seq += query_seq
                        wt_seq += ref_seq
                        
                        # Find mutations by comparing sequences directly
                        for idx, (q_aa, r_aa) in enumerate(zip(query_seq, ref_seq)):
                            if q_aa != r_aa:
                                pos = query_start + idx
                                mutations.append((pos, r_aa, q_aa))
                        
                        i += 3
                        continue
            
            i += 1
        
        return {
            'wt_sequence': wt_seq,
            'mut_sequence': mut_seq,
            'mutations': mutations
        }
    except Exception as e:
        print(f"Error parsing alignment file {alignment_path}: {e}")
        return None

def read_fasta_sequence(fasta_path):
    """
    Read a FASTA file and return the sequence (without header).
    """
    if not os.path.exists(fasta_path):
        return None
    
    try:
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
        
        sequence = ""
        for line in lines:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line
        
        return sequence
    except Exception as e:
        print(f"Error reading FASTA {fasta_path}: {e}")
        return None

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
    Parses the input file (CSV or Excel) and yields mutation details with alignment data.
    Returns: gene_name, pdb_id, genome_id, mutations_list, wt_sequence, mut_sequence
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
    
    # Group by genome_id and gene to avoid processing same combination multiple times
    processed = set()

    for _, row in filtered_df.iterrows():
        gene_raw = str(row['Gene Name'])
        variant_raw = str(row['Variant'])

        # Extract Gene Name and Genome ID
        # Format: emrA_CP107114.1_emrA or gyrA_CP107114.1_gyrA
        gene_name = None
        genome_id = None
        
        parts = gene_raw.split('_')
        
        # Look for genome accession (contains CP, NZ_CP, etc. with numbers)
        for part in parts:
            if re.match(r'(NZ_)?[A-Z]+\d+', part):  # Matches CP107114.1, NZ_CP107117, etc.
                genome_id = part
                break
        
        # Extract gene name - check if it's in our PDB map
        for key in GENE_PDB_MAP:
            if key in gene_raw.lower():
                gene_name = key
                break
        
        if not gene_name:
            if parts and parts[0] in GENE_PDB_MAP:
                gene_name = parts[0]
            else:
                # Skip genes not in PDB map
                continue

        # Skip if no genome ID found or already processed
        if not genome_id:
            continue
            
        combo_key = f"{genome_id}_{gene_name}"
        if combo_key in processed:
            continue
        processed.add(combo_key)
        
        # Get PDB ID
        pdb_id = GENE_PDB_MAP.get(gene_name)
        if not pdb_id:
            continue
        
        # Look for alignment file
        alignment_filename = f"{gene_name}_{genome_id}_{gene_name}.txt"
        alignment_path = os.path.join(ALIGNMENT_DIR, alignment_filename)
        
        if not os.path.exists(alignment_path):
            print(f"Warning: Alignment file not found: {alignment_path}")
            continue
        
        # Parse alignment to get mutations
        alignment_data = parse_alignment_file(alignment_path)
        if not alignment_data or not alignment_data['mutations']:
            print(f"Info: No mutations found in alignment for {gene_name} {genome_id}")
            continue
        
        # Get mutant sequence from FASTA file
        fasta_filename = f"{genome_id}_{gene_name}.faa"
        fasta_path = os.path.join(FASTA_DIR, fasta_filename)
        
        mut_sequence = None
        if os.path.exists(fasta_path):
            mut_sequence = read_fasta_sequence(fasta_path)
        
        if not mut_sequence:
            mut_sequence = alignment_data['mut_sequence']
        
        yield {
            'gene_name': gene_name,
            'pdb_id': pdb_id,
            'genome_id': genome_id,
            'mutations': alignment_data['mutations'],  # List of (pos, wt_aa, mut_aa)
            'wt_sequence': alignment_data['wt_sequence'],
            'mut_sequence': mut_sequence
        }

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

def create_mutant_pdb_by_threading(pdb_path, mutations, output_pdb):
    """
    Create a mutant PDB by threading the mutated sequence onto the structure.
    Uses PyMOL to substitute residues and adds REMARK lines to highlight mutations.
    mutations: list of (position, wt_aa, mut_aa) tuples
    """
    # Create PML script to perform mutations using PyMOL's residue editing
    temp_pml = output_pdb.replace('.pdb', '_thread.pml')
    pdb_path_safe = pdb_path.replace('\\', '/')
    output_pdb_safe = output_pdb.replace('\\', '/')
    
    # Create mutation annotation string
    mutation_info = "; ".join([f"{wt}{pos}{mut}" for pos, wt, mut in mutations])
    
    pml_content = f"""# Mutation threading script
load {pdb_path_safe}, protein

# Select and show the structure
hide everything
show cartoon, protein
show sticks, protein

"""
    
    # Generate PyMOL commands to mutate each position
    for pos, wt_aa, mut_aa in mutations:
        mut_aa_three = aa_single_to_three(mut_aa)
        pml_content += f"""# Mutate position {pos}: {wt_aa} -> {mut_aa}
alter protein and resi {pos}, resn='{mut_aa_three}'

"""
    
    # Create temp output with mutation remarks
    pml_content += f"""# Rebuild the structure and save
rebuild

# Save to temporary file first
save {output_pdb_safe}.tmp, protein
quit
"""
    
    with open(temp_pml, 'w') as f:
        f.write(pml_content)
    
    # We'll add remarks after PyMOL creates the file
    # Store mutation info in a way that post_process_pdb can access it
    temp_pml_mutations = temp_pml.replace('.pml', '.mutations')
    with open(temp_pml_mutations, 'w') as f:
        for pos, wt_aa, mut_aa in mutations:
            f.write(f"{pos},{wt_aa},{mut_aa}\n")
    
    return temp_pml

def add_mutation_remarks_to_pdb(pdb_path, mutations):
    """
    Add REMARK lines to PDB file highlighting mutation sites for easier visualization.
    Inserts remarks before the ATOM records.
    """
    if not os.path.exists(pdb_path):
        return False
    
    try:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()
        
        # Find where to insert remarks (before first ATOM line)
        insert_idx = 0
        for idx, line in enumerate(lines):
            if line.startswith('ATOM'):
                insert_idx = idx
                break
        
        # Create mutation remarks
        remarks = []
        remarks.append("REMARK   MUTATION SITES\n")
        remarks.append("REMARK   ================\n")
        
        for pos, wt_aa, mut_aa in mutations:
            wt_aa_three = aa_single_to_three(wt_aa)
            mut_aa_three = aa_single_to_three(mut_aa)
            remark = f"REMARK   Position {pos}: {wt_aa_three} -> {mut_aa_three} ({wt_aa}{pos}{mut_aa})\n"
            remarks.append(remark)
        
        remarks.append("REMARK   ================\n")
        
        # Insert remarks at the appropriate location
        new_lines = lines[:insert_idx] + remarks + lines[insert_idx:]
        
        with open(pdb_path, 'w') as f:
            f.writelines(new_lines)
        
        return True
    except Exception as e:
        print(f"Warning: Could not add mutation remarks to {pdb_path}: {e}")
        return False

def add_mutation_remarks_to_pdb(pdb_path, mutations):
    """
    Add REMARK lines to PDB file highlighting mutation sites for easier visualization.
    Inserts remarks before the ATOM records.
    """
    if not os.path.exists(pdb_path):
        return False
    
    try:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()
        
        # Find where to insert remarks (before first ATOM line)
        insert_idx = 0
        for idx, line in enumerate(lines):
            if line.startswith('ATOM'):
                insert_idx = idx
                break
        
        # Create mutation remarks
        remarks = []
        remarks.append("REMARK   MUTATION SITES\n")
        remarks.append("REMARK   ================\n")
        
        for pos, wt_aa, mut_aa in mutations:
            wt_aa_three = aa_single_to_three(wt_aa)
            mut_aa_three = aa_single_to_three(mut_aa)
            remark = f"REMARK   Position {pos}: {wt_aa_three} -> {mut_aa_three} ({wt_aa}{pos}{mut_aa})\n"
            remarks.append(remark)
        
        remarks.append("REMARK   ================\n")
        
        # Insert remarks at the appropriate location
        new_lines = lines[:insert_idx] + remarks + lines[insert_idx:]
        
        with open(pdb_path, 'w') as f:
            f.writelines(new_lines)
        
        return True
    except Exception as e:
        print(f"Warning: Could not add mutation remarks to {pdb_path}: {e}")
        return False

def create_comparison_pml_script(gene, wt_pdb_path, mut_pdb_path, mutations, base_filename, interactive=False):
    """
    Creates a PML script for side-by-side comparison of wild-type and mutant structures.
    mutations: list of (position, wt_aa, mut_aa) tuples
    """
    pml_filename = f"{base_filename}.pml"
    png_filename = f"{base_filename}.png"
    
    wt_pdb_safe = wt_pdb_path.replace('\\', '/')
    mut_pdb_safe = mut_pdb_path.replace('\\', '/')
    
    # Get mutation positions for selection
    mut_positions = "+".join([str(pos) for pos, _, _ in mutations])
    
    # Create mutation labels text
    mutation_labels = ", ".join([f"{wt}{pos}{mut}" for pos, wt, mut in mutations])
    
    pml_content = f"""# Side-by-side comparison: Wild-type vs Mutant
# Mutations: {mutation_labels}

# Load structures
load {wt_pdb_safe}, wild_type
load {mut_pdb_safe}, mutant

# Align structures
align mutant, wild_type

# Basic visualization
hide everything
show cartoon, wild_type
show cartoon, mutant

# Color schemes
color lightblue, wild_type
color lightgreen, mutant
bg_color white

# Highlight mutation sites
select wt_mutations, wild_type and resi {mut_positions}
select mut_mutations, mutant and resi {mut_positions}

show spheres, wt_mutations
show spheres, mut_mutations

color yellow, wt_mutations
color firebrick, mut_mutations

# Set up side-by-side view
set grid_mode, 1
set grid_slot, 1, wild_type
set grid_slot, 2, mutant

# Zoom to show complete structures
zoom complete

# Center on mutations for better view
center wt_mutations

# Rotate for better 3D perspective
rotate y, 15

"""

    # Add labels for each mutation
    label_y_offset = 0
    for pos, wt, mut in mutations:
        pml_content += f"""# Label mutation {wt}{pos}{mut}
label wild_type and resi {pos} and name CA, "{wt}{pos}"
label mutant and resi {pos} and name CA, "{mut}{pos}"
set label_color, yellow, wild_type and resi {pos}
set label_color, red, mutant and resi {pos}
"""

    pml_content += f"""
# Render high-quality image
png {png_filename}, width=2400, height=1200, ray=1
"""
    
    if not interactive:
        pml_content += "quit\n"
        
    with open(pml_filename, 'w') as f:
        f.write(pml_content)
    return pml_filename

def find_pymol_executable(custom_path=None):
    """
    Attempts to find the PyMOL executable.
    """
    if custom_path:
        if os.path.exists(custom_path):
            return custom_path
        print(f"Error: Custom PyMOL path '{custom_path}' not found.")
        return None

    # Check PATH
    exe = shutil.which("pymol")
    if exe: return exe
    
    # Check common Windows paths
    common_paths = [
        r"C:\Program Files\PyMOL\PyMOLWin.exe",
        r"C:\Program Files (x86)\PyMOL\PyMOLWin.exe",
        os.path.expanduser(r"~\AppData\Local\Schrodinger\PyMOL2\PyMOLWin.exe"),
        os.path.expanduser(r"~\AppData\Local\Schrodinger\PyMOL2\pymol.exe")
    ]
    
    for path in common_paths:
        if os.path.exists(path):
            return path
            
    return None

def main():
    parser = argparse.ArgumentParser(description="Automate PyMOL mutation visualizations with side-by-side comparison.")
    parser.add_argument("input_file", nargs='?', help="Path to input CSV/Excel file.")
    parser.add_argument("-i", "--interactive", action="store_true", help="Keep PyMOL open after loading (do not quit).")
    parser.add_argument("--pymol-path", help="Path to PyMOL executable.")
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
            input_file = "Input/AMR_Full_Analysis.xlsx"

    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return

    output_dir = "PyMOL_Visuals"
    mutant_pdb_dir = os.path.join(PDB_CACHE_DIR, "Mutants")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(mutant_pdb_dir, exist_ok=True)

    # Check for PyMOL
    pymol_exe = find_pymol_executable(args.pymol_path)
    if not pymol_exe:
        print("Warning: PyMOL executable not found.")
        print("Please ensure PyMOL is in your PATH or provide it via --pymol-path.")
        print("Scripts will be generated but not executed (no PNGs).")
    else:
        print(f"Using PyMOL at: {pymol_exe}")

    print(f"Processing {input_file}...")
    print(f"Looking for alignment files in: {ALIGNMENT_DIR}")
    print(f"Looking for FASTA files in: {FASTA_DIR}")
    
    count = 0
    for mutation_data in parse_mutations(input_file):
        gene = mutation_data['gene_name']
        pdb_id = mutation_data['pdb_id']
        genome_id = mutation_data['genome_id']
        mutations = mutation_data['mutations']
        wt_sequence = mutation_data['wt_sequence']
        mut_sequence = mutation_data['mut_sequence']
        
        if not pdb_id:
            print(f"Skipping {gene}: No PDB ID mapped.")
            continue

        # Get wild-type PDB file (cached)
        wt_pdb_path = get_pdb_file(pdb_id)
        if not wt_pdb_path:
            print(f"Skipping {gene}: Could not download PDB {pdb_id}.")
            continue

        # Create mutation string for filename
        mutation_str = "_".join([f"{wt}{pos}{mut}" for pos, wt, mut in mutations])
        
        print(f"\n> Processing {genome_id} {gene}: {mutation_str}")
        print(f"  Found {len(mutations)} mutation(s)")
        
        # Generate mutant PDB filename
        mutant_pdb_filename = f"{genome_id}_{gene}_mutant.pdb"
        mutant_pdb_path = os.path.join(mutant_pdb_dir, mutant_pdb_filename)
        
        # Create base filename for outputs
        base_name = os.path.join(output_dir, f"{genome_id}_{gene}_{mutation_str}")
        
        # Step 1: Create mutation threading script
        mutate_pml = create_mutant_pdb_by_threading(wt_pdb_path, mutations, mutant_pdb_path)
        
        if pymol_exe:
            # Step 2: Run mutagenesis to create mutant structure
            print(f"  Creating mutant structure...")
            try:
                subprocess.run([pymol_exe, "-c", "-q", mutate_pml], check=True)
                time.sleep(1)
                
                # Clean up mutagenesis script
                if os.path.exists(mutate_pml):
                    try:
                        os.remove(mutate_pml)
                    except:
                        pass
                
                # Check if mutant PDB was created
                if not os.path.exists(mutant_pdb_path):
                    print(f"  Error: Mutant PDB not created. Skipping visualization.")
                    continue
                
                # Step 3: Add mutation remarks to the PDB file for clarity
                add_mutation_remarks_to_pdb(mutant_pdb_path, mutations)
                
                # Step 4: Create side-by-side comparison visualization
                print(f"  Generating comparison visualization...")
                comparison_pml = create_comparison_pml_script(
                    gene, wt_pdb_path, mutant_pdb_path, mutations, base_name, interactive=args.interactive
                )
                
                # Step 4: Run comparison visualization
                if not args.interactive:
                    subprocess.run([pymol_exe, "-c", "-q", comparison_pml], check=True)
                    time.sleep(1)
                    
                    # Clean up comparison script
                    if os.path.exists(comparison_pml):
                        try:
                            os.remove(comparison_pml)
                        except:
                            pass
                    
                    print(f"  ✓ Created: {base_name}.png")
                else:
                    print(f"  [Interactive] Created {comparison_pml} (Manual open required)")
                
                count += 1
                
            except subprocess.CalledProcessError as e:
                print(f"  Error running PyMOL: {e}")
            except Exception as e:
                print(f"  Error: {e}")
        else:
            print(f"  [Dry Run] Would create mutant and comparison for {gene} {genome_id}")
            count += 1

    print(f"\n{'='*60}")
    print(f"Done! Processed {count} gene-genome combinations.")
    print(f"Output directory: {output_dir}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
