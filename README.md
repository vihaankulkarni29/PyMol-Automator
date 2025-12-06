# PyMOL Automator

A Python tool to automate the generation of 3D protein structure visualizations using PyMOL. It reads mutation data from a CSV/Excel file, maps genes to PDB structures, and generates PyMOL scripts to highlight specific mutations.

## Features

-   **Automated Visualization**: Generates `.pml` scripts and PNG images for batches of mutations.
-   **Smart Mapping**: Automatically maps gene names (e.g., `gyrA`, `emrE`) to their corresponding PDB structures.
-   **PDB Caching**: Downloads PDB files once and caches them locally to save bandwidth and improve speed.
-   **Interactive Mode**: Option to generate scripts that stay open for manual inspection.
-   **Headless Support**: Can run in the background to generate images without opening the PyMOL GUI.

## Prerequisites

-   **Python 3.x**
-   **PyMOL**: Must be installed. The script attempts to auto-detect it in standard locations (e.g., `C:\Program Files\PyMOL`, `AppData`). If not found, you can specify the path manually.

## Installation

1.  Clone the repository:
    ```bash
    git clone https://github.com/vihaankulkarni29/PyMol-Automator.git
    cd PyMol-Automator
    ```
2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage

### Basic Usage (Headless Image Generation)
Place your input file (e.g., `AMR_Full_Analysis.xlsx`) in the `Input/` directory or pass it as an argument.

```bash
python pymol_auto_visualizer.py Input/AMR_Full_Analysis.xlsx
```
This will generate PNG images in the `PyMOL_Visuals/` directory.

### Specifying PyMOL Path
If the script cannot find PyMOL automatically, use the `--pymol-path` argument:

```bash
python pymol_auto_visualizer.py --pymol-path "C:\Path\To\PyMOLWin.exe"
```

### Interactive Mode (Manual Inspection)
To generate scripts that you can open in PyMOL without them immediately closing:

```bash
python pymol_auto_visualizer.py -i
```
This will create `.pml` files in `PyMOL_Visuals/`. You can then double-click these files to view the mutation in 3D.

## Directory Structure

-   `pymol_auto_visualizer.py`: Main script.
-   `Input/`: Place your mutation data files here.
-   `PyMOL_Visuals/`: Generated output (images and scripts).
-   `PDB_Cache/`: Stores downloaded PDB structure files.

## Input Format

The input file (CSV or Excel) must contain:
-   **Gene Name**: e.g., `gyrA`, `parC_...`
-   **Variant**: e.g., `S83L`, `I219V`

## License
MIT
