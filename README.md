## Sq

### Overview
Sq computes the static structure factor S(q) from structure files using the Dynasor toolkit. It provides two CLI commands: `sq` to compute and save S(q) data, and `fwhm` to calculate the Full Width at Half Maximum of the structure factor peak.

### Features
- **Compute S(q)** from structure files (supports multiple formats via ASE)
- **Calculate FWHM** of the first peak in S(q)
- **Spherical averaging** over q-points
- **Configurable resolution** via `q_max`, `max_points`, `q_width`, and `nq`
- **Optional plotting** for FWHM visualization
- **CLI-first** workflow with sensible defaults

### Requirements
- **Python**: >= 3.11
- **Dependencies** (installed automatically): `ase`, `click`, `dynasor`, `gsd`, `numpy`, `matplotlib` (for plotting with `fwhm` command)

### Installation
- Using `pip`:
```bash
pip install .
```
- Using `uv` (recommended for speed and reproducibility):
```bash
uv pip install .
```

This will install two console scripts: `sq` and `fwhm`.

### Quick Start

#### Compute S(q)
Given a structure file (e.g., `test/tmp_traj.xyz`), run:
```bash
sq test/tmp_traj.xyz -o sq.dat
```
This writes two columns to `sq.dat`:
- Column 1: q (1/length units consistent with your input cell)
- Column 2: S(q)

#### Calculate FWHM
To calculate the FWHM of the first peak:
```bash
fwhm test/tmp_traj.xyz -o fwhm.dat
```

### CLI Usage

#### `sq` Command
```bash
sq STRUCTURE_FILE [options]
```

**Arguments:**
- **STRUCTURE_FILE**: Path to a structure file (supports formats readable by ASE, e.g., `.xyz`, `.cif`, `.vasp`, etc.)

**Options:**
- `-q, --q_max FLOAT`          Maximum q to sample (default: 4.0)
- `-m, --max_points INTEGER`   Maximum number of q-points for spherical sampling (default: 20000)
- `-w, --q_width FLOAT`        Gaussian smearing width in q for averaging (default: 0.02)
- `-n, --nq INTEGER`           Number of linearly spaced q points for the output (default: 1000)
- `-f, --format TEXT`          Format of the structure file (default: `extxyz`)
- `-o, --output_file TEXT`     Output filename (default: `sq.dat`)

**Example:**
```bash
sq input.xyz \
  --q_max 6.0 \
  --max_points 40000 \
  --q_width 0.03 \
  --nq 1500 \
  --format extxyz \
  --output_file my_sq.dat
```

#### `fwhm` Command
```bash
fwhm STRUCTURE_FILE [options]
```

**Arguments:**
- **STRUCTURE_FILE**: Path to a structure file (supports formats readable by ASE)

**Options:**
- `-q, --q_max FLOAT`              Maximum q to sample (default: 4.0)
- `-m, --max_points INTEGER`       Maximum number of q-points for spherical sampling (default: 20000)
- `-w, --q_width FLOAT`            Gaussian smearing width in q for averaging (default: 0.02)
- `-n, --nq INTEGER`               Number of linearly spaced q points (default: 1000)
- `-sf, --structure_format TEXT`   Format of the structure file (default: `extxyz`)
- `-s, --skip_front INTEGER`      Number of points to skip from the front when finding peak (default: 100)
- `-o, --output_file TEXT`         Output filename (default: `fwhm.dat`)
- `-p, --plot`                     Plot S(q) with FWHM markers and save as PNG (flag)

**Example:**
```bash
fwhm input.xyz \
  --q_max 6.0 \
  --skip_front 150 \
  --plot \
  --output_file my_fwhm.dat
```

This will:
- Calculate S(q) and find the FWHM of the first peak
- Print FWHM, left x, and right x values to console
- Save the values to `my_fwhm.dat`
- Generate a plot `my_fwhm.dat.png` showing S(q) with vertical lines marking the FWHM boundaries

### How It Works (brief)
Under the hood, both commands:
- Parse your structure file with ASE (supports multiple formats)
- Write a two-frame temporary trajectory `tmp_traj.xyz` (required by Dynasor)
- Build spherical q-points up to `q_max`
- Compute static structure factors with Dynasor
- Apply spherical averaging with Gaussian smearing of width `q_width`
- Normalize by number of atoms

For `sq`: Outputs `q, S(q)` to a text file.

For `fwhm`: Finds the maximum in S(q) (after skipping initial points), calculates the Full Width at Half Maximum, and optionally plots the result.

### Input Notes
- For physically meaningful S(q), ensure your structure file includes a realistic periodic cell.
- If your structure lacks cell info, results may not reflect a periodic bulk system.
- Supported formats include: `extxyz`, `xyz`, `cif`, `vasp`, and other formats supported by ASE.

### Output

#### `sq` Command Output
- Plain text file with two columns: `q  S(q)`
- Number of rows equals `nq`

#### `fwhm` Command Output
- Plain text file with three lines:
  - `FWHM: <value> A^-1`
  - `Left x: <value> A^-1`
  - `Right x: <value> A^-1`
- If `--plot` is used, also generates a PNG file showing S(q) with vertical lines marking FWHM boundaries

### Examples and Test Data
- A sample output file is provided at `test/sq.dat` (1001 lines).
- Example input can be placed at `test/tmp_traj.xyz` or use your own file.

### API (Python)
You can also call the underlying functions directly:
```python
from structure_factor.structure2sf import sq_of_structure, fwhm_structure_file

# Compute S(q)
q, Sq = sq_of_structure(
    structure_file="input.xyz",
    q_max=4.0,
    max_points=20000,
    q_width=0.02,
    nq=1000,
    structure_format="extxyz"
)

# Calculate FWHM
fwhm, xl, xr = fwhm_structure_file(
    structure_file="input.xyz",
    q_max=4.0,
    max_points=20000,
    q_width=0.02,
    nq=1000,
    structure_format="extxyz",
    skip_front=100
)
```

### Development
- Create a virtual environment and install in editable mode:
```bash
uv venv
source .venv/bin/activate
uv pip install -e .
```
- Run the CLI from source:
```bash
python -m structure_factor.structure2sf input.xyz -o sq.dat
```
  or simply use the installed `sq` and `fwhm` entrypoints.

### License
Specify your license here (e.g., MIT, Apache-2.0). If adding a license file, update this section accordingly.
