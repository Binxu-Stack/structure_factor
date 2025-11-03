## Sq

### Overview
Sq computes the static structure factor S(q) from an XYZ structure using the Dynasor toolkit. It exposes a simple CLI `sq` to read an input `.xyz`, compute spherically averaged S(q) up to a chosen `q_max`, and write `q  S(q)` data to a text file.

### Features
- **Compute S(q)** from a single-frame XYZ file
- **Spherical averaging** over q-points
- **Configurable resolution** via `q_max`, `max_points`, `q_width`, and `nq`
- **CLI-first** workflow with sensible defaults

### Requirements
- **Python**: >= 3.11
- **Dependencies** (installed automatically): `ase`, `click`, `dynasor`, `gsd`, `numpy`

### Installation
- Using `pip`:
```bash
pip install .
```
- Using `uv` (recommended for speed and reproducibility):
```bash
uv pip install .
```

This will install a console script named `sq`.

### Quick Start
Given an input XYZ file (e.g., `test/tmp_traj.xyz`), run:
```bash
sq test/tmp_traj.xyz -o sq.dat
```
This writes two columns to `sq.dat`:
- Column 1: q (1/length units consistent with your input cell)
- Column 2: S(q)

### CLI Usage
```bash
sq XYZ_FILE [options]
```

#### Arguments
- **XYZ_FILE**: Path to an input `.xyz` file. Should contain one structure with a valid periodic cell in the header if you expect meaningful q-space results.

#### Options
- `-q, --q_max FLOAT`          Define the maximum q to sample (default: 4.0)
- `-m, --max_points INTEGER`   Maximum number of q-points for spherical sampling (default: 20000)
- `-w, --q_width FLOAT`        Gaussian smearing width in q for averaging (default: 0.02)
- `-n, --nq INTEGER`           Number of linearly spaced q points for the output (default: 1000)
- `-o, --output_file TEXT`     Output filename (default: `sq.dat`)

Example with custom parameters:
```bash
sq input.xyz \
  --q_max 6.0 \
  --max_points 40000 \
  --q_width 0.03 \
  --nq 1500 \
  --output_file my_sq.dat
```

### How It Works (brief)
Under the hood, `sq`:
- Parses your XYZ with ASE
- Writes a two-frame temporary trajectory `tmp_traj.xyz` (required by Dynasor)
- Builds spherical q-points up to `q_max`
- Computes static structure factors with Dynasor
- Applies spherical averaging with Gaussian smearing of width `q_width`
- Normalizes by number of atoms and writes `q, S(q)`

### Input Notes
- For physically meaningful S(q), ensure your XYZ includes a realistic periodic cell.
- If your XYZ lacks cell info, results may not reflect a periodic bulk system.

### Output
- Plain text file with two columns: `q  S(q)`
- Number of rows equals `nq`

### Examples and Test Data
- A sample output file is provided at `test/sq.dat` (1001 lines).
- Example input can be placed at `test/tmp_traj.xyz` or use your own file.

### API (Python)
You can also call the underlying function directly:
```python
from structure_factor.gsd2sf import sq_of_xyz

q, Sq = sq_of_xyz(
    xyz_file="input.xyz",
    q_max=4.0,
    max_points=20000,
    q_width=0.02,
    nq=1000,
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
python -m structure_factor.gsd2sf input.xyz -o sq.dat
```
  or simply use the installed `sq` entrypoint.

### License
Specify your license here (e.g., MIT, Apache-2.0). If adding a license file, update this section accordingly.
