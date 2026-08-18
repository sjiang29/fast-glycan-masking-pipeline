# Fast Glycan MAN5 Masking

A fast conformer-library approach for generating ensembles of **Man5-glycosylated antigen structures** for antibody docking.

This project replaces the computationally expensive **Rosetta GlycanTreeModeler** ensemble-generation step with a reusable Man5 conformer library. The library is generated once and can then be rapidly placed onto one or more N-glycosylation sites.

The resulting glycosylated antigen ensembles can be used for downstream antibody docking and glycan masking studies.

---

## Workflow

Original workflow:

```text
Protein
    │
    ▼
Rosetta FastDesign
    │
    ▼
SimpleGlycosylateMover
    │
    ▼
GlycanTreeModeler
    │
    ▼
Glycosylated antigen ensemble
    │
    ▼
Antibody docking
```

Fast-library workflow:

```text
Reference Man5 glycoprotein
          │
          ▼
Generate reusable Man5 library
          │
          ▼
FastDesign if needed
          │
          ▼
Place Man5 conformers
          │
          ▼
Protein/glycan clash filtering
          │
          ▼
Glycosylated antigen ensemble
          │
          ▼
Antibody docking
```

The Man5 conformer library only needs to be generated once and can be reused for different proteins and glycosylation sites.

---

# Repository Structure

```text
fast-glycan-man5-masking/
│
├── README.md
├── DOMAIN_KNOWLEDGE.md
├── LICENSE
├── requirements.txt
├── .gitignore
│
├── configs/
│   └── man5_sampling.yaml
│
├── src/
│   └── fast_glycan_masking/
│       ├── __init__.py
│       ├── rotamer_library_generator.py
│       ├── glycan_rotamer_generator.py
│       ├── glycan_library_placer.py
│       └── fast_glycan_masking.py
│
├── scripts/
│   ├── generate_man5_library.sh
│   └── place_multiple_glycans.sh
│
├── examples/
│   ├── reference_glycoprotein.pdb
│   ├── designed_antigen.pdb
│   └── README.md
│
├── libraries/
├── outputs/
└── tests/
```

---

# Code Structure

### `rotamer_library_generator.py`

Generic geometry utilities used by the glycan-specific code.

Includes functions for:

- bond detection
- bond-graph analysis
- dihedral calculation
- torsion rotation
- rotation matrices
- PDB writing

This file provides the underlying geometry engine and is imported by `glycan_rotamer_generator.py`.

---

### `glycan_rotamer_generator.py`

Generates the reusable Man5 conformer library from a reference glycoprotein.

It:

- identifies the Man5 glycan attached to a selected ASN,
- reconstructs glycan connectivity from PDB `LINK` records,
- samples the ASN → first GlcNAc attachment torsion,
- samples downstream glycan φ/ψ/ω torsions,
- rejects glycan self-clashes,
- writes the reusable `.pdb` and `.npz` libraries.

---

### `glycan_library_placer.py`

Places library conformers onto one or more target ASN residues.

It:

- validates N-X-S/T sequons,
- aligns each Man5 conformer to the target ASN,
- filters protein–glycan clashes,
- filters glycan–glycan clashes,
- supports simultaneous multi-site placement,
- writes complete glycosylated protein structures.

---

### `fast_glycan_masking.py`

Wrapper for engineered glycosylation sites.

It:

- identifies mutations needed to create an N-X-S/T sequon,
- generates a Rosetta resfile,
- runs the existing FastDesign workflow,
- validates the designed sites,
- calls `glycan_library_placer.py` for final glycan placement.

Native and newly engineered glycosylation sites can be modeled together.

---

# Installation

Clone the repository:

```bash
git clone https://github.com/sjiang29/fast-glycan-man5-masking.git
cd fast-glycan-man5-masking
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Add the source directory to the Python path:

```bash
export PYTHONPATH="$PWD/src:${PYTHONPATH:-}"
```

The package uses relative imports, so run modules from the repository root using:

```bash
python -m fast_glycan_masking.<module>
```

Do not run files inside `src/fast_glycan_masking/` directly.

---

# Stage 1 — Generate the Man5 Library

This step normally only needs to be performed once.

```bash
python -m fast_glycan_masking.glycan_rotamer_generator \
    --pdb examples/Glyc_Des_head_6uig.pdb \
    --chain A \
    --resid 200 \
    --config configs/man5_sampling.yaml \
    --n-conformers 10000 \
    --max-attempts 200000 \
    --seed 2026 \
    --out-prefix libraries/Man5_library
```

Outputs:

```text
libraries/Man5_library.npz
libraries/Man5_library.pdb
```

`--n-conformers` specifies the number of accepted glycan conformers stored in the reusable library.

`--max-attempts` specifies the maximum number of sampling attempts. It should normally be larger than `--n-conformers` because self-clashing structures are rejected.

---

## Quick Test

```bash
python -m fast_glycan_masking.glycan_rotamer_generator \
    --pdb examples/Glyc_Des_head_6uig.pdb \
    --chain A \
    --resid 200 \
    --config configs/man5_sampling.yaml \
    --n-conformers 100 \
    --max-attempts 5000 \
    --out-prefix outputs/Man5_library_test
```

---

## Inspect Glycan Topology

To inspect the detected glycan tree and torsions without generating a library:

```bash
python -m fast_glycan_masking.glycan_rotamer_generator \
    --pdb examples/Glyc_Des_head_6uig.pdb \
    --chain A \
    --resid 200 \
    --config configs/man5_sampling.yaml \
    --report-only
```

---

# Stage 2 — Prepare Glycosylation Sites

A target N-linked glycosylation site must satisfy:

```text
N-X-S/T
```

where `X` cannot be Proline.

If the sequon already exists, FastDesign is not required.

If the sequon does not exist, use the existing Rosetta FastDesign workflow or the `fast_glycan_masking.py` wrapper to introduce the required mutations before placement.

---

# Stage 3 — Place Man5 Glycans

## Single Site

```bash
python -m fast_glycan_masking.glycan_library_placer \
    --library libraries/Man5_library.npz \
    --protein examples/designed_antigen.pdb \
    --site A:87 \
    --n-models 100 \
    --out-dir outputs/A87 \
    --prefix Glyc_Des_A87
```

---

## Multiple Sites

Repeat `--site` for each glycosylation position:

```bash
python -m fast_glycan_masking.glycan_library_placer \
    --library libraries/Man5_library.npz \
    --protein examples/designed_antigen.pdb \
    --site A:87 \
    --site A:200 \
    --n-models 100 \
    --out-dir outputs/A87_A200 \
    --prefix Glyc_Des_A87_A200
```

For every output model, the program:

1. selects one library conformer for each site,
2. aligns each glycan to its target ASN,
3. checks protein–glycan clashes,
4. checks glycan–glycan clashes,
5. writes only clash-free models.

---

# Stage 4 — FastDesign + Placement

For a newly engineered site:

```bash
python -m fast_glycan_masking.fast_glycan_masking \
    --protein examples/antigen.pdb \
    --library libraries/Man5_library.npz \
    --new-site A:200 \
    --n-models 100 \
    --out-dir outputs/A200 \
    --prefix Glyc_Des_A200
```

Native sites can also be included:

```bash
python -m fast_glycan_masking.fast_glycan_masking \
    --protein examples/antigen.pdb \
    --library libraries/Man5_library.npz \
    --native-site A:100 \
    --native-site A:120 \
    --new-site A:200 \
    --n-models 100 \
    --out-dir outputs/A100_A120_A200 \
    --prefix Glyc_multi
```

This produces models containing glycans at all requested sites.

---

# Library Files

## `Man5_library.pdb`

Multi-model PDB containing generated Man5 conformers.

Useful for:

- PyMOL
- ChimeraX
- visualization
- manual inspection
- debugging

## `Man5_library.npz`

Binary NumPy library used during placement.

It stores information including:

- coordinates
- atom names
- residue names
- residue numbers
- elements
- glycan topology
- source ASN attachment frame
- sampled glycan torsions
- sampled ASN–GlcNAc attachment angles

The `.npz` format allows much faster reuse than repeatedly parsing large multi-model PDB files.

---

# Sampling Configuration

Sampling behavior is controlled by:

```text
configs/man5_sampling.yaml
```

Example:

```yaml
attachment:
  mode: template_jitter
  sd: 20.0

default:
  phi:
    mode: template_jitter
    sd: 15.0

  psi:
    mode: template_jitter
    sd: 15.0

  omega:
    mode: template_jitter
    sd: 15.0
```

Linkage-specific distributions can override these defaults.

See `DOMAIN_KNOWLEDGE.md` for the scientific interpretation of the attachment and glycosidic torsions.

---

# Clash Filtering

The recommended placement mode uses scaled van der Waals radii:

```bash
--clash-mode vdw \
--vdw-scale 0.70
```

A distance-based mode is also available:

```bash
--clash-mode distance \
--clash-cutoff 2.0
```

If no library conformers survive at a site, inspect:

- local protein accessibility,
- ASN orientation,
- attachment sampling,
- glycan conformational sampling,
- clash-filter settings.

---

# Current Limitations

- Currently focused on **N-linked Man5**.
- No Rosetta minimization is performed after library placement.
- Intended primarily for generating ensembles for glycan masking and antibody docking.
- Sampling quality depends on the torsion distributions provided in the YAML configuration.

---

# Future Work

Potential extensions include:

- Man9
- hybrid glycans
- complex glycans
- post-placement structural refinement
- GPU-accelerated clash detection
- torsion distributions derived from experimental glycan structures

---

# Acknowledgements

This project was developed in the Meiler Laboratory at Vanderbilt University.

The original generic rotamer-generation framework for small molecules was developed by a laboratory colleague and adapted here for branched glycan conformational sampling and glycan masking applications.