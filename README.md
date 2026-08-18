# Fast Glycan Masking Pipeline

Fast placement of **Man5 N-linked glycans** onto one or more protein sites using a reusable conformer library.

This pipeline was designed as a faster alternative to repeatedly running Rosetta `GlycanTreeModeler`. It retains **Rosetta FastDesign** for engineering valid N-glycosylation sequons when needed, while replacing repeated glycan conformational modeling with a pre-generated Man5 rotamer library.

The updated version samples flexibility at **two levels**:

1. **ASN → first GlcNAc attachment orientation**
2. **GlcNAc/Man → downstream glycan φ/ψ/ω torsions**

This allows the glycan root orientation as well as its downstream branches to adapt to different protein environments.

---

# 1. Overview

The pipeline consists of three main stages:

```text
Reference glycoprotein containing Man5
                │
                ▼
     glycan_rotamer_generator
                │
                │ Sample:
                │ • ASN→GlcNAc attachment torsion
                │ • glycan φ/ψ
                │ • ω for 1→6 linkages
                │ • reject glycan self-clashes
                ▼
        Man5_library.npz
                │
                ▼
      glycan_library_placer
                │
                │ • align library to target ASN
                │ • protein–glycan clash filtering
                │ • glycan–glycan clash filtering
                ▼
     Glycosylated protein models
```

For newly engineered glycosylation sites:

```text
Protein
   │
   ▼
Rosetta FastDesign
   │
   │ Create N-X-S/T sequon
   │ Optimize ASN/local side chain
   ▼
Designed protein
   │
   ▼
Man5 library placement
   │
   ▼
Glycosylated models
```

---

# 2. Package Structure

Recommended directory structure:

```text
fast_man5_pipeline/
│
├── fast_glycan_masking/
│   ├── __init__.py
│   ├── glycan_rotamer_generator.py
│   ├── glycan_library_placer.py
│   ├── rotamer_library_generator.py
│   └── fast_glycan_masking.py
│
├── configs/
│   └── man5_sampling.yaml
│
├── examples/
│   ├── Glyc_Des_head_6uig.pdb
│   └── Des_head_6uig_cut.pdb
│
├── Man5_library.npz
├── template_FastDesign.xml
├── fast_design.sh
└── README.md
```

The `fast_glycan_masking/` directory must contain an `__init__.py` file because the scripts use package-relative imports.

Commands should therefore be run from the project root using:

```bash
python -m fast_glycan_masking.<module>
```

rather than running the individual Python files directly.

---

# 3. Dependencies

The Python components require:

```text
Python 3
NumPy
PyYAML
```

Install them, for example, with:

```bash
pip install numpy pyyaml
```

Rosetta is additionally required if `fast_glycan_masking.py` is used to create new N-X-S/T sequons with FastDesign.

---

# 4. Man5 Conformer Library Generation

## Purpose

`glycan_rotamer_generator.py` generates a reusable Man5 conformer library from a reference Rosetta glycoprotein PDB.

The reference PDB must contain:

* an ASN carrying the Man5 glycan;
* appropriate `LINK` records describing the ASN–glycan and sugar–sugar connectivity;
* the complete Man5 structure.

The generator identifies the glycan attached to the requested ASN and reconstructs the glycan tree from the `LINK` records.

---

# 5. What Is Sampled?

The updated generator samples two types of conformational flexibility.

## 5.1 ASN → First GlcNAc Attachment

The previous implementation retained the ASN→first-sugar geometry directly from the template. Therefore, all library conformers emerged from the ASN in essentially the same root orientation.

The updated implementation additionally samples the attachment torsion:

```text
CG(ASN) – ND2(ASN) – C1(GlcNAc) – O5(GlcNAc)
```

Rotation occurs around the:

```text
ND2 – C1
```

bond.

Conceptually:

```text
                 glycan
                    |
                    |
ASN–CG–ND2────────C1
          ↖       ↗
             ↻
      attachment torsion
```

Changing this torsion rotates the entire glycan relative to ASN.

Importantly, the implementation does **not** randomly perturb the covalent geometry.

The following remain fixed:

```text
ND2–C1 bond length
CG–ND2–C1 bond angle
ND2–C1–O5 bond angle
```

Only the attachment dihedral is sampled.

This allows different root orientations without introducing arbitrary distortion into the ASN–GlcNAc linkage.

---

## 5.2 Downstream Glycan φ/ψ/ω

Sugar–sugar torsions continue to be sampled independently.

For each glycosidic linkage, the generator identifies:

```text
φ
ψ
```

and, for `1→6` linkages:

```text
ω
```

The sampled torsions are applied to the complete downstream glycan subtree.

Therefore, each library member can differ in both:

```text
ASN → root GlcNAc orientation
              +
downstream Man5 conformation
```

For example:

```text
Conformer 1

ASN ─── GlcNAc
              \
               GlcNAc ─ Man
                         / \
                       Man Man


Conformer 2

ASN ─────── GlcNAc
             |
             GlcNAc
                 \
                  Man
                 /  \
               Man  Man
```

---

# 6. Sampling Configuration

Sampling distributions are controlled through a YAML configuration file.

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

## Attachment sampling

A recommended initial setting is:

```yaml
attachment:
  mode: template_jitter
  sd: 20.0
```

This samples around the attachment orientation observed in the reference glycoprotein.

For example, if the reference attachment torsion is approximately:

```text
30°
```

sampling is centered around that orientation with a standard deviation of:

```text
20°
```

This is generally preferable to unrestricted rotation because the reference structure provides a chemically reasonable starting geometry.

A larger value such as:

```yaml
attachment:
  mode: template_jitter
  sd: 30.0
```

can be tested if a target site has severe steric restrictions.

---

## Sampling modes

Three sampling modes are supported.

### Template jitter

```yaml
mode: template_jitter
sd: 20.0
```

Samples around the torsion observed in the reference structure.

An optional offset can also be supplied:

```yaml
mode: template_jitter
offset: 10.0
sd: 20.0
```

### Normal

```yaml
mode: normal
mean: 70.0
sd: 15.0
```

Samples from a normal distribution centered at an explicitly specified absolute angle.

### Discrete

```yaml
mode: discrete
values: [-60, -30, 0, 30, 60]
weights: [1, 2, 3, 2, 1]
```

Samples from explicitly specified absolute torsion values.

For most applications, `template_jitter` is the recommended starting point.

---

# 7. Generate the Man5 Library

Example:

```bash
python -m fast_glycan_masking.glycan_rotamer_generator \
    --pdb examples/Glyc_Des_head_6uig.pdb \
    --chain A \
    --resid 200 \
    --config configs/man5_sampling.yaml \
    --n-conformers 10000 \
    --out-prefix Man5_library
```

Here:

```text
--pdb
```

is the reference glycoprotein containing Man5.

```text
--chain A
--resid 200
```

identify the ASN carrying the reference glycan.

```text
--n-conformers 10000
```

requests 10,000 clash-free conformers.

The older option:

```text
--n-models
```

is retained as a backward-compatible alias for `--n-conformers`.

---

# 8. Library Generation Outputs

The generator writes:

```text
Man5_library.pdb
Man5_library.npz
```

## `Man5_library.pdb`

A multi-model PDB containing the accepted Man5 conformations.

This file is useful for visualization and inspection.

## `Man5_library.npz`

The NumPy archive is the library used by the placement program.

It contains information including:

```text
coords
atom_names
residue_numbers
residue_names
chains
elements
attachment_frame
root_residue_number
linkage information
torsion names
torsion types
template angles
sampled angles
attachment angles
```

The attachment frame is derived from:

```text
ASN CB
ASN CG
ASN ND2
```

of the reference glycoprotein.

This frame allows each library conformer to be transferred onto another ASN.

---

# 9. Self-Clash Filtering During Library Generation

After the attachment and downstream torsions are sampled, the resulting Man5 structure is checked for internal clashes.

Example options include:

```bash
--clash-cutoff 2.0
--exclude-bonds 3
```

Atom pairs separated by a small number of covalent bonds are excluded from the self-clash calculation.

Conformations containing unacceptable internal clashes are rejected.

Generation continues until either:

```text
requested number of conformers
```

is reached or:

```text
--max-attempts
```

is exceeded.

---

# 10. Place Man5 at a Target Site

Once the library has been generated, it can be reused at many target positions.

Example:

```bash
python -m fast_glycan_masking.glycan_library_placer \
    --library Man5_library.npz \
    --protein examples/Des_head_6uig_cut.pdb \
    --site A:200 \
    --n-models 10 \
    --out-dir outputs/A200 \
    --prefix Glyc_Des_A200
```

---

# 11. N-X-S/T Requirement

Every target position must already form a valid N-linked glycosylation sequon:

```text
N-X-S/T
```

where:

```text
X != Proline
```

For example:

```text
N-A-T
N-G-S
N-V-T
```

are valid.

The placer validates the sequon before attempting placement.

If the requested site is not valid, use `fast_glycan_masking.py` to create the sequon using Rosetta FastDesign.

---

# 12. How Placement Works

For every target ASN, the placer extracts:

```text
CB
CG
ND2
```

The same three atoms were stored from the source ASN when the library was generated.

A Kabsch alignment is calculated between:

```text
source ASN CB/CG/ND2
          ↓
target ASN CB/CG/ND2
```

The complete glycan conformer is then transformed using that rotation and translation.

Importantly, the conformer's sampled ASN→GlcNAc orientation is preserved during this transformation.

Therefore, the target site is tested against conformers containing different:

```text
root attachment orientations
+
downstream Man5 conformations
```

rather than only different downstream Man5 conformations.

---

# 13. Protein–Glycan Clash Filtering

After alignment, every library conformer is checked against the target protein.

Two clash modes are supported.

## van der Waals mode

Recommended:

```bash
--clash-mode vdw \
--vdw-scale 0.70
```

For two atoms with van der Waals radii:

```text
r1
r2
```

a clash occurs approximately when:

```text
distance < scale × (r1 + r2)
```

A smaller `vdw-scale` is more permissive.

For example:

```text
0.70
```

is more permissive than:

```text
0.80
```

The current pipeline uses `0.70` as a practical default.

## Distance mode

Alternatively:

```bash
--clash-mode distance \
--clash-cutoff 2.0
```

This treats atom pairs closer than the specified distance as clashes.

---

# 14. Why Attachment Sampling Matters

Consider a site where the template-derived root orientation points directly into the protein:

```text
            protein
          █████████
ASN ───── NAG██████
              \
               Man5
```

If only downstream φ/ψ/ω are sampled, changing the outer Man residues may not resolve the collision near the first GlcNAc.

With attachment sampling, another conformer may instead have:

```text
                 NAG ─ Man5
                /
ASN ───────────
        ███████
        protein
```

The library therefore has an additional way to escape local steric constraints.

This can be especially important when the same reusable Man5 library is transferred to many structurally different protein sites.

---

# 15. Multi-Site Glycan Placement

Multiple sites can be supplied by repeating `--site`.

Example:

```bash
python -m fast_glycan_masking.glycan_library_placer \
    --library Man5_library.npz \
    --protein examples/protein.pdb \
    --site A:100 \
    --site A:120 \
    --site A:200 \
    --n-models 100 \
    --out-dir outputs/multisite \
    --prefix Glyc_multi
```

Each output model independently selects one library conformer for each site.

For example:

```text
Model 1
A100 → conformer 238
A120 → conformer 8174
A200 → conformer 1230

Model 2
A100 → conformer 6201
A120 → conformer 931
A200 → conformer 7452
```

The program checks both:

```text
protein ↔ glycan clashes
```

and:

```text
glycan ↔ glycan clashes
```

before accepting the combination.

Thus native and engineered glycosylation sites can be modeled simultaneously.

---

# 16. Existing Glycans

By default, existing carbohydrate residues in the input protein are removed before placement.

This allows native and newly engineered glycosylation sites to be rebuilt consistently from the same reusable library.

The behavior corresponds conceptually to the previous workflow's:

```text
strip_existing = 1
```

---

# 17. Placement Outputs

For a command such as:

```bash
--out-dir outputs/A200
--prefix Glyc_Des_A200
```

the output directory contains individual glycosylated PDB models such as:

```text
Glyc_Des_A200_0001.pdb
Glyc_Des_A200_0002.pdb
...
```

as well as placement metadata.

The updated placement metadata includes the selected attachment orientation for each successful model/site.

The NPZ output includes:

```text
selected_library_indices
selected_attachment_angles
site_labels
attachment_rmsd
coords
```

The manifest also records selected attachment angles in degrees.

This is useful when comparing easy and difficult glycosylation sites.

For example, if successful conformers at one site consistently use attachment orientations around a particular range, this may indicate that the local protein environment strongly constrains the root glycan orientation.

---

# 18. Creating New Glycosylation Sites with FastDesign

`fast_glycan_masking.py` combines Rosetta FastDesign with fast Man5 placement.

Example:

```bash
python -m fast_glycan_masking.fast_glycan_masking \
    --protein examples/protein.pdb \
    --library Man5_library.npz \
    --new-site A:200 \
    --n-models 100 \
    --out-dir outputs/A200 \
    --prefix Glyc_Des_A200
```

For a requested new site at position `i`, the wrapper ensures:

```text
i     = N
i + 1 != P
i + 2 = S or T
```

If necessary:

```text
position i     → N
position i+1   → A if originally P
position i+2   → T if not already S/T
```

A Rosetta resfile is generated automatically and FastDesign is run.

The designed structure is then passed to the library placer.

---

# 19. Native + New Sites

Native sites can be supplied separately using:

```text
--native-site
```

while engineered positions use:

```text
--new-site
```

For example:

```bash
python -m fast_glycan_masking.fast_glycan_masking \
    --protein examples/protein.pdb \
    --library Man5_library.npz \
    --native-site A:100 \
    --native-site A:120 \
    --new-site A:200 \
    --n-models 100 \
    --out-dir outputs/A200 \
    --prefix Glyc_Des_A200
```

The conceptual result is:

```text
WT native sites:
A100
A120

requested new masking site:
A200

Final modeled sites:
A100 + A120 + A200
```

All three glycans are sampled and placed in each final model.

This preserves native glycan masking positions while allowing additional engineered sites to be evaluated.

---

# 20. Recommended Workflow

## Step 1 — Generate the library once

```bash
python -m fast_glycan_masking.glycan_rotamer_generator \
    --pdb examples/Glyc_Des_head_6uig.pdb \
    --chain A \
    --resid 200 \
    --config configs/man5_sampling.yaml \
    --n-conformers 10000 \
    --out-prefix Man5_library
```

The resulting library contains variation in:

```text
ASN→GlcNAc attachment
+
Man5 φ/ψ/ω
```

## Step 2 — Reuse the library

For each target protein/site:

```bash
python -m fast_glycan_masking.glycan_library_placer \
    --library Man5_library.npz \
    --protein examples/Des_head_6uig_cut.pdb \
    --site A:200 \
    --n-models 10 \
    --out-dir outputs/A200 \
    --prefix Glyc_Des_A200
```

## Step 3 — Use FastDesign when necessary

If the target does not already contain an N-X-S/T sequon, run the wrapper:

```bash
python -m fast_glycan_masking.fast_glycan_masking \
    --protein examples/protein.pdb \
    --library Man5_library.npz \
    --new-site A:200 \
    --n-models 100
```

---

# 21. Diagnosing Failed Sites

A useful output line is:

```text
[info] A200: attachment RMSD=... A; XXXX/10000 pass protein clash filter
```

There are several possible cases.

## Many conformers survive

For example:

```text
4500/10000 pass protein clash filter
```

The site is geometrically accessible and placement should be straightforward.

## Few conformers survive

For example:

```text
20/10000 pass protein clash filter
```

The site is strongly constrained.

Inspect:

* selected attachment angles;
* glycan conformations;
* local protein environment;
* `vdw-scale`.

## Zero conformers survive

```text
0/10000 pass protein clash filter
```

The program will report that no library conformers pass the protein clash filter.

Possible causes include:

* the site is buried;
* the ASN side-chain orientation is unfavorable;
* attachment sampling is too narrow;
* the clash threshold is too strict;
* the local backbone environment cannot accommodate Man5;
* the designed ASN requires additional protein-side-chain optimization.

A useful diagnostic experiment is to increase:

```yaml
attachment:
  mode: template_jitter
  sd: 30.0
```

and regenerate the library.

If this produces surviving conformers where the narrower library produced none, root orientation was likely an important limiting factor.

Do not automatically keep increasing attachment flexibility indefinitely; the goal is to sample plausible conformations rather than force every requested site to accept a glycan.

---

# 22. Reproducibility

Both generation and placement use NumPy random number generators with a user-controlled seed.

For example:

```bash
--seed 2026
```

Using the same:

```text
input
configuration
library
seed
```

allows the stochastic sampling to be reproduced.

---

# 23. Important Design Principle

The responsibilities of the pipeline are intentionally separated:

```text
Rosetta FastDesign
│
├── protein sequence
├── N-X-S/T creation
└── local ASN/protein geometry

        ↓

Glycan rotamer generator
│
├── ASN→GlcNAc attachment torsion
├── glycan φ/ψ
├── glycan ω
└── glycan self-clash filtering

        ↓

Glycan library placer
│
├── transfer conformers to target ASN
├── protein–glycan clash filtering
├── glycan–glycan clash filtering
└── multi-site model generation
```

This avoids repeatedly performing expensive full glycan modeling for every candidate masking position.

---

# 24. Summary

The updated fast glycan masking workflow expands the original reusable-library approach by introducing flexibility at the protein–glycan interface.

Previously:

```text
ASN→GlcNAc orientation = fixed
glycan φ/ψ/ω           = sampled
```

Updated:

```text
ASN→GlcNAc torsion     = sampled
glycan φ/ψ/ω           = sampled
```

while:

```text
bond lengths           = fixed
bond angles            = fixed
```

The resulting library therefore captures both **root orientation** and **downstream Man5 conformational diversity**.

The library can then be reused across many candidate masking positions, with each target site selecting conformers that survive local protein–glycan and glycan–glycan clash filtering.
