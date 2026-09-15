# Fast Glycan Masking Library from Existing Rosetta Models

`fast_glycan_masking_from_existing_models` builds a reusable N-glycan conformer
library from glycosylated protein structures that were already generated and
relaxed with Rosetta.

This pipeline is intended for deep glycan masking datasets organized into
`pos_*` folders. Instead of generating thousands of glycan conformers again
from a single template, it reuses the conformational information already
contained in the expensive Rosetta calculations.

---

# 1. Usage

## 1.1 Expected input directory

The input directory should contain folders named `pos_0`, `pos_1`, ...,
`pos_N`.

For example:

```text
all_positions/
├── pos_0/
│   ├── Glyc_score.sc
│   ├── Glyc_Des_head_..._0001.pdb
│   ├── Glyc_Des_head_..._0002.pdb
│   └── ...
├── pos_1/
│   ├── Glyc_score.sc
│   ├── Glyc_Des_head_..._0001.pdb
│   └── ...
├── pos_2/
│   └── ...
├── pos_3/
│   ├── Glyc_score.sc
│   ├── Glyc_Des_head_..._0001.pdb
│   └── ...
└── ...
```

`pos_0` is the wild-type reference used for the Rosetta energy baseline.

For `pos_i`, `i` is the amino-acid index where the additional N-glycan was
introduced during deep glycan masking.

Only PDB files beginning with `Glyc` are treated as glycosylated Rosetta
models.

There is **no assumption that a folder contains exactly 100 models**:

- 0 `Glyc*.pdb` files: skip the position.
- 83 models: use all 83.
- 100 models: use all 100.
- 147 models: use all 147.

Each folder should contain `Glyc_score.sc`. The pipeline uses the
`total_score` column and matches its model/description field to the actual PDB
files present in that folder.

---

## 1.2 Recommended: use the shell script

Edit:

```bash
run_existing_models_library.sh
```

The main settings are near the top:

```bash
MODE="empirical"

INPUT_DIR="/path/to/all_pos_folders"
OUT_PREFIX="libraries/Man5_empirical"

REU_DELTA="0"

NATIVE_POSITIONS=""
INCLUDE_WT_GLYCANS="false"
```

Then run:

```bash
bash run_existing_models_library.sh
```

The script supports three modes:

```bash
MODE="analyze"
MODE="empirical"
MODE="expanded"
```

### Analyze

Use:

```bash
MODE="analyze"
```

This only compares the Rosetta energies of all `pos_*` folders with `pos_0`.
No glycan library is generated.

The output is:

```text
position_scores.csv
```

This is recommended as the first run on a new dataset.

### Empirical library

Use:

```bash
MODE="empirical"
```

This performs the energy analysis, selects acceptable positions, extracts all
N-linked glycans from their Rosetta models, aligns them into a common ASN
attachment frame, and creates the empirical conformer library.

For example:

```bash
INPUT_DIR="/data/EEEV_DGS"
OUT_PREFIX="libraries/EEEV_Man5_empirical"
REU_DELTA="0"
```

### Expanded library

Use:

```bash
MODE="expanded"
```

This first builds the empirical library and then creates additional conformers
from each extracted Rosetta glycan.

Configure:

```bash
SAMPLES_PER_GLYCAN="4"
ATTACHMENT_SIGMA_DEG="20"
GLYCOSIDIC_SIGMA_DEG="20"
```

The original empirical conformer is always retained.

For large Rosetta datasets, building and examining the empirical library first
is recommended before enabling expansion.

---

## 1.3 Run directly with Python

### Analyze positions

```bash
python -m fast_glycan_masking_from_existing_models analyze \
    --input-dir /path/to/all_pos_folders \
    --reu-delta 0 \
    --out position_scores.csv
```

### Build an empirical library

```bash
python -m fast_glycan_masking_from_existing_models build \
    --input-dir /path/to/all_pos_folders \
    --reu-delta 0 \
    --out-prefix libraries/Man5_empirical
```

If the known native N-glycan positions are 141 and 200:

```bash
python -m fast_glycan_masking_from_existing_models build \
    --input-dir /path/to/all_pos_folders \
    --reu-delta 0 \
    --native-positions 141,200 \
    --out-prefix libraries/Man5_empirical
```

The native-position list is useful for provenance/validation. Glycan discovery
itself is based on the actual PDB connectivity.

### Include WT glycans

By default, `pos_0` supplies the WT energy baseline but its glycans are not
added to the conformer library.

To include them:

```bash
--include-wt-glycans
```

### Build an expanded library

```bash
python -m fast_glycan_masking_from_existing_models build \
    --input-dir /path/to/all_pos_folders \
    --reu-delta 0 \
    --samples-per-glycan 4 \
    --attachment-sigma-deg 20 \
    --glycosidic-sigma-deg 20 \
    --out-prefix libraries/Man5_expanded
```

### Limit final library size

For very large datasets:

```bash
--max-library-size 50000
```

A reproducible random subset is selected using `--seed`.

### Write a multi-model PDB

The NPZ library is always written.

If a human-readable multi-model PDB is also wanted:

```bash
--write-pdb
```

This is disabled by default because a library containing tens of thousands of
glycans can produce a very large PDB.

---

## 1.4 Output files

For:

```bash
--out-prefix libraries/Man5_empirical
```

the main outputs are:

```text
libraries/
├── Man5_empirical.npz
├── Man5_empirical_positions.csv
├── Man5_empirical_provenance.csv
├── Man5_empirical_manifest.json
└── Man5_empirical_skipped.tsv       # only when failures occur
```

If `--write-pdb` is supplied:

```text
Man5_empirical.pdb
```

### `*_positions.csv`

Contains position-level Rosetta filtering information, including:

```text
position
n_pdb
n_score_rows
n_matched
mean_total_score
wt_mean_total_score
delta_reu
selected
status
```

### `*_provenance.csv`

Each conformer retains its source information:

```text
source_folder
source_position
source_pdb
source_model
source_total_score

position_mean_total_score
wt_mean_total_score
position_delta_reu

attachment_chain
attachment_resid
glycan_role

conformer_origin
parent_empirical_index
```

`glycan_role` can distinguish the newly engineered site from native N-linked
sites.

`conformer_origin` distinguishes:

```text
rosetta_empirical
torsion_expanded
```

This prevents experimentally/physically sampled Rosetta conformers from being
mixed conceptually with conformers generated later by the fast sampler.

---

# 2. Pipeline Logic

## 2.1 Why this pipeline exists

The original fast glycan masking library generator starts from a glycosylated
template structure and creates diversity by changing glycan torsions.

The deep glycan masking calculations already produced a different source of
information: thousands of glycosylated proteins that were explicitly modeled
and relaxed by Rosetta.

Those calculations are expensive, but once they exist they can be reused.

The new strategy is therefore:

```text
existing Rosetta glycoproteins
              ↓
identify energetically acceptable positions
              ↓
extract their relaxed glycan conformations
              ↓
put all glycans into one common coordinate frame
              ↓
build reusable empirical conformer library
              ↓
optional additional torsion sampling
```

---

## 2.2 Step 1 — Establish the WT Rosetta baseline

`pos_0` is treated as the WT reference.

The pipeline finds all actual:

```text
Glyc*.pdb
```

files in `pos_0` and matches them to entries in:

```text
pos_0/Glyc_score.sc
```

Only matched models are used.

The WT baseline is:

```text
WT_mean = mean(total_score of matched pos_0 models)
```

This prevents stale score-file entries for missing PDBs from affecting the
calculation.

---

## 2.3 Step 2 — Evaluate every glycan-masking position

For every `pos_i`, the same procedure is performed.

For example:

```text
pos_581/
├── Glyc_score.sc
├── Glyc_..._0001.pdb
├── Glyc_..._0002.pdb
├── ...
└── Glyc_..._0137.pdb
```

All 137 valid models are used.

The mean score is:

```text
mean_i = mean(total_score of matched models in pos_i)
```

The energy difference relative to WT is:

```text
delta_REU_i = mean_i - WT_mean
```

The default criterion is:

```text
delta_REU_i < 0
```

which is equivalent to:

```text
mean_i < WT_mean
```

The threshold is configurable:

```bash
--reu-delta -5
```

would require:

```text
delta_REU_i < -5
```

---

## 2.4 Step 3 — Handle positions without Rosetta glycan models

Some deep glycan masking positions may fail to produce acceptable glycosylated
Rosetta structures.

For example:

```text
pos_57/
```

may contain no files beginning with `Glyc`.

The pipeline records:

```text
n_pdb = 0
status = NO_GLYCAN_MODELS
selected = false
```

and skips the position.

No fixed number of expected models is assumed.

---

## 2.5 Step 4 — Extract all N-linked glycans

For each PDB from a selected position, the pipeline examines PDB `LINK`
records.

An N-linked glycan attachment is identified through connectivity equivalent to:

```text
ASN ND2  →  sugar C1
```

The code does not rely on a single carbohydrate residue name such as `NAG`.

This is important for Rosetta output because carbohydrate naming can differ,
for example `Glc`, `Man`, or other Rosetta-compatible names.

Once the ASN-to-sugar attachment is found, the carbohydrate `LINK` graph is
followed to recover the complete glycan tree.

---

## 2.6 Step 5 — Collect engineered AND native glycans

Suppose a model from:

```text
pos_581
```

contains:

```text
ASN 141 → native glycan
ASN 200 → native glycan
ASN 581 → engineered glycan
```

all three glycan conformations are useful Rosetta-relaxed structures and are
therefore extracted.

Their provenance remains separate:

```text
A:141  native
A:200  native
A:581  engineered
```

This increases the amount of physically sampled conformational information
available to the empirical library.

---

## 2.7 Step 6 — Normalize every glycan to a common ASN frame

Glycans originate from many different positions on the protein.

Their raw Cartesian coordinates therefore cannot simply be combined.

For each glycan, the local attachment ASN provides a coordinate frame:

```text
CB
 \
  CG
   \
    ND2 ─ C1 ─ glycan
```

The three ASN atoms:

```text
CB
CG
ND2
```

are Kabsch-aligned to one canonical ASN attachment frame.

The same rigid-body transformation is applied to the complete glycan.

Conceptually:

```text
glycan at ASN 1   ─┐
glycan at ASN 141 ─┼── align ASN CB/CG/ND2 ──→ common frame
glycan at ASN 200 ─┤
glycan at ASN 581 ─┘
```

This removes differences caused simply by where the glycan was located on the
protein.

Importantly, the glycan's Rosetta-generated orientation **relative to its own
ASN attachment site is preserved**.

Aligning only one N atom would not be sufficient because one point defines
translation but not orientation. The `CB/CG/ND2` frame defines both.

---

## 2.8 Step 7 — Validate common glycan topology

All conformers stored in one coordinate array must correspond atom-for-atom.

The first successfully extracted glycan establishes the canonical topology and
atom signature.

Subsequent glycans must have compatible:

```text
glycan topology
residue identities
atom identities
atom correspondence
```

The carbohydrate tree is placed into deterministic root-first ordering before
comparison.

A glycan with incompatible topology is skipped and reported in:

```text
*_skipped.tsv
```

Use:

```bash
--strict
```

to stop immediately instead.

---

## 2.9 Step 8 — Build the empirical library

After normalization, the extracted coordinates are stored in:

```text
Man5_empirical.npz
```

These conformers are directly derived from the existing Rosetta models.

They are labeled:

```text
conformer_origin = rosetta_empirical
```

This is the primary output of the pipeline.

For a large deep glycan scanning dataset, this empirical library may already
contain thousands or tens of thousands of conformations.

---

## 2.10 Step 9 — Optional conformer expansion

Additional conformers can be generated from each empirical Rosetta glycan.

The expansion step samples rotations around:

1. the ASN ND2–root sugar C1 attachment axis;
2. carbohydrate glycosidic linkage axes.

The goal is to introduce additional torsional diversity without deliberately
changing normal covalent bond lengths or bond angles.

For example:

```text
20,000 empirical glycans

--samples-per-glycan 4

20,000 original
+
80,000 generated
=
100,000 total conformers
```

Generated structures are labeled:

```text
conformer_origin = torsion_expanded
```

and retain a link to their empirical parent through:

```text
parent_empirical_index
```

---

## 2.11 Recommended workflow

For a new dataset, use the pipeline in stages.

### Stage A — Analyze

```bash
MODE="analyze"
bash run_existing_models_library.sh
```

Inspect which positions pass the Rosetta energy criterion.

### Stage B — Build empirical library

```bash
MODE="empirical"
bash run_existing_models_library.sh
```

Inspect the number and diversity of Rosetta-derived conformers.

### Stage C — Expand only if needed

```bash
MODE="expanded"
bash run_existing_models_library.sh
```

Additional sampling is useful only if the empirical Rosetta library does not
already provide sufficient conformational coverage.

A large number of conformers is not automatically better. Because many
Rosetta models may be geometrically similar, a later diversity-analysis or
clustering step can be used to measure redundancy and select representative
conformations before increasing library size further.
