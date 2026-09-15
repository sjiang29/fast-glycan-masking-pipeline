# fast_glycan_masking_from_existing_models

Build a reusable N-glycan conformer library from thousands of already-relaxed
Rosetta glycoprotein models.

## Expected layout

```text
ROOT/
  pos_0/
    Glyc_score.sc
    Glyc_...._0001.pdb
    ...
  pos_1/
    Glyc_score.sc
    Glyc_...._0001.pdb
    ...
  pos_2/
    # zero Glyc*.pdb is allowed; this position is skipped
  ...
```

`pos_0` is the WT energy baseline. `pos_N` means the engineered N-glycan is at
protein sequence/PDB residue index N. All `Glyc*.pdb` models present are used;
there is no assumption of exactly 100 models.

Only score rows whose `description` matches an existing `Glyc*.pdb` are used.
The Rosetta column used for filtering is `total_score`.

## 1. Analyze positions

```bash
python -m fast_glycan_masking_from_existing_models analyze \
  --input-dir /path/to/ROOT \
  --reu-delta 0 \
  --out position_scores.csv
```

A position is selected when:

`mean_total_score(pos_i) - mean_total_score(pos_0) < reu_delta`

Default `reu_delta=0`.

## 2. Build empirical library

```bash
python -m fast_glycan_masking_from_existing_models build \
  --input-dir /path/to/ROOT \
  --out-prefix libraries/Man5_empirical \
  --native-positions 200
```

All ASN-linked glycans in each selected Rosetta model are extracted, including
the engineered glycan and native glycans. Each glycan is independently aligned
to one canonical ASN attachment frame using ASN CB/CG/ND2. This removes global
protein orientation while preserving the Rosetta glycan orientation relative
to its attachment ASN.

`pos_0` defines the WT score baseline and is excluded from the glycan library by
default. Add `--include-wt-glycans` to include its native glycans.

## 3. Optional expansion

```bash
python -m fast_glycan_masking_from_existing_models build \
  --input-dir /path/to/ROOT \
  --out-prefix libraries/Man5_expanded \
  --samples-per-glycan 4 \
  --attachment-sigma-deg 20 \
  --glycosidic-sigma-deg 20
```

The empirical Rosetta conformer is always retained. Additional conformers are
created by rotations around the ASN-ND2/root-C1 attachment bond and glycosidic
LINK bond axes. Bond lengths and bond angles are not deliberately distorted.

For a large empirical library, first inspect its diversity before expanding it;
thousands of Rosetta models may already provide enough conformational coverage.

## Outputs

For `--out-prefix libraries/Man5_empirical`:

- `Man5_empirical.npz` — reusable coordinate library + provenance
- `Man5_empirical_positions.csv` — position-level REU filtering
- `Man5_empirical_provenance.csv` — source PDB/site/score for every conformer
- `Man5_empirical_manifest.json` — run parameters and summary
- `Man5_empirical_skipped.tsv` — invalid/mismatched models, if any
- `Man5_empirical.pdb` — optional with `--write-pdb`

## Notes

The extractor follows PDB `LINK` records beginning with `ASN ND2 -> sugar C1`.
This intentionally avoids hard-coding Rosetta carbohydrate residue names such
as `Glc`/`Man`.

All glycans in one library must have the same topology and atom signature.
Models with a different glycan topology are skipped unless `--strict` is used.
