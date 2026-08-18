# Domain Knowledge for Fast Man5 Glycan Masking

This document summarizes the scientific concepts behind the fast Man5 conformer-library approach.

For installation and command-line usage, see `README.md`.

---

# 1. Why Use a Glycan Conformer Library?

Glycans are highly flexible molecules.

A Man5 glycan attached to a protein can adopt many conformations because several glycosidic bonds contain rotatable torsions.

Traditional Rosetta glycan modeling generates and optimizes these conformations separately for each protein and glycosylation site.

For glycan masking studies, the main goal is different:

> Generate a sufficiently diverse set of physically reasonable glycan conformations to represent the steric space accessible to the glycan during antibody docking.

The conformer-library approach therefore generates glycan conformations once and reuses them across different protein sites.

---

# 2. N-Linked Glycosylation Sequon

Canonical N-linked glycosylation occurs at:

```text
N-X-S/T
```

where:

```text
N = Asparagine
X = any amino acid except Proline
S/T = Serine or Threonine
```

The glycan is attached to the side-chain nitrogen `ND2` of Asparagine.

The first carbohydrate residue is GlcNAc.

Conceptually:

```text
Protein — ASN — GlcNAc — GlcNAc — Man — ...
```

---

# 3. Two Levels of Glycan Flexibility

The current library samples two major levels of conformational flexibility.

## ASN → First GlcNAc Orientation

The protein–glycan interface is represented using the torsion:

```text
CG(ASN) – ND2(ASN) – C1(GlcNAc) – O5(GlcNAc)
```

Rotation occurs around the:

```text
ND2 – C1
```

bond.

This changes the orientation of the entire glycan relative to the ASN.

Conceptually:

```text
              glycan
                |
ASN-CG-ND2-----C1
          ↻
```

The previous implementation transferred the root glycan orientation directly from the reference structure.

That can be restrictive because a root orientation that works at one protein position may point toward the protein surface at another position.

Sampling the attachment torsion allows alternative orientations to be tested.

---

# 4. What Is Kept Fixed at the Attachment?

The pipeline does not randomly distort the covalent ASN–GlcNAc connection.

The following remain fixed during attachment-torsion sampling:

```text
ND2–C1 bond length
CG–ND2–C1 bond angle
ND2–C1–O5 bond angle
```

The main variable is the dihedral:

```text
CG–ND2–C1–O5
```

This provides rotational flexibility without arbitrarily changing basic covalent geometry.

---

# 5. Glycosidic φ and ψ Angles

For sugar–sugar glycosidic bonds, conformational flexibility is primarily described using the glycosidic torsions:

```text
φ
ψ
```

Different linkage types have different preferred conformational regions.

Examples include:

```text
β1→4
α1→3
α1→6
```

The preferred φ/ψ distributions are therefore linkage-dependent rather than identical for every glycosidic bond.

The pipeline determines the linkage type from the glycan connectivity and samples the corresponding distribution defined in:

```text
configs/man5_sampling.yaml
```

---

# 6. ω Angle for 1→6 Linkages

A `1→6` linkage contains an additional rotatable bond compared with many `1→3` or `1→4` linkages.

Therefore, it requires an additional torsion:

```text
ω
```

The library samples:

```text
φ
ψ
ω
```

for these linkages.

For other supported linkages, φ and ψ are sampled without ω.

---

# 7. Why Not Sample All Angles Uniformly?

Uniformly rotating every glycosidic bond over the entire `0–360°` range would generate many conformations that are energetically unfavorable or rarely observed.

Instead, the pipeline supports sampling around linkage-specific preferred conformations.

A typical configuration may use:

```yaml
default:
  phi:
    mode: template_jitter
    sd: 15.0

  psi:
    mode: template_jitter
    sd: 15.0
```

More specific distributions can be defined for individual linkage types.

This allows the conformer library to emphasize plausible conformational regions without performing full Rosetta energy optimization.

---

# 8. Template-Jitter Sampling

A simple sampling strategy is:

```yaml
mode: template_jitter
sd: 15.0
```

If the torsion in the reference glycan is:

```text
θtemplate
```

the new conformations are sampled around that value.

Conceptually:

```text
θsample ≈ θtemplate ± variation
```

The same strategy can be used for the ASN→GlcNAc attachment torsion:

```yaml
attachment:
  mode: template_jitter
  sd: 20.0
```

This retains information from a physically reasonable reference glycan while introducing conformational diversity.

---

# 9. Why Attachment Sampling Is Important

Consider two surface positions.

At one site, the reference attachment orientation may point into solvent:

```text
Protein — ASN
              \
               GlcNAc — Man5
```

The glycan can be placed successfully.

At another site, the same relative orientation may point toward the protein:

```text
              protein surface
                 ███████
Protein — ASN — GlcNAc██
```

Changing only downstream Man5 φ/ψ/ω angles may not solve a collision involving the root GlcNAc.

Attachment-torsion sampling allows the whole glycan to rotate:

```text
              GlcNAc — Man5
             /
Protein — ASN
       █████ protein
```

This increases the probability of finding a clash-free root orientation.

---

# 10. ASN Alignment During Placement

The reusable glycan library is generated using one reference ASN.

The library stores the reference ASN frame using:

```text
CB
CG
ND2
```

For a new glycosylation position, the placer obtains the same atoms from the target ASN.

It then aligns:

```text
source ASN CB/CG/ND2
          ↓
target ASN CB/CG/ND2
```

and applies the corresponding rigid-body transformation to the glycan coordinates.

Because different conformers already contain different attachment and downstream torsion values, these differences are preserved after placement.

---

# 11. Protein–Glycan Clash Filtering

Not every physically reasonable free-glycan conformation will fit at every protein position.

After placement, each conformer is therefore checked against the protein.

A conformer is rejected when glycan atoms overlap too strongly with protein atoms.

The recommended mode uses scaled van der Waals radii:

```text
threshold = scale × (r1 + r2)
```

where:

```text
r1 = van der Waals radius of glycan atom
r2 = van der Waals radius of protein atom
```

A smaller scale is more permissive.

For example:

```text
0.70
```

permits closer contacts than:

```text
0.80
```

The clash filter is intended as a fast steric screen rather than a full molecular-mechanics energy calculation.

---

# 12. Glycan–Glycan Clash Filtering

When multiple glycans are placed simultaneously, individually acceptable conformers can still collide with one another.

For example:

```text
         glycan A
            \
Protein ------\
               X
Protein ------/
            /
         glycan B
```

The multi-site placement procedure therefore performs two filters:

```text
protein ↔ glycan
glycan ↔ glycan
```

Only combinations passing both checks are written.

---

# 13. Native and Engineered Glycosylation Sites

A protein may already contain native glycosylation sites.

If a new masking site is added, the desired model may contain:

```text
native site 1
+
native site 2
+
engineered site
```

These should be considered simultaneously because neighboring glycans can influence the steric environment available to each other.

The multi-site placement implementation therefore chooses an independent library conformer for every requested site and evaluates the complete combination.

---

# 14. Role of Rosetta FastDesign

The conformer library is intended to replace the expensive repeated glycan ensemble-generation step, not all Rosetta modeling.

Rosetta FastDesign still has an important role when introducing new masking sites.

It can:

```text
create the N-X-S/T sequon
build the ASN side chain
optimize the local protein environment
```

The responsibilities are therefore separated:

```text
FastDesign
    │
    ├── protein sequence
    └── local protein/ASN geometry
          │
          ▼
Man5 conformer library
    │
    ├── ASN→GlcNAc orientation
    └── glycan φ/ψ/ω
          │
          ▼
Clash filtering
```

---

# 15. What the Pipeline Is Designed to Represent

This pipeline does not attempt to identify one exact minimum-energy glycan structure.

Instead, it generates an **ensemble**.

For antibody docking, the important question is often:

> What regions of the antigen surface can reasonably be occupied or sterically masked by the glycan?

Therefore, conformational diversity can be more useful than selecting a single optimized glycan structure.

---

# 16. Interpretation of a Failed Site

If:

```text
0 / 10000
```

library conformers survive at a target site, several explanations are possible:

- the site is buried,
- the local ASN orientation is unfavorable,
- the root GlcNAc cannot fit,
- attachment sampling is too narrow,
- downstream glycan sampling is insufficient,
- the clash criterion is too restrictive,
- local protein remodeling may be required.

A failed site therefore does not automatically mean there is an error in the code.

It may indicate that the requested location is geometrically unfavorable for Man5 attachment under the sampled conformational assumptions.

---

# 17. Scope of the Current Model

The current implementation is intended primarily for:

```text
N-linked Man5
+
glycan masking
+
antibody docking
```

It should not be interpreted as a replacement for:

- molecular dynamics,
- free-energy calculations,
- full glycan energy minimization,
- detailed prediction of glycan population distributions.

Its purpose is efficient generation of structurally plausible glycosylated antigen ensembles for downstream structural analysis.