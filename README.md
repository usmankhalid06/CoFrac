# CoFrac

**Cross-objective blind source separation for interpreting ESM-2 protein language model representations.**

CoFrac treats the per-residue activations of ESM-2 as a blind source separation (BSS) problem and asks whether a recovered latent feature belongs to the *representation* or is an artifact of one decomposition. It runs six decompositions (sparse, independence-based, cross-layer) plus a decomposition-free control, and validates the recovered structure against held-out Swiss-Prot annotations — with the disulfide / non-disulfide cysteine split as the functional anchor.

Code accompanying the paper *"Cross-objective blind source separation reveals functional sub-features in ESM-2 protein language model representations."*

---

## Quick start

```bash
# 1. Regenerate the activation file (~4 GB, needs internet)
python script_Extract_esm2_activations_multilayer_12.py      # -> esm2_35M_multilayer.mat
```
```matlab
% 2. Run the six decompositions at the operating layer (MATLAB)
script1_SAE_only_SingleLayer_9
script2_LSICA_only_SingleLayer_9
script3_ACSD_only_SingleLayer_9
script4_ICA_SingleLayer_9
script5_sgBACES_only_SingleLayer_9
script6_gICA_only_SingleLayer_9                              % -> per-atom CSVs + result tables
```
```bash
# 3. Render the structural figures
python script_protein_Plot_SAE.py
python script_protein_Plot_ICA.py
python script_protein_Plot_LSICA.py
python script_convert_html_png_12.py                          # -> struct_atom*.html / .png
```

---

## The activation file is not in this repo

The main data file **`esm2_35M_multilayer.mat`** (ESM-2-35M activations, ~4 GB, MATLAB v7.3/HDF5) is **too large for GitHub** and is not included. Regenerate it in one step with `script_Extract_esm2_activations_multilayer_12.py` (Step 1 above). This needs internet access to `huggingface.co` (ESM-2 weights) and `rest.uniprot.org` (protein set + annotations).

The exact 297-protein set is pinned in **`protein_metadata_12.tsv`** (accession, name, family, function, and full sequence for every protein), so the input can be reproduced independently of later UniProt updates.

---

## What's in here

### Step 1 — Data extraction (Python)

| File | Purpose |
|------|---------|
| `script_Extract_esm2_activations_multilayer_12.py` | Regenerates `esm2_35M_multilayer.mat`: downloads ESM-2-35M, pulls the reviewed human Swiss-Prot slice (`organism_id:9606`), runs each sequence through the model once, and saves per-residue activations for all 12 layers plus residue-level labels (active site, binding site, transmembrane, disulfide bond, DNA binding). |
| `script_Download_protein_metadata_12.py` | Fetches metadata for the 297 accessions; writes `protein_metadata_12.tsv`. |
| `protein_metadata_12.tsv` | The pinned 297-protein set (accession + sequence + annotations). |

### Step 2 — Decompositions and analysis (MATLAB)

Method implementations:

| File | Method |
|------|--------|
| `my_SAE.m` | TopK sparse autoencoder (ℓ₀) |
| `my_ACSD.m` | Adaptive Consistent Sequential Dictionary Learning (ℓ₁) |
| `LSICA.m` | Sparse orthogonal component analysis |
| `sgBACES.m` | Sparse group BACES (cross-layer, common + layer-specific) |
| `gICA.m` | Group ICA (cross-layer independence) |

> **ICA dependency.** There is no `ICA.m` in this repo. The ICA driver (`script4_ICA_SingleLayer_9.m`) and `gICA.m` rely on the external **FastICA** package for MATLAB (`fastica` / `fpica`). Download it from [https://research.ics.aalto.fi/ica/fastica/](https://research.ics.aalto.fi/ica/fastica/) and add it to your MATLAB path before running Step 2. We use the symmetric approach with the `tanh` contrast function.

Driver scripts (operating layer 9, K = 200):

| File | Runs |
|------|------|
| `script1_SAE_only_SingleLayer_9.m` | SAE + census, cysteine partition, chance-partition null |
| `script2_LSICA_only_SingleLayer_9.m` | LSICA |
| `script3_ACSD_only_SingleLayer_9.m` | ACSD |
| `script4_ICA_SingleLayer_9.m` | FastICA (spatial) |
| `script5_sgBACES_only_SingleLayer_9.m` | sgBACES |
| `script6_gICA_only_SingleLayer_9.m` | Group ICA |

Each driver reads `esm2_35M_multilayer.mat`, fits its decomposition, names atoms by amino-acid identity and Swiss-Prot function, and writes the per-atom activation CSVs (`atom*_SAE.csv`, `atom*_ICA.csv`, `atom*_LSICA.csv`, …) used by the plotting scripts.

### Step 3 — Structural figures (Python)

| File | Purpose |
|------|---------|
| `script_protein_Plot_SAE.py` | Maps the SAE transmembrane / disulfide / DNA-binding atoms onto AlphaFold structures (EBI AlphaFold API); firing residues in orange, cysteines as sticks in the disulfide row. Writes one `struct_atom*_SAE.html` per atom. |
| `script_protein_Plot_ICA.py` | Same, for ICA atoms. |
| `script_protein_Plot_LSICA.py` | Same, for LSICA atoms. |
| `script_convert_html_png_12.py` | Renders the `struct_atom*.html` viewers to PNG (headless Chrome / Selenium). |

---

## Requirements

**Python** (extraction and plotting)
```
pip install torch transformers numpy requests hdf5storage pandas py3Dmol selenium
```
`script_convert_html_png_12.py` also needs Chrome + a matching ChromeDriver on the PATH.

**MATLAB** (analysis) — R2021b or newer recommended (for `matfile` partial loading of the v7.3 `.mat`). The ICA and group-ICA drivers additionally require the external **FastICA** package ([download](https://research.ics.aalto.fi/ica/fastica/)) on the MATLAB path.

---

## Reproducibility notes

- **Data set.** The extractor queries UniProt with `(reviewed:true) AND (organism_id:9606)` and stops at the residue cap, yielding the same 297 human proteins (150,283 residues) used in the paper. Because it queries the live UniProt database, a future release could shift the set slightly; `protein_metadata_12.tsv` pins the exact accessions and sequences for a byte-exact reproduction.
- **Larger-model appendix.** The ESM-2-150M replication uses the same protein set re-encoded with a larger model (`MODEL_NAME = facebook/esm2_t30_150M_UR50D`, 30 layers, layer 18). Add that extractor variant to reproduce the appendix.
- **ICA.** FastICA is run from a single initialization; individual enrichment magnitudes can vary across restarts, while the qualitative two-pole cysteine partition is stable.
- **Family redundancy.** The protein set is not controlled for protein-family redundancy; `protein_metadata_12.tsv` includes the `protein_families` field for anyone wishing to quantify or control for it.

---

## Citation

A citation entry will be added here once the paper is published.
