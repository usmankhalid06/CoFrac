CoFrac: Cross-objective blind source separation for interpreting ESM-2 representations

Code and data-generation scripts for the paper
"Cross-objective blind source separation reveals functional sub-features in ESM-2 protein language model representations."

CoFrac frames the per-residue activations of the ESM-2 protein language model as a blind source separation (BSS) problem and tests whether a recovered latent partition belongs to the representation rather than to any single decomposition objective. The pipeline extracts activations in Python, runs six decompositions plus a decomposition-free control in MATLAB, and renders the structural figures back in Python.


Important: the activation file is not included

The main data file, esm2_35M_multilayer.mat (ESM-2-35M activations, ~4 GB, MATLAB v7.3/HDF5), is too large to host on GitHub and is not part of this repository. It must be regenerated locally in one step with the extractor script (details below). Regeneration needs internet access to huggingface.co (to download the ESM-2 weights) and rest.uniprot.org (to download the protein set and annotations).

The exact 297-protein set used in the paper is pinned in protein_metadata_12.tsv (accession, name, family, function, and full sequence for every protein), so the input set can be reproduced independently of any later UniProt updates.


Repository contents

Step 1 — Data extraction (Python)

FilePurposescript_Extract_esm2_activations_multilayer_12.pyRegenerates esm2_35M_multilayer.mat. Downloads ESM-2-35M, pulls the reviewed human Swiss-Prot slice (organism_id:9606), runs each sequence through the model once, and saves per-residue activations for all 12 layers plus residue-level labels (active site, binding site, transmembrane, disulfide bond, DNA binding).script_Download_protein_metadata_12.pyFetches metadata (name, family, function, sequence) for the 297 accessions; writes protein_metadata_12.tsv. Useful for the family-redundancy check discussed in the paper.protein_metadata_12.tsvThe pinned 297-protein set (accession + sequence + annotations).

Step 2 — Decompositions and analysis (MATLAB)

Method implementations:

FileMethodmy_SAE.mTopK sparse autoencoder (ℓ₀)my_ACSD.mAdaptive Consistent Sequential Dictionary Learning (ℓ₁)LSICA.mSparse orthogonal component analysissgBACES.mSparse group BACES (cross-layer, common + layer-specific)gICA.mGroup ICA (cross-layer independence)

Driver scripts (operating layer 9, K = 200):

FileRunsscript1_SAE_only_SingleLayer_9.mSAE + census, cysteine partition, chance-partition nullscript2_LSICA_only_SingleLayer_9.mLSICAscript3_ACSD_only_SingleLayer_9.mACSDscript4_ICA_SingleLayer_9.mFastICA (spatial)script5_sgBACES_only_SingleLayer_9.msgBACESscript6_gICA_only_SingleLayer_9.mGroup ICA

Each driver reads esm2_35M_multilayer.mat, fits its decomposition, names atoms by amino-acid identity and Swiss-Prot function, and writes the per-atom activation CSVs (atom*_SAE.csv, atom*_ICA.csv, atom*_LSICA.csv, …) used by the plotting scripts.

Step 3 — Structural figures (Python)

FilePurposescript_protein_Plot_SAE.pyMaps the SAE transmembrane / disulfide / DNA-binding atoms onto AlphaFold structures (fetched from the EBI AlphaFold API), firing residues in orange, cysteines as sticks in the disulfide row. Writes one struct_atom*_SAE.html per atom.script_protein_Plot_ICA.pySame for ICA atoms.script_protein_Plot_LSICA.pySame for LSICA atoms.script_convert_html_png_12.pyRenders the struct_atom*.html viewers to PNG (headless Chrome / Selenium).


Requirements

Python (extraction and plotting)

pip install torch transformers numpy requests hdf5storage pandas py3Dmol selenium

script_convert_html_png_12.py also needs Chrome + a matching ChromeDriver on the PATH.

MATLAB (analysis) — R2021b or newer recommended (for matfile partial loading of the v7.3 .mat).


How to reproduce, end to end

bash# 1. Regenerate the activation file (~4 GB, needs internet)
python script_Extract_esm2_activations_multilayer_12.py
#    -> esm2_35M_multilayer.mat

matlab% 2. In MATLAB, run each decomposition at the operating layer
script1_SAE_only_SingleLayer_9
script2_LSICA_only_SingleLayer_9
script3_ACSD_only_SingleLayer_9
script4_ICA_SingleLayer_9
script5_sgBACES_only_SingleLayer_9
script6_gICA_only_SingleLayer_9
%    -> per-atom CSVs + console tables (Tables 3 and 4 of the paper)

bash# 3. Render the structural figures
python script_protein_Plot_SAE.py
python script_protein_Plot_ICA.py
python script_protein_Plot_LSICA.py
python script_convert_html_png_12.py
#    -> struct_atom*.html and .png


Reproducibility notes


Data set. The extractor queries UniProt with (reviewed:true) AND (organism_id:9606) and stops at the residue cap, which yields the same 297 human proteins (150,283 residues) used in the paper. Because it queries the live UniProt database, a future UniProt release could shift the set slightly; protein_metadata_12.tsv pins the exact accessions and sequences used here for anyone who needs a byte-exact reproduction.
Larger-model appendix. The ESM-2-150M replication (Appendix A) uses the same protein set re-encoded with a larger model. It is produced by the 150M variant of the extractor (MODEL_NAME = facebook/esm2_t30_150M_UR50D, 30 layers, layer 18); add that script to the repo if you want the appendix reproduced.
ICA. FastICA is run from a single initialization; individual enrichment magnitudes can vary across restarts, while the qualitative two-pole cysteine partition is stable (see Appendix A.2).
Family redundancy. The protein set is not controlled for protein-family redundancy (see Limitations); protein_metadata_12.tsv includes the protein_families field for anyone wishing to quantify or control for it.


Citation

If you use this code, please cite the paper (see the repository release / CITATION once available).
