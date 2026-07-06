import sys
import numpy as np
import requests
import hdf5storage                       # writes MATLAB v7.3 (HDF5) files

# ============================ CONFIG =======================================
MODEL_NAME   = "facebook/esm2_t12_35M_UR50D"  # 12 layers, dim 480
LAYERS       = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]   # transformer layers to extract (1..12)
N_PROTEINS   = 400        # how many Swiss-Prot proteins to pull (<= 500 per request)
MAX_RESIDUES = 150000     # cap total rows; v7.3 handles big files (180k x 480 x 12 ~ 4.1 GB)
ORGANISM_ID  = 9606       # 9606 = human; None for all organisms
MAX_LEN      = 1022       # ESM-2 sequence-length cap (longer sequences truncated)
OUTPUT       = "esm2_35M_multilayer.mat"
SAVE_DTYPE   = np.float32

FEATURE_TYPES = ["Active site", "Binding site", "Transmembrane",
                 "Disulfide bond", "DNA binding"]
# ===========================================================================


def fetch_swissprot(n, organism_id):
    query = "(reviewed:true)"
    if organism_id:
        query += f" AND (organism_id:{organism_id})"
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": query, "format": "json", "size": min(n, 500)}
    r = requests.get(url, params=params,
                     headers={"User-Agent": "dl-diagnostic/1.0"}, timeout=120)
    r.raise_for_status()
    return r.json().get("results", [])[:n]


def residue_labels_for(features, seq_len, type_to_col):
    K = len(type_to_col)
    lab = np.zeros((seq_len, K), dtype=np.uint8)
    for f in features:
        ftype = (f.get("type") or "").strip().lower()
        col = type_to_col.get(ftype)
        if col is None:
            continue
        loc = f.get("location", {}) or {}
        s = (loc.get("start", {}) or {}).get("value")
        e = (loc.get("end", {}) or {}).get("value")
        if s is None or e is None:
            continue
        s, e = int(s), int(e)
        if "disulfide" in ftype:
            positions = [p for p in (s, e) if 1 <= p <= seq_len]
        else:
            positions = range(max(1, s), min(seq_len, e) + 1)
        for p in positions:
            lab[p - 1, col] = 1
    return lab


def main():
    type_to_col = {t.lower(): i for i, t in enumerate(FEATURE_TYPES)}

    print(f"Fetching up to {N_PROTEINS} Swiss-Prot proteins "
          f"(organism_id={ORGANISM_ID}) ...")
    entries = fetch_swissprot(N_PROTEINS, ORGANISM_ID)
    print(f"  got {len(entries)} entries")

    print(f"Loading {MODEL_NAME} ...")
    import torch
    from transformers import AutoTokenizer, AutoModel
    tok = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = AutoModel.from_pretrained(MODEL_NAME, output_hidden_states=True)
    model.eval()
    n_layers = model.config.num_hidden_layers
    bad = [l for l in LAYERS if not (1 <= l <= n_layers)]
    if bad:
        sys.exit(f"LAYERS {bad} out of range; model has {n_layers} layers.")
    print(f"  model has {n_layers} layers; extracting {LAYERS}")

    X_parts, idx_parts, pos_parts, lab_parts, accessions = [], [], [], [], []
    total = 0

    for entry in entries:
        acc = entry.get("primaryAccession")
        seq = (entry.get("sequence", {}) or {}).get("value")
        if not acc or not seq:
            continue
        seq = seq[:MAX_LEN]
        L = len(seq)

        enc = tok(seq, return_tensors="pt", truncation=True, max_length=MAX_LEN + 2)
        with torch.no_grad():
            out = model(**enc)
        hs = out.hidden_states           # tuple length (n_layers+1); hs[l] is layer l

        # stack the requested layers into (L_residues, dim, n_selected_layers)
        per_layer = []
        for l in LAYERS:
            resid = hs[l][0][1:1 + L].cpu().numpy()   # (L, dim), strip <cls>/<eos>
            per_layer.append(resid.astype(SAVE_DTYPE))
        block = np.stack(per_layer, axis=2)           # (L, dim, n_layers)
        Lr = block.shape[0]                           # actual residue count

        pidx = len(accessions) + 1
        accessions.append(acc)
        X_parts.append(block)
        idx_parts.append(np.full((Lr, 1), pidx, dtype=np.int32))
        pos_parts.append(np.arange(1, Lr + 1, dtype=np.int32).reshape(-1, 1))
        lab_parts.append(residue_labels_for(entry.get("features", []), Lr, type_to_col))

        total += Lr
        if len(accessions) % 25 == 0:
            print(f"  {len(accessions)} proteins, {total} residues")
        if total >= MAX_RESIDUES:
            print(f"  reached MAX_RESIDUES={MAX_RESIDUES}; stopping")
            break

    if not X_parts:
        sys.exit("No activations collected -- check network / query.")

    X = np.ascontiguousarray(np.vstack(X_parts))     # (N, dim, n_layers), single
    protein_idx = np.vstack(idx_parts)
    residue_pos = np.vstack(pos_parts)
    labels = np.vstack(lab_parts)
    N, d, Ln = X.shape

    # ---- SAVE as MATLAB v7.3 (HDF5): no 2 GB cap, partial-loadable by layer ----
    print(f"\nSaving {OUTPUT}  (MATLAB v7.3 / HDF5) ...")
    mdict = {
        "X": X,
        "layers_used": np.array(LAYERS, dtype=np.int32).reshape(1, -1),
        "protein_idx": protein_idx,
        "residue_pos": residue_pos,
        "labels": labels,
        "label_names": np.array(FEATURE_TYPES, dtype=object),
        "accessions": np.array(accessions, dtype=object),
    }
    hdf5storage.savemat(OUTPUT, mdict, format="7.3",
                        matlab_compatible=True, store_python_metadata=False)

    print(f"  X: {N} residues x {d} dims x {Ln} layers  ({X.dtype})")
    print(f"  layers_used: {LAYERS}")
    print(f"  proteins: {len(accessions)}")
    print("  label coverage (residues positive per concept):")
    for i, name in enumerate(FEATURE_TYPES):
        print(f"    {name:<16} {int(labels[:, i].sum())}")
    approx_gb = N * d * Ln * 4 / 1e9
    print(f"  (X size ~ {approx_gb:.2f} GB; v7.3 handles it, MATLAB reads per-layer)")
    print("\nIn MATLAB:  mf=matfile('%s'); Xl=double(mf.X(:,:,k));" % OUTPUT)


if __name__ == "__main__":
    main()