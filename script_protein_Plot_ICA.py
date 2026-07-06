import requests, pandas as pd, os
import py3Dmol
ATOMS = {22: "Transmembrane", 48: "Disulfide bond", 150: "DNA binding"}
PREFERRED = {22: None, 48: None, 150: None}
def load_af_for_atom(atom):
    df = pd.read_csv(f"atom{atom}_ICA.csv")
    ranked = df.groupby("accession")["activation"].sum().sort_values(ascending=False).index.tolist()
    pref = PREFERRED.get(atom)
    order = ([pref] + [a for a in ranked if a != pref]) if pref else ranked
    for cand in order:
        local = f"AF-{cand}-model.pdb"
        if os.path.exists(local) and "ATOM" in open(local).read():
            return cand, local, df
        try:
            meta = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{cand}", timeout=30)
        except Exception as e:
            print(f"  skip {cand}: request error {e}"); continue
        if meta.status_code != 200:
            print(f"  skip {cand}: API {meta.status_code}"); continue
        js = meta.json()
        if not js:
            print(f"  skip {cand}: no AF entry"); continue
        pdb_url = js[0].get("pdbUrl")
        if not pdb_url:
            print(f"  skip {cand}: no pdbUrl"); continue
        r = requests.get(pdb_url, timeout=60)
        if r.status_code == 200 and "ATOM" in r.text:
            open(local, "w").write(r.text); return cand, local, df
        print(f"  skip {cand}: download {r.status_code}")
    return None, None, None
def firing_residues(pdbfile, df, acc, thr=0.5):
    ca = [l for l in open(pdbfile) if l.startswith("ATOM") and l[12:16].strip()=="CA"]
    resnum = [int(l[22:26]) for l in ca]
    resname = {int(l[22:26]): l[17:20].strip() for l in ca}
    sub = df[df.accession==acc]; act = dict(zip(sub.resnum.astype(int), sub.activation.astype(float)))
    mx = max(act.values()) if len(act) else 1
    fire = [r for r in resnum if act.get(r,0.0)/(mx or 1) >= thr]
    return fire, [r for r in fire if resname.get(r)=="CYS"]
def make_view(atom, concept):
    acc, pdbfile, df = load_af_for_atom(atom)
    if acc is None:
        print(f"atom {atom} ({concept}): no model"); return None, acc
    print(f"atom {atom} ({concept}) -> {acc}")
    fire, fire_cys = firing_residues(pdbfile, df, acc)
    view = py3Dmol.view(width=700, height=700)
    view.addModel(open(pdbfile).read(), "pdb")
    view.setStyle({"cartoon": {"color": "lightgrey"}})
    if fire:
        view.addStyle({"resi":[str(r) for r in fire]}, {"cartoon":{"color":"orange"}})
    if "Disulfide" in concept and fire_cys:
        view.addStyle({"resi":[str(r) for r in fire_cys]},
                      {"stick":{"colorscheme":"yellowCarbon","radius":0.3}})
    view.zoomTo(); view.setBackgroundColor("white")
    return view, acc
# --- one HTML per atom (no titles, no toolbar) ---
for atom, concept in ATOMS.items():
    view, acc = make_view(atom, concept)
    if view is None:
        continue
    inner = view._make_html()
    html = (f'<html><head><meta charset="utf-8">'
            f'<script src="https://3dmol.org/build/3Dmol-min.js"></script>'
            f'<style>body{{margin:0;overflow:hidden}} '
            f'.viewer_3Dmoljs_menu,.mol-container .icon{{display:none!important}}</style>'
            f'</head><body>{inner}</body></html>')
    fn = f"struct_atom{atom}_ICA.html"
    open(fn, "w", encoding="utf-8").write(html)
    print(f"  wrote {fn}  (atom {atom} -> {acc})")