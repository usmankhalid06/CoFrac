import time, requests

# === The 297 accessions from your .mat (hardcoded; nothing is read from disk) ===
ACC = [
    "A0A0C5B5G6", "A0A1B0GTW7", "A0JNW5", "A0JP26", "A0PK11", "A1A4S6", "A1A519", "A1L190",
    "A1L3X0", "A1X283", "A2A2Y4", "A2RU14", "A2RUB6", "A2RUC4", "A3KN83", "A4D1B5",
    "A4GXA9", "A5D8V7", "A5PLL7", "A6BM72", "A6H8Y1", "A6NCS4", "A6NER3", "A6NFY7",
    "A6NGG8", "A6NI61", "A6NKB5", "A6NNB3", "A6QL63", "A7E2V4", "A7MCY6", "A7MD48",
    "A8MQ03", "A8MW99", "A9UHW6", "B1AK53", "B1AL88", "B2RUY7", "B3KU38", "B6A8C7",
    "B7U540", "C9JLW8", "C9JRZ8", "D3W0D1", "E0CX11", "O00115", "O00116", "O00159",
    "O00161", "O00165", "O00168", "O00214", "O00237", "O00254", "O00268", "O00291",
    "O00300", "O00322", "O00329", "O00330", "O00337", "O00400", "O00409", "O00422",
    "O00444", "O00453", "O00462", "O00478", "O00487", "O00505", "O00506", "O00507",
    "O00560", "O00591", "O00622", "O00624", "O00635", "O00712", "O00746", "O00748",
    "O00762", "O14497", "O14503", "O14508", "O14519", "O14524", "O14548", "O14561",
    "O14613", "O14618", "O14628", "O14638", "O14640", "O14654", "O14657", "O14744",
    "O14756", "O14757", "O14763", "O14788", "O14804", "O14813", "O14843", "O14907",
    "O14908", "O14910", "O14925", "O14926", "O14933", "O14939", "O14949", "O14958",
    "O14960", "O14967", "O14975", "O14979", "O14994", "O15067", "O15078", "O15111",
    "O15116", "O15119", "O15120", "O15127", "O15156", "O15162", "O15169", "O15230",
    "O15232", "O15245", "O15247", "O15259", "O15260", "O15263", "O15265", "O15266",
    "O15269", "O15297", "O15318", "O15370", "O15389", "O15399", "O15400", "O15446",
    "O15484", "O15488", "O15511", "O15520", "O15533", "O15540", "O15541", "O15551",
    "O15554", "O42043", "O43143", "O43155", "O43166", "O43172", "O43181", "O43196",
    "O43237", "O43264", "O43294", "O43310", "O43323", "O43395", "O43439", "O43464",
    "O43474", "O43521", "O43525", "O43526", "O43529", "O43543", "O43548", "O43556",
    "O43566", "O43581", "O43593", "O43598", "O43603", "O43609", "O43715", "O43719",
    "O43761", "O43808", "O43818", "O43823", "O43825", "O43826", "O43827", "O43854",
    "O43865", "O43866", "O43868", "O43897", "O43903", "O43913", "O43914", "O60235",
    "O60245", "O60264", "O60266", "O60281", "O60284", "O60307", "O60315", "O60381",
    "O60393", "O60449", "O60469", "O60481", "O60502", "O60519", "O60548", "O60568",
    "O60635", "O60636", "O60667", "O60678", "O60682", "O60687", "O60706", "O60729",
    "O60741", "O60759", "O60784", "O60814", "O60831", "O60879", "O60882", "O60884",
    "O60909", "O60921", "O60930", "O60938", "O60939", "O60941", "O60942", "O75019",
    "O75022", "O75072", "O75077", "O75132", "O75151", "O75159", "O75175", "O75177",
    "O75190", "O75191", "O75309", "O75319", "O75334", "O75339", "O75344", "O75365",
    "O75367", "O75369", "O75376", "O75379", "O75390", "O75396", "O75438", "O75443",
    "O75449", "O75469", "O75475", "O75486", "O75503", "O75554", "O75558", "O75564",
    "O75594", "O75610", "O75618", "O75628", "O75629", "O75638", "O75648", "O75674",
    "O75679", "O75752", "O75762", "O75781", "O75807", "O75815", "O75821", "O75828",
    "O75838", "O75864", "O75879", "O75882", "O75884", "O75897", "O75928", "O75935",
    "O75952",
]
OUT = "protein_metadata_12.tsv"
CH  = 100   # accessions per UniProt request
# ================================================================================
assert len(ACC) == 297, f"expected 297, got {len(ACC)}"
print(f"Fetching metadata for {len(ACC)} proteins from UniProt...")

fields = "accession,id,protein_name,protein_families,cc_function,sequence"
base   = "https://rest.uniprot.org/uniprotkb/stream"
header, body, got = None, [], set()
for i in range(0, len(ACC), CH):
    chunk = ACC[i:i+CH]
    q = " OR ".join(f"accession:{a}" for a in chunk)
    r = requests.get(base, params={"format": "tsv", "fields": fields, "query": q}, timeout=120)
    r.raise_for_status()
    lines = [ln for ln in r.text.splitlines() if ln.strip()]
    if header is None and lines:
        header = lines[0]
    for ln in lines[1:]:
        body.append(ln)
        got.add(ln.split("\t")[0])
    print(f"  {min(i+CH, len(ACC))}/{len(ACC)} requested, {len(body)} rows so far")
    time.sleep(0.3)

with open(OUT, "w", encoding="utf-8") as f:
    f.write(header + "\n" + "\n".join(body) + "\n")

missing = [a for a in ACC if a not in got]
print(f"\nwrote {OUT}: {len(body)} proteins")
if missing:
    print(f"{len(missing)} accessions returned no row (secondary/obsolete IDs): {missing}")
    print("-> paste these into https://www.uniprot.org/id-mapping to recover them.")