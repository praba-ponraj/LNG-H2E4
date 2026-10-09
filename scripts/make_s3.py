from pathlib import Path
from collections import Counter
import csv, hashlib, json

ROOT = Path(__file__).resolve().parents[1]
F = ROOT / "results/foldseek"
OUT = ROOT / "results/tables"
NAMES = {
    "9bk5A": "LNG", "xm2h": "XM2H", "am2m": "AM2M",
    "4r8u_LF": "Pol IV LF", "1k1s_LF": "Dbh LF"
}

def read(name):
    with (F / name).open() as h:
        return list(csv.reader(h, delimiter="\t"))

def write(name, header, rows):
    with (OUT / name).open("w", newline="") as h:
        w = csv.writer(h)
        w.writerow(header)
        w.writerows(rows)

def accession(target):
    return target.split("-")[1]

mapping = {
    r[0]: r[1] for r in read("cath_target_superfamilies.tsv")[1:]
}

cath = []
for stem, name in NAMES.items():
    counts = Counter(
        mapping[r[1]] for r in read(stem + "_cath.tsv")
        if float(r[3]) >= 0.5
    )
    for sf, count in counts.most_common():
        cath.append([
            name, sf, count, sum(counts.values()), stem + "_cath.tsv"
        ])

control = {
    accession(r[1]) for r in read("4r8u_LF_afdbsp.tsv")
    if float(r[3]) >= 0.6
}

af, lng = [], []
for stem in ["4r8u_LF", "xm2h", "am2m", "9bk5A"]:
    source = stem + "_afdbsp.tsv"
    rows = read(source)
    qualifying = [r for r in rows if float(r[3]) >= 0.5]
    ids = {accession(r[1]) for r in qualifying}
    af.append([
        NAMES[stem], len(rows), len(qualifying),
        sum(float(r[3]) >= 0.6 for r in rows),
        f"{max(float(r[3]) for r in rows):.4f}",
        len(ids), len(ids & control), source
    ])
    if stem == "9bk5A":
        for r in qualifying:
            lng.append([
                accession(r[1]), r[1], f"{float(r[3]):.4f}",
                f"{float(r[6]):.3f}",
                "Yes" if accession(r[1]) in control else "No",
                source, rows.index(r) + 1
            ])

selected = []
for database, source, targets in [
    ("PDB100", "9bk5A_pdb.tsv",
     ["4r8u-assembly2_B", "1k1s-assembly1_A"]),
    ("CATH50", "9bk5A_cath.tsv",
     ["af_Q9VNX1_325_451_3.30.1490.100"])
]:
    rows = read(source)
    ranked = sorted(rows, key=lambda r: float(r[3]), reverse=True)
    for target in targets:
        r = next(r for r in rows if r[1] == target)
        selected.append([
            database, target, len(rows), ranked.index(r) + 1,
            f"{float(r[3]):.4f}", f"{float(r[4]):.4f}",
            f"{float(r[6]):.3f}", source, rows.index(r) + 1
        ])

assert [r[3] for r in selected] == [197, 590, 5]
assert len(lng) == 7
assert sum(r[4] == "Yes" for r in lng) == 3

write("Table_S3A_CATH_superfamily_summary.csv",
      ["Query", "CATH superfamily", "Hits qTM >=0.5",
       "Total qualifying hits", "Source TSV"], cath)

write("Table_S3B_AFDB_search_summary.csv",
      ["Query", "Returned rows", "Rows qTM >=0.5", "Rows qTM >=0.6",
       "Maximum qTM", "Unique accessions qTM >=0.5",
       "Shared with control qTM >=0.6", "Source TSV"], af)

write("Table_S3C_selected_LNG_search_hits.csv",
      ["Database", "Target identifier", "Returned rows",
       "Rank by descending qTM", "qTM", "tTM", "Probability",
       "Source TSV", "Original TSV row"], selected)

write("Table_S3D_LNG_AFDB_accessions.csv",
      ["UniProt accession", "Target identifier", "qTM", "Probability",
       "Shared with control qTM >=0.6", "Source TSV",
       "Original TSV row"], lng)

sources = (
    list(F.glob("*_cath.tsv")) + list(F.glob("*_afdbsp.tsv")) +
    [F / "9bk5A_pdb.tsv", F / "cath_target_superfamilies.tsv"]
)
checks = {
    p.name: hashlib.sha256(p.read_bytes()).hexdigest()
    for p in sources
}
(OUT / "source_checksums.json").write_text(
    json.dumps(checks, indent=2) + "\n"
)

print("Verified and wrote all four Table S3 CSV files.")
