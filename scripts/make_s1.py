import csv
import re
from pathlib import Path

root = Path(__file__).resolve().parents[1]
folder = root / "results" / "tables"
report = (folder / "dssp_topology_details.txt").read_text()

names = {"PolIV_LF": "Pol IV LF", "Dbh_LF": "Dbh LF"}
elements, pairs = [], []
protein = None

for line in report.splitlines():
    match = re.match(r"^(\S+)\s+\(\d+-\d+, \d+ residues\)", line)
    if match:
        protein = names.get(match[1], match[1])
        continue

    match = re.match(r"^\s+(E\d|H\d)\s+(\d+)-(\d+)\s+\((\d+) res\)", line)
    if match:
        elements.append([protein, *match.groups()])
        continue

    match = re.match(
        r"^\s+E(\d+) - E(\d+)\s+(antiparallel|parallel)\s+"
        r"\((\d+) contacts, closest ([\d.]+) A, cos=([+\-\d.]+)\)",
        line,
    )
    if match:
        a, b, sense, count, distance, cosine = match.groups()
        pairs.append([protein, f"E{a}-E{b}", sense, count, distance, cosine])

if len(elements) != 30 or len(pairs) != 15:
    raise SystemExit(
        f"Unexpected totals: {len(elements)} elements, {len(pairs)} pairings"
    )

outputs = [
    ("S1A_secondary_structure_elements.csv",
     ["Protein", "Element", "Start residue", "End residue", "Residue count"],
     elements),
    ("S1B_strand_pairings.csv",
     ["Protein", "Strand pair", "Sense", "CA contacts <5.5 A",
      "Minimum CA distance (A)", "Axis cosine"],
     pairs),
]

for filename, header, rows in outputs:
    path = folder / filename
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"Wrote {path.name}: {len(rows)} rows")
