#!/usr/bin/env python3
"""Count frozen CATH50 hits by superfamily at query TM-score >= 0.5."""
from pathlib import Path
from collections import Counter
import csv

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "results/foldseek"
THRESHOLD = 0.5

with (FOLDER / "cath_target_superfamilies.tsv").open() as handle:
    mapping = {
        row["target"]: row["superfamily"]
        for row in csv.DictReader(handle, delimiter="\t")
    }

summary = []
for stem in ["9bk5A", "xm2h", "am2m", "4r8u_LF", "1k1s_LF"]:
    counts = Counter()
    with (FOLDER / f"{stem}_cath.tsv").open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if float(row[3]) >= THRESHOLD:
                counts[mapping[row[1]]] += 1

    total = sum(counts.values())
    print(f"\n{stem}: {total} hits at qTM >= {THRESHOLD}")
    for sf, count in counts.most_common():
        print(f"  {sf}: {count}")
        summary.append([stem, sf, count, total, THRESHOLD])

output = ROOT / "results/tables/cath_superfamily_counts.csv"
with output.open("w") as handle:
    writer = csv.writer(handle)
    writer.writerow(["Query", "Superfamily", "Hits", "Total qualifying hits", "qTM threshold"])
    writer.writerows(summary)

print(f"\nSaved {output}")
