#!/usr/bin/env python3
"""Build Supplementary Table S2 from the frozen TM-align table.

Appends the natural-to-natural baseline row if absent, then writes
S2 as CSV and Markdown with readable protein names.
"""
import csv
import os
import re

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TSV = ROOT / "results/tables/final_structure_comparison.tsv"
BASELINE_TXT = ROOT / "results/tmalign/tmalign_4r8u_LF_vs_1k1s_LF.txt"
OUT_CSV = ROOT / "results/tables/S2_structural_comparisons.csv"
OUT_MD = ROOT / "results/tables/S2_structural_comparisons.md"

NAMES = {
    "9bk5A.pdb":    ("LNG", "RFdiffusion", "domain"),
    "xm2h.pdb":     ("XM2H", "SCUBA", "domain"),
    "am2m.pdb":     ("AM2M", "SCUBA", "domain"),
    "4r8u_LF.pdb":  ("Pol IV LF", "natural", "domain"),
    "1k1s_LF.pdb":  ("Dbh LF", "natural", "domain"),
    "4r8u.pdb":     ("Pol IV", "natural", "full chain"),
    "1k1s.pdb":     ("Dbh", "natural", "full chain"),
}


def parse_baseline(path):
    """Pull aligned length, RMSD, Seq_ID and both TM-scores from TM-align output."""
    text = open(path).read()
    m = re.search(
        r"Aligned length=\s*(\d+),\s*RMSD=\s*([\d.]+),.*?=\s*([\d.]+)", text)
    tms = re.findall(r"TM-score=\s*([\d.]+)", text)
    if not m or len(tms) < 2:
        raise SystemExit(f"Could not parse {path}")
    return ["4r8u_LF.pdb", "1k1s_LF.pdb", m.group(1), m.group(2),
            m.group(3), tms[0], tms[1]]


def main():
    for p in (TSV, BASELINE_TXT):
        if not os.path.exists(p):
            raise SystemExit(f"Missing: {p}  (check the repository inputs)")

    with open(TSV) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r]
    header, data = rows[0], rows[1:]

    have_baseline = any(
        {r[0], r[1]} == {"4r8u_LF.pdb", "1k1s_LF.pdb"} for r in data)
    if not have_baseline:
        new = parse_baseline(BASELINE_TXT)
        data.append(new)
        with open(TSV, "a") as fh:
            fh.write("\t".join(new) + "\n")
        print(f"Appended baseline row -> {TSV}")
    else:
        print("Baseline row already present; not re-appending.")
    print(f"{len(data)} comparison rows.")

    unknown = {r[i] for r in data for i in (0, 1)} - set(NAMES)
    if unknown:
        print("WARNING: unnamed structures, add to NAMES:", sorted(unknown))

    out_header = ["Query", "Query origin", "Query extent",
                  "Target", "Target origin", "Target extent",
                  "Aligned length", "RMSD (A)", "Seq ID (%)",
                  "TM (query-normalized)", "TM (target-normalized)"]

    def fmt(r):
        q = NAMES.get(r[0], (r[0], "", ""))
        t = NAMES.get(r[1], (r[1], "", ""))
        return [q[0], q[1], q[2], t[0], t[1], t[2], r[2],
                f"{float(r[3]):.2f}", f"{100 * float(r[4]):.1f}",
                f"{float(r[5]):.3f}", f"{float(r[6]):.3f}"]

    # natural-natural first, then designs vs natural domains, then the rest
    def rank(r):
        q, t = NAMES.get(r[0], ("", "", ""))[1], NAMES.get(r[1], ("", "", ""))[1]
        if q == "natural" and t == "natural":
            return 0
        if t == "natural":
            return 1
        return 2

    body = [fmt(r) for r in sorted(data, key=lambda r: (rank(r), r[0], r[1]))]

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(out_header)
        w.writerows(body)

    widths = [max(len(h), *(len(row[i]) for row in body))
              for i, h in enumerate(out_header)]
    with open(OUT_MD, "w") as fh:
        fh.write("| " + " | ".join(h.ljust(widths[i])
                 for i, h in enumerate(out_header)) + " |\n")
        fh.write("|" + "|".join("-" * (w + 2) for w in widths) + "|\n")
        for row in body:
            fh.write("| " + " | ".join(c.ljust(widths[i])
                     for i, c in enumerate(row)) + " |\n")

    print(f"Wrote {OUT_CSV}\nWrote {OUT_MD}")


if __name__ == "__main__":
    main()
