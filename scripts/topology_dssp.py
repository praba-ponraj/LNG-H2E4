#!/usr/bin/env python3
"""
topology_dssp.py - H2E4 topology comparison using DSSP hydrogen-bond
secondary structure assignment.

DSSP is run on the complete chain; any residue range is applied afterwards,
so H-bond context is never truncated.

Usage:
    python topology_dssp.py LABEL=file.pdb:CHAIN[:START-END] ...
"""

import sys
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb
from biotite.application.dssp import DsspApp

PAIR_CUTOFF = 5.5
MIN_PAIRS = 2

# DSSP code -> simplified class
HELIX = {"H", "G", "I"}
STRAND = {"E", "B"}


def load(spec):
    label, rest = spec.split("=", 1)
    parts = rest.split(":")
    path, chain = parts[0], parts[1]
    rng = None
    if len(parts) > 2:
        lo, hi = (int(x) for x in parts[2].split("-"))
        rng = (lo, hi)

    arr = pdb.PDBFile.read(path).get_structure(model=1)
    arr = arr[struc.filter_amino_acids(arr) & (arr.chain_id == chain)]
    arr = arr[arr.element != "H"]
    if arr.array_length() == 0:
        raise SystemExit(f"{label}: nothing selected from {path} chain {chain}")
    return label, arr, rng


def dssp_segments(arr, rng):
    """Run full-chain DSSP and map assignments by explicit residue identifiers."""
    import subprocess
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp:
        inp = Path(tmp) / "coordinates.pdb"
        out = Path(tmp) / "assignments.dssp"
        pf = pdb.PDBFile()
        pf.set_structure(arr)
        pf.write(str(inp))
        subprocess.run(
            ["mkdssp", "--output-format=dssp", str(inp), str(out)],
            check=True, capture_output=True, text=True
        )
        assignments = {}
        started = False
        for line in out.read_text().splitlines():
            if line.lstrip().startswith("#"):
                started = True
                continue
            if not started or len(line) < 17 or line[13] == "!":
                continue
            try:
                key = (line[11], int(line[5:10]), line[10].strip())
            except ValueError:
                continue
            assignments[key] = line[16].strip() or "C"

    starts = struc.get_residue_starts(arr)
    keys = [
        (str(arr.chain_id[i]), int(arr.res_id[i]), str(arr.ins_code[i]).strip())
        for i in starts
    ]
    missing = [key for key in keys if key not in assignments]
    if missing:
        print(f"DSSP omitted residues (assigned C for segmentation): {missing}")

    ca = arr[arr.atom_name == "CA"]
    coords_by_key = {
        (str(c), int(r), str(ins).strip()): xyz
        for c, r, ins, xyz in zip(ca.chain_id, ca.res_id, ca.ins_code, ca.coord)
    }
    keys = [
        key for key in keys
        if rng is None or rng[0] <= key[1] <= rng[1]
    ]
    codes = np.array([assignments.get(key, "C") for key in keys])
    res_ids = np.array([key[1] for key in keys])
    simple = [
        "a" if code in HELIX else "b" if code in STRAND else "c"
        for code in codes
    ]

    segs = []
    start = 0
    for i in range(1, len(keys) + 1):
        boundary = i == len(keys)
        if not boundary:
            boundary = (
                simple[i] != simple[start]
                or keys[i][0] != keys[i-1][0]
                or keys[i][1] != keys[i-1][1] + 1
            )
        if boundary:
            kind = simple[start]
            selected = keys[start:i]
            minimum = 4 if kind == "a" else 3
            if kind in ("a", "b") and len(selected) >= minimum:
                xyz = np.array([
                    coords_by_key[key] for key in selected
                    if key in coords_by_key
                ])
                if len(xyz) >= 2:
                    segs.append((kind, selected[0][1], selected[-1][1], xyz))
            start = i
    return segs, codes, res_ids


def topology_string(segs):
    out, e, h = [], 0, 0
    for kind, _, _, _ in segs:
        if kind == "b":
            e += 1
            out.append(f"E{e}")
        else:
            h += 1
            out.append(f"H{h}")
    return "-".join(out)


def strand_axis(coords):
    v = coords[-1] - coords[0]
    norm = np.linalg.norm(v)
    return v / norm if norm > 1e-6 else v


def sheet_graph(strands):
    n = len(strands)
    pairs = {}
    for i in range(n):
        for j in range(i + 1, n):
            ci, cj = strands[i][3], strands[j][3]
            d = np.linalg.norm(ci[:, None, :] - cj[None, :, :], axis=-1)
            npairs = int((d < PAIR_CUTOFF).sum())
            if npairs >= MIN_PAIRS:
                cos = float(np.dot(strand_axis(ci), strand_axis(cj)))
                pairs[(i, j)] = ("parallel" if cos > 0 else "antiparallel",
                                 npairs, cos, float(d.min()))
    return pairs


def spatial_order(n, pairs):
    adj = {i: [] for i in range(n)}
    for (i, j) in pairs:
        adj[i].append(j)
        adj[j].append(i)
    ends = [i for i in range(n) if len(adj[i]) == 1]
    if not ends:
        return None
    start = min(ends)
    order, seen, cur = [start], {start}, start
    while True:
        nxt = [k for k in adj[cur] if k not in seen]
        if not nxt:
            break
        cur = nxt[0]
        order.append(cur)
        seen.add(cur)
    return order if len(order) == n else None


def report(label, arr, rng):
    segs, codes, res_ids = dssp_segments(arr, rng)
    strands = [s for s in segs if s[0] == "b"]
    helices = [s for s in segs if s[0] == "a"]

    span = f"{res_ids.min()}-{res_ids.max()}" if len(res_ids) else "empty"
    print(f"\n{'=' * 66}\n{label}  ({span}, {len(res_ids)} residues)\n{'=' * 66}")
    print(f"topology        : {topology_string(segs)}")
    print(f"composition     : {len(helices)} helices, {len(strands)} strands "
          f"-> H{len(helices)}E{len(strands)}")

    print("\nelements:")
    e = h = 0
    for kind, lo, hi, _ in segs:
        if kind == "b":
            e += 1
            name = f"E{e}"
        else:
            h += 1
            name = f"H{h}"
        print(f"  {name:4s} {lo:>4d}-{hi:<4d}  ({hi - lo + 1} res)")

    print("\nDSSP string:")
    s = "".join(codes)
    for i in range(0, len(s), 60):
        print(f"  {int(res_ids[i]):>4d}  {s[i:i+60]}")

    pairs = sheet_graph(strands)
    if not pairs:
        print("\nno strand pairs within cutoff")
        return

    print("\nsheet pairing:")
    for (i, j), (sense, npairs, cos, dmin) in sorted(pairs.items()):
        print(f"  E{i+1} - E{j+1}  {sense:<13s} "
              f"({npairs} contacts, closest {dmin:.1f} A, cos={cos:+.2f})")

    order = spatial_order(len(strands), pairs)
    if order:
        print(f"\nstrand order in sheet : "
              f"{' - '.join(f'E{i+1}' for i in order)}")
        senses = []
        for a, b in zip(order, order[1:]):
            senses.append("A" if pairs[(min(a,b), max(a,b))][0] == "antiparallel"
                          else "P")
        print(f"sense along sheet     : {' '.join(senses)}  "
              f"({'all antiparallel' if set(senses) == {'A'} else 'mixed'})")
    else:
        print("\nstrand order: could not linearize")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    for spec in sys.argv[1:]:
        label, arr, rng = load(spec)
        report(label, arr, rng)
    print()


if __name__ == "__main__":
    main()
