# Structural convergence of LNG and SCUBA designs on a natural H2E4 topology

Analysis supporting a Protein Science Research Note comparing
the RFdiffusion-derived α-cobratoxin binder LNG with the
SCUBA architecture-directed designs XM2H and AM2M, and the
natural little-finger domains of Pol IV and Dbh.

All five experimental structures share the principal element
order E1–H1–E2–E3–H2–E4 and an antiparallel sheet with spatial
strand order E1–E4–E2–E3. XM2H and AM2M derive from the same
SCUBA sketch and are treated as one architecture-directed
design observation.

Shared topology does not imply identical geometry or establish
that all designs belong to the same natural superfamily.
The structural comparisons and database searches characterize
these distinctions separately.

## Repository contents

- structures/: analyzed coordinate inputs, original PDB entries
  and a description of analyzed extents.
- scripts/: topology analysis, table generation, CATH hit counting,
  supplementary figure generation and PyMOL rendering.
- results/tables/: frozen structural comparisons, topology summary,
  DSSP details, Tables S1A, S1B and S2, and CATH hit counts.
- results/tmalign/: raw TM-align outputs.
- results/foldseek/: complete saved hit lists, target-to-superfamily
  mapping, search commands and database version records.
- figures/: main figure source and supplementary figure, with
  captions and figure-specific documentation.

## Supplementary material

- Figure S1: structure-guided sequence comparison of all five proteins.
- Table S1A: principal secondary-structure element boundaries.
- Table S1B: strand-pair geometry and contact counts.
- Table S2: 16 pairwise structural comparisons, including both
  TM-score normalizations, aligned lengths and sequence identities.

Table legends are in results/tables/Supplementary_table_legends.md.

## Reproduce summaries

Run from the repository root:

python scripts/make_s1.py
python scripts/make_s2.py
python scripts/count_cath.py
python scripts/make_figure_s1.py --root .

These commands use the saved analysis outputs. They do not
rerun TM-align or Foldseek searches.

Figure S1 requires matplotlib. The other summary scripts use
the Python standard library. Topology analysis requires NumPy,
Biotite and mkdssp; the recorded environment used NumPy 2.5.1,
Biotite 1.7.1 and mkdssp 4.6.1.

Main-figure ribbon rendering requires PyMOL. See figures/README.md.
The final main figure was assembled manually.

## Interpretation and provenance

Sequence identity is measured over TM-align-aligned residues.
Query- and target-normalized TM-scores describe the same
alignment using different normalization lengths.

AAA in topology.tsv means that all three adjacent strand pairs
in the spatial sheet order are antiparallel.

CATH summaries count returned hits at query-normalized TM-score
>= 0.5. Counts describe the saved search results and their
superfamily distribution, rather than formal classification
of the designed proteins.

Figure S1 combines four pairwise alignments around Pol IV LF.
Insertion placement is for display and does not establish
equivalence between residues in other rows.

See the folder READMEs for coordinate provenance, search
parameters, database versions and visualization procedures.

## Citation and archive

Citation metadata, licensing and the release-specific Zenodo
DOI will be finalized before the archived release.
