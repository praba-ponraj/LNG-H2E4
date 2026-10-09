# Figure S1 reproduction package

This package adds Figure S1 to the LNG-H2E4 repository. Extract its contents into the existing repository root, preserving the folder structure.

## Reproduce

From the repository root:

```bash
python scripts/make_figure_s1.py --root .
```

Requires Python 3 and matplotlib. The script uses the supplied coordinate files, four frozen TM-align outputs and corrected Table S1A boundaries. No structural alignment is rerun.

## Outputs

- `figures/Figure_S1.pdf`: vector figure for submission.
- `figures/Figure_S1.png`: preview on an opaque white background.
- `figures/Figure_S1.svg`: editable vector figure.
- `figures/Figure_S1_caption.txt`: proposed supplementary caption.
- `figures/Figure_S1_alignment.fasta`: displayed reference-anchored alignment.
- `figures/Figure_S1_residue_mapping.csv`: display columns mapped to PDB residue numbers and elements.
- `figures/Figure_S1_validation.json`: input sequence and alignment validation.

## Interpretation

The four pairwise alignments share Pol IV LF as their reference. Insertion slots are padded and left-aligned for display only; this is not an evolutionary multiple-sequence alignment. Pairwise identities must be taken from the frozen TM-align results, not recomputed from this display. Original residue numbering is retained. The secondary-structure bars use the corrected DSSP boundaries, including Dbh E1 at residues 246–256.

Repository documentation, licensing and database provenance are provided in the main README, LICENSE files and results/foldseek/README.md.
