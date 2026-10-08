# Coordinate inputs

| Analyzed file | Protein | Source | Analyzed extent |
|---|---|---|---|
| 9bk5A.pdb | LNG | 9BK5:A | Residues 8–85; 78 residues |
| xm2h.pdb | XM2H | 7DKK | Chain A; 89 residues |
| am2m.pdb | AM2M | 7DKO | Chain A; 90 residues |
| 4r8u_LF.pdb | Pol IV little-finger domain | 4R8U:A | Residues 233–341; 109 residues |
| 1k1s_LF.pdb | Dbh little-finger domain | 1K1S:A | Residues 240–341; 102 residues |

Lengths count coordinate residues in the analyzed chain or
residue selection, including connecting loops.

The analyzed files preserve the inputs used for the frozen
TM-align and Foldseek comparisons. The specified LNG extent
is also the selection used for topology analysis.

Original downloaded entries are retained as 9bk5.pdb,
7dkk.pdb, 7dko.pdb, 4r8u.pdb and 1k1s.pdb.

1k1s_A_coordinates.pdb contains the amino-acid coordinates
of chain A, model 1, with hydrogen atoms removed. It was
written to permit direct DSSP processing without the original
entry's problematic metadata.

DSSP assignments were mapped using chain, residue number and
insertion code. The corrected full-chain Dbh assignment omits
residue A35, which lies outside the little-finger domain.
