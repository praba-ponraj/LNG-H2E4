# Frozen Foldseek searches

Foldseek version currently installed: 10.941cd33.
The original search commands were recovered from shell history.
The executable version at the time of each search was not separately logged.

## Databases

- CATH50: v4.3.0
- AlphaFold/Swiss-Prot: v6, 2025_03
- PDB100: PDB_DATE 250101

Complete database version files are retained in database_versions/.
The FOLDSEEK_COMMIT recorded in pdb.version describes database provenance,
not necessarily the executable used for these searches.

## TSV columns

The files have no header. Columns are:
1. query
2. target
3. alntmscore
4. qtmscore
5. ttmscore
6. lddt
7. prob
8. evalue

CATH superfamily counts use qtmscore >= 0.5.

## Original commands

CATH50, for 9bk5A, xm2h, am2m, 4r8u_LF and 1k1s_LF:

foldseek easy-search QUERY.pdb foldseek_db/cath50 QUERY_cath.tsv tmp --format-output "query,target,alntmscore,qtmscore,ttmscore,lddt,prob,evalue" -e 10

AlphaFold/Swiss-Prot, for 9bk5A, xm2h, am2m and 4r8u_LF:

foldseek easy-search QUERY.pdb foldseek_db/afdb_sp QUERY_afdbsp.tsv tmp --format-output "query,target,alntmscore,qtmscore,ttmscore,lddt,prob,evalue" -e 10

PDB100, for 9bk5A and xm2h:

foldseek easy-search QUERY.pdb foldseek_db/pdb QUERY_pdb.tsv tmp --format-output "query,target,alntmscore,qtmscore,ttmscore,lddt,prob,evalue" --exhaustive-search 1 -e 10

QUERY denotes each listed input stem. Searches used the original
PDB inputs; their analyzed extents are documented separately.
