# Figures

Figure1.pptx is the manually assembled main figure.
Figure_S1.pdf is the supplementary sequence comparison.
Its caption and reproducibility files are provided alongside it.

## Main-figure ribbon rendering

Run from the repository root:

pymol -cq scripts/figure1_final.pml
pymol -cq scripts/f1a_render.pml

The first script superposes structures using PyMOL cealign,
sets the shared orientation, and saves the session.
The second renders individual ribbons on opaque white,
using the adjusted camera and clipping planes.

Strands are coloured by element position; helices are grey
with transparency. Connecting loops and terminal extensions
are hidden for display. Coordinate inputs are not truncated
by these visualization commands.

Final cropping, labels, topology schematic and panel assembly
were performed manually in Figure1.pptx.
The rendering scripts regenerate component images, not the
complete manually assembled figure.

Dbh is the protein represented by PDB 1K1S.
