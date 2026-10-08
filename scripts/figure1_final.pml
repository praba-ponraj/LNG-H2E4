reinitialize
set ray_opaque_background, 0
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_flat_sheets, 1
set ray_shadows, 0
bg_color white

load structures/9bk5A.pdb,   LNG
load structures/xm2h.pdb,    XM2H
load structures/am2m.pdb,    AM2M
load structures/4r8u_LF.pdb, PolIV
load structures/1k1s_LF.pdb, Dbh

hide everything
show cartoon
color grey80

color marine,       LNG   and resi 9-13
color grey40,       LNG   and resi 16-28
color splitpea,     LNG   and resi 34-42
color yelloworange, LNG   and resi 45-53
color grey40,       LNG   and resi 57-74
color firebrick,    LNG   and resi 78-83

color marine,       XM2H  and resi 2-8
color grey40,       XM2H  and resi 13-30
color splitpea,     XM2H  and resi 36-42
color yelloworange, XM2H  and resi 45-51
color grey40,       XM2H  and resi 57-76
color firebrick,    XM2H  and resi 82-88

color marine,       AM2M  and resi 2-7
color grey40,       AM2M  and resi 13-28
color splitpea,     AM2M  and resi 34-42
color yelloworange, AM2M  and resi 45-53
color grey40,       AM2M  and resi 57-77
color firebrick,    AM2M  and resi 84-89

color marine,       PolIV and resi 243-254
color grey40,       PolIV and resi 257-278
color splitpea,     PolIV and resi 286-293
color yelloworange, PolIV and resi 298-304
color grey40,       PolIV and resi 310-324
color firebrick,    PolIV and resi 330-338

color marine,       Dbh  and resi 246-256
color grey40,       Dbh  and resi 259-275
color splitpea,     Dbh  and resi 283-290
color yelloworange, Dbh  and resi 294-301
color grey40,       Dbh  and resi 308-323
color firebrick,    Dbh  and resi 331-339

cealign PolIV, LNG
cealign PolIV, XM2H
cealign PolIV, AM2M
cealign PolIV, Dbh

hide cartoon, not (ss H or ss S)
hide cartoon, PolIV and not resi 243-338
hide cartoon, Dbh  and not resi 246-339

set cartoon_transparency, 0.65, LNG   and resi 16-28+57-74
set cartoon_transparency, 0.65, XM2H  and resi 13-30+57-76
set cartoon_transparency, 0.65, AM2M  and resi 13-28+57-77
set cartoon_transparency, 0.65, PolIV and resi 257-278+310-324
set cartoon_transparency, 0.65, Dbh  and resi 259-275+308-323

set_view (\
    -0.216248170,    0.426443726,    0.878281951,\
    -0.617454767,   -0.756561458,    0.215322644,\
     0.756297827,   -0.495736063,    0.426914960,\
     0.000009924,    0.000227161, -126.878303528,\
     0.964183569,   17.261358261,   46.333900452,\
   104.216819763,  149.545715332,  -20.000000000 )

python
from pymol import cmd
objs = ["LNG","XM2H","AM2M","PolIV","Dbh"]
front = cmd.get_view()
for obj in objs:
    cmd.disable("all"); cmd.enable(obj); cmd.set_view(front)
    cmd.ray(1100, 1300); cmd.png("figures/f1_front_%s.png" % obj, dpi=600)
cmd.disable("all"); cmd.enable("all"); cmd.set_view(front)
cmd.turn("y", 180)
back = cmd.get_view()
for obj in objs:
    cmd.disable("all"); cmd.enable(obj); cmd.set_view(back)
    cmd.ray(1100, 1300); cmd.png("figures/f1_back_%s.png" % obj, dpi=600)
cmd.enable("all"); cmd.set_view(front)
cmd.save("figures/figure1_final.pse")
python end
