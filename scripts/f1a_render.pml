load figures/figure1_final.pse

bg_color white
set ray_opaque_background, 1
set antialias, 2
set cartoon_fancy_helices, 1

set_view (\
    -0.216248170,    0.426443726,    0.878281951,\
    -0.617454767,   -0.756561458,    0.215322644,\
     0.756297827,   -0.495736063,    0.426914960,\
     0.000009924,    0.000227161, -150.000000000,\
     0.964183569,   17.261358261,   46.333900452,\
   104.216819763,  199.545715332,  -20.000000000 )

python
from pymol import cmd
v = cmd.get_view()
for o in ["LNG", "XM2H", "AM2M", "PolIV", "Dbh"]:
    cmd.disable("all")
    cmd.enable(o)
    cmd.set_view(v)
    cmd.png("figures/f1a_%s_clean.png" % o,
            width=1100, height=1500, dpi=600, ray=1)
python end
