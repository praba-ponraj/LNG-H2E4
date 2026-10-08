#!/usr/bin/env python3
"""Combine frozen TM-align pairs about Pol IV LF; draw verified SSE boundaries.
Run: python scripts/make_figure_s1.py --root .
Requires Python 3, matplotlib. No alignment or DSSP recalculation.
"""
import argparse,csv,json,math,re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
AA=dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),'ARNDCQEGHILKMFPSTWYV'))
ORDER=['LNG','XM2H','AM2M','Pol IV LF','Dbh LF']
FILES=dict(zip(ORDER,['9bk5A.pdb','xm2h.pdb','am2m.pdb','4r8u_LF.pdb','1k1s_LF.pdb']))
PAIRS={'LNG':'tmalign_9bk5A_vs_4r8u_LF.txt','XM2H':'tmalign_xm2h_vs_4r8u_LF.txt','AM2M':'tmalign_am2m_vs_4r8u_LF.txt','Dbh LF':'tmalign_4r8u_LF_vs_1k1s_LF.txt'}
COLORS={'E1':'#0066CC','E2':'#8CB300','E3':'#E8C547','E4':'#A51E22','H1':'#999999','H2':'#666666'}
def records(path):
    result=[];seen=set()
    for line in path.read_text().splitlines():
        if line.startswith('ENDMDL'): break
        if not line.startswith('ATOM') or line[12:16].strip()!='CA' or line[21]!='A' or line[16] not in ' A':continue
        key=(int(line[22:26]),line[26].strip())
        if key in seen:continue
        seen.add(key);result.append((AA.get(line[17:20],'X'),key))
    return result

def run(root):
    structures=root/'structures';tables=root/'results/tables';raw=root/'results/tmalign';out=root/'figures';out.mkdir(parents=True,exist_ok=True)
    rec={name:records(structures/f) for name,f in FILES.items()}
    seq={name:''.join(a for a,k in rr) for name,rr in rec.items()}
    assert [len(seq[n]) for n in ORDER]==[78,89,90,109,102], {n:len(s) for n,s in seq.items()}
    sse={n:{} for n in ORDER}
    with (tables/'S1A_secondary_structure_elements.csv').open() as fh:
        for row in csv.DictReader(fh):
            n=row['Protein'];lo=int(row['Start residue']);hi=int(row['End residue'])
            for i,(aa,(r,ins)) in enumerate(rec[n]):
                if lo<=r<=hi:sse[n][i]=row['Element']
    anchor='Pol IV LF';L=len(seq[anchor]);maps={};insertions={}
    for name,f in PAIRS.items():
        lines=(raw/f).read_text().splitlines();j=next(i for i,x in enumerate(lines) if 'denotes residue pairs' in x)
        a,b=lines[j+1].strip(),lines[j+3].strip()
        if name=='Dbh LF': ref,mov=a,b
        else: mov,ref=a,b
        assert ref.replace('-','')==seq[anchor], f+' reference sequence mismatch'
        assert mov.replace('-','')==seq[name], f+' moving sequence mismatch'
        mapping={};ins={};ri=mi=0
        for r,m in zip(ref,mov):
            if r=='-':
                if m!='-':ins.setdefault(ri,[]).append(mi)
            else: mapping[ri]=mi if m!='-' else None;ri+=1
            if m!='-':mi+=1
        assert ri==L and mi==len(seq[name])
        maps[name]=mapping;insertions[name]=ins
    widths={i:max([len(ins.get(i,[])) for ins in insertions.values()]+[0]) for i in range(L+1)}
    rows={n:[] for n in ORDER};columns=[]
    for slot in range(L+1):
        for k in range(widths[slot]):
            columns.append(('insertion',slot,k))
            for n in ORDER:
                ids=insertions.get(n,{}).get(slot,[])
                rows[n].append(ids[k] if k<len(ids) else None)
        if slot<L:
            columns.append(('reference',slot,0))
            for n in ORDER:rows[n].append(slot if n==anchor else maps[n].get(slot))
    chars={n:''.join('-' if i is None else seq[n][i] for i in rows[n]) for n in ORDER}
    for n in ORDER:assert chars[n].replace('-','')==seq[n]
    (out/'Figure_S1_alignment.fasta').write_text(''.join('>'+n.replace(' ','_')+'\n'+chars[n]+'\n' for n in ORDER))
    with (out/'Figure_S1_residue_mapping.csv').open('w',newline='') as fh:
        w=csv.writer(fh);w.writerow(['Alignment column','Protein','Residue number','Insertion code','Amino acid','Element'])
        for n in ORDER:
            for c,i in enumerate(rows[n],1):
                if i is not None:w.writerow([c,n,*rec[n][i][1],rec[n][i][0],sse[n].get(i,'loop')])
    total=len(columns);chunk=50;blocks=math.ceil(total/chunk)
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(blocks,1,figsize=(9.5,blocks*2.15+.65),squeeze=False)
    for block,ax in enumerate(axes[:,0]):
        lo=block*chunk;hi=min(total,lo+chunk)
        ax.set_xlim(-13,chunk+4);ax.set_ylim(-.9,5.25);ax.axis('off')
        ax.text(-12.5,4.85,f'Alignment columns {lo+1}–{hi}',fontsize=9,color='#555555')
        for rn,n in enumerate(ORDER):
            y=4-rn
            ax.text(-1.5,y,n,ha='right',va='center',fontsize=9,fontweight='bold' if n=='Pol IV LF' else 'normal')
            for c in range(lo,hi):
                i=rows[n][c];x=c-lo
                ax.text(x,y,chars[n][c],fontfamily='DejaVu Sans Mono',ha='center',va='center',fontsize=8.3,color='#AAAAAA' if i is None else '#222222')
                if i is not None:
                    elem=sse[n].get(i)
                    if elem:
                        ax.add_patch(Rectangle((x-.46,y+.22),.92,.10,color=COLORS[elem],linewidth=0))
            # One label per element per row/block; alignment gaps may interrupt bars.
            for elem in COLORS:
                positions=[c-lo for c in range(lo,hi) if rows[n][c] is not None and sse[n].get(rows[n][c])==elem]
                if len(positions)>=3:
                    ax.text((min(positions)+max(positions))/2,y+.40,elem,ha='center',fontsize=6.5,color='#444444')
            ids=[i for i in rows[n][lo:hi] if i is not None]
            if ids:ax.text(chunk+1,y,f'{rec[n][ids[0]][1][0]}–{rec[n][ids[-1]][1][0]}',fontsize=7,va='center',color='#555555')
        ax.plot([-12.5,chunk+3],[4.62,4.62],color='#CCCCCC',lw=.5)
    fig.subplots_adjust(left=.035,right=.985,top=.97,bottom=.07,hspace=.16)
    legend='  '.join(f'{e}' for e in ['E1','E2','E3','E4','H1','H2'])
    for j,e in enumerate(['E1','E2','E3','E4','H1','H2']):
        x=.30+j*.07
        fig.text(x,.025,'━',color=COLORS[e],fontsize=14,va='center');fig.text(x+.022,.025,e,fontsize=8,va='center')
    for ext in ['pdf','svg','png']:fig.savefig(out/f'Figure_S1.{ext}',dpi=300,facecolor='white',transparent=False)
    plt.close(fig)
    (out/'Figure_S1_validation.json').write_text(json.dumps({'sequence_lengths':{n:len(s) for n,s in seq.items()},'alignment_columns':total,'validated_against_PDB':True,'reference':anchor,'pair_files':PAIRS,'insertion_policy':'Left-justified within each reference gap slot; insertion residues across rows are not structurally aligned.'},indent=2))
    print('Validated five PDB sequences and four frozen alignments;',total,'display columns.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);args=p.parse_args();run(args.root)
