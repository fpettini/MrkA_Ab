"""Densita' del complesso in s2: globale (massa da topologia / volume) e dell'acqua di bulk (lontano dalla proteina).
Uso: python densita.py topol.top min_c3.gro equil_s1.gro equil_s2.gro equil_s2.xtc"""
import sys, re, numpy as np, mdtraj as md
from scipy.spatial import cKDTree
NA = 6.02214076e23
top, g0, g1, g2, xtc = sys.argv[1:6]; TDA = float(sys.argv[6])
# massa della proteina da topol.top ([ atoms ] del moleculetype della proteina, colonna massa)
sez=None; mprot=0.0; nprot=0
for r in open(top):
    r0=r.split(';')[0].strip()
    if not r0 or r0.startswith('#'): continue
    m=re.match(r'\[\s*(\S+)\s*\]',r0)
    if m: sez=m.group(1); continue
    if sez=='atoms': c=r0.split(); mprot+=float(c[7]); nprot+=1
mol={}
for r in open(top).read().split('[ molecules ]')[1].splitlines():
    c=r.split(';')[0].split()
    if len(c)==2: mol[c[0]]=mol.get(c[0],0)+int(c[1])
MW, MK, MCL = 16.0+2*1.008, 39.1, 35.45
nw, nk, ncl = mol['SOL'], mol['K'], mol['CL']
mtot = mprot + nw*MW + nk*MK + ncl*MCL
print(f"massa: proteina {mprot:.1f} Da ({nprot} atomi), acque {nw} x {MW:.3f}, K {nk}, Cl {ncl}; totale {mtot/1000:.1f} kDa")
def vol_gro(p):
    b=[float(x) for x in open(p).read().splitlines()[-1].split()]
    return b[0]*b[1]*b[2]   # nm3: matrice triangolare inferiore, volume = prodotto della diagonale
def rho(v): return mtot/NA/(v*1e-21)   # g/cm3
for nome,p in (('min_c3 / s1 (NVT)',g0),('fine s1',g1),('fine s4',g2)):
    v=vol_gro(p); print(f"{nome:18s} volume {v:9.2f} nm3  densita' globale {rho(v):.4f} g/cm3")
# traiettoria s2: volume per frame e acqua di bulk
t = md.load(xtc, top=g2)
V = t.unitcell_volumes  # nm3
print(f"s2 xtc: {t.n_frames} frame, {t.time[0]:.0f}-{t.time[-1]:.0f} ps; volume medio ultimi 2.5 ns {V[t.time>=TDA].mean():.2f} nm3 -> densita' globale {rho(V[t.time>=TDA].mean()):.4f} g/cm3")
q=len(V)//5
print("densita' globale per quinti di s2: " + " ".join(f"{rho(V[i*q:(i+1)*q].mean()):.4f}" for i in range(5)))
topo=t.topology
prot=np.array([a.index for a in topo.atoms if a.residue.name not in ('SOL','K','CL','HOH') and a.element.symbol!='H' and a.name!='MW'])
ow=np.array([a.index for a in topo.atoms if a.residue.is_water and a.element.symbol=='O'], dtype=int)
ioni=np.array([a.index for a in topo.atoms if a.residue.name in ('K','CL')])
mion=np.array([MK if topo.atom(i).residue.name=='K' else MCL for i in ioni])
print(f"atomi pesanti proteina {len(prot)}, OW {len(ow)}, ioni {len(ioni)}")
rng=np.random.default_rng(1)
def avvolgi(x, box):
    a,b,c = box
    s=np.floor(x[:,2]/c[2]); x=x-s[:,None]*c
    s=np.floor(x[:,1]/b[1]); x=x-s[:,None]*b
    s=np.floor(x[:,0]/a[0]); x=x-s[:,None]*a
    return x
shift=np.array([[i,j,k] for i in (-1,0,1) for j in (-1,0,1) for k in (-1,0,1)])
ris={1.0:[],1.2:[],1.5:[]}
fr=[f for f in range(t.n_frames) if t.time[f]>=TDA][::2]
for f in fr:
    box=t.unitcell_vectors[f]; x=t.xyz[f]
    P=avvolgi(x[prot].astype(np.float64),box); W=avvolgi(x[ow].astype(np.float64),box); I=avvolgi(x[ioni].astype(np.float64),box)
    imm=(P[None,:,:]+(shift@box)[:,None,:]).reshape(-1,3)
    albero=cKDTree(imm)
    dw,_=albero.query(W); di,_=albero.query(I)
    L=np.array([box[0][0],box[1][1],box[2][2]]); npt=1_000_000
    pts=rng.random((npt,3))*L; dp,_=albero.query(pts)
    vbox=L.prod()
    for R in ris:
        vb=vbox*(dp>R).mean()
        mw=(dw>R).sum()*MW; mi=mion[di>R].sum()
        ris[R].append(((mw)/NA/(vb*1e-21), (mw+mi)/NA/(vb*1e-21), (dw>R).sum(), vb))
for R,v in ris.items():
    v=np.array(v)
    print(f"bulk (OW a piu' di {R} nm dalla proteina, {len(fr)} frame ultima meta'): acqua {v[:,0].mean():.4f} +- {v[:,0].std(ddof=1)/np.sqrt(len(v)):.4f} g/cm3, "
          f"soluzione (acqua+ioni) {v[:,1].mean():.4f}; ~{v[:,2].mean():.0f} acque in {v[:,3].mean():.0f} nm3")
