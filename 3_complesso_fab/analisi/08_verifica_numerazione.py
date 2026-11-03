# Verifica della numerazione dei residui nei .gro (criteri di analisi, punto 9; 2 nov 2026).
# Uso, dalla cartella complesso_fab:  python3 08_verifica_numerazione.py libero/run1/run1.gro complesso/run1/run1.gro ...
# Confronta ogni catena proteica del gro con mrka-fab_completo.pdb (sequenza, primo e ultimo numero) e stampa
# i nomi dei residui di controllo. Nei gro le Cys a ponte sono CYX e i terminali hanno prefisso N/C (NGLY, CCYX).
import sys, hashlib
T={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
   'CYX':'C','HIE':'H','HID':'H','HIP':'H','LYN':'K','ASH':'D','GLH':'E','CYM':'C'}
def one(n):
    if n in T: return T[n]
    if len(n)==4 and n[0] in 'NC' and n[1:] in T: return T[n[1:]]
    return None
def pdb_chains(p):
    ch={}
    for l in open(p):
        if l.startswith('ATOM') and l[17:20] in T:
            ch.setdefault(l[21],{})[int(l[22:26])]=T[l[17:20]]
    return ch
def gro_chains(p):
    L=open(p).read().splitlines()
    n=int(L[1]); res=[]
    for l in L[2:2+n]:
        k=(int(l[0:5]),l[5:10].strip())
        if not res or res[-1][0]!=k[0] or res[-1][1]!=k[1]: res.append(k)
    prot=[r for r in res if one(r[1])]
    chains=[];cur=[]
    for r in prot:
        if cur and r[0]!=cur[-1][0]+1: chains.append(cur);cur=[]
        cur.append(r)
    if cur: chains.append(cur)
    altri={}
    for r in res:
        if not one(r[1]): altri[r[1]]=altri.get(r[1],0)+1
    return chains,altri
ref=pdb_chains('mrka-fab_completo.pdb')
for p in sys.argv[1:]:
    h=hashlib.sha256(open(p,'rb').read()).hexdigest()[:12]
    chains,altri=gro_chains(p)
    print(f'== {p} (sha256 {h}): {len(chains)} catene proteiche; non proteina: {altri}')
    for i,c in enumerate(chains):
        seq=''.join(one(r[1]) for r in c)
        tag=[k for k,v in ref.items() if len(v)==len(c) and ''.join(v[x] for x in sorted(v))==seq and min(v)==c[0][0] and max(v)==c[-1][0]]
        print(f'  catena {i+1}: {c[0][0]}{c[0][1]}-{c[-1][0]}{c[-1][1]}, {len(c)} residui, uguale a pdb catena {tag or "NESSUNA"}')
        d={r[0]:r[1] for r in c}
        if i==0:
            print('    A87',d.get(87),' A94',d.get(94),' A43',d.get(43),' A42',d.get(42),' A40',d.get(40),' A116',d.get(116),' A203',d.get(203),' A217',d.get(217))
            print('    epitopo',[f"{x}{d.get(x)}" for x in [42,86,88,90,91,92,93,94,96,167,168]])
        else:
            print('    ',{x:d.get(x) for x in [1,22,23,87,88,92,95,96,101,102,134,194,200,204,214,220,223] if x in d})
