"""POST HOC, 12 nov 2026 (dopo i risultati finali): descrizione per finestra (50-300 e 300-500 ns) dei contatti del cristallo
con la catena pesante e con la leggera e della punta 90-93 nelle 5 repliche del complesso. NON e' un criterio.
Il taglio a 300 ns viene dal controllo del 10 nov (a occhio, sulla replica 2): e' descrittivo, non un tempo misurato.
Legge le cache di 10_serie.py (contatti per fotogramma) e diagnostica_cerniera.json di 17 (un fotogramma ogni 10)."""
import os, sys, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # comune_fab.py sta nella cartella sopra
import comune_fab as C
W = sys.argv[1]; pdb = f"{W}/mrka-fab_completo.pdb"
_, rn = C.leggi_pdb(pdb)
cris = sorted(C.contatti_cristallo(pdb, 4.0))
H = [p for p in cris if p[1] == "H"]; L = [p for p in cris if p[1] == "L"]
print(f"contatti del cristallo a 4.0 A: {len(cris)} ({len(H)} con la catena pesante, {len(L)} con la leggera)")

def fraz(d, w, cut=4.0):
    fr = np.where(w)[0]
    sel = np.isin(d["c_f"], fr) & (d["c_d"] <= cut)
    ch, n = np.unique((d["c_a"][sel].astype(np.int64) * 2 + d["c_c"][sel]) * 1000 + d["c_k"][sel], return_counts=True)
    return {(int(x // 2000), "LH"[int(x // 1000 % 2)], int(x % 1000)): float(m) / len(fr) for x, m in zip(ch, n)}

dj = json.load(open(f"{W}/risultati/diagnostica_cerniera.json"))
print("\nmedia delle frazioni dei contatti del cristallo (pesante | leggera) e punta 90-93 (RMSD dei Ca, A) per finestra")
print("replica   finestra   pesante  leggera   tip|Fv  tip|core  rot Fv")
tab = {}
for r in ["run1", "run2", "run3", "run4", "run5"]:
    z = np.load(f"{W}/cache/serie_complesso_{r}.npz", allow_pickle=True)
    d = {k: z[k] for k in ("t_ns", "c_f", "c_a", "c_c", "c_k", "c_d")}
    t = d["t_ns"]; x = dj["repliche"][r]; tj = np.array(x["t_ns"]); s = x["serie"]
    for a, b in ((50, 300), (300, 500)):
        w = (t >= a) & ((t < b) if b < 500 else (t <= b)); wj = (tj >= a) & ((tj < b) if b < 500 else (tj <= b))
        f = fraz(d, w); tab[(r, a)] = f
        print(f"{r:8s}  {a:3d}-{b:3d}    {np.mean([f.get(p, 0.0) for p in H]):.2f}     {np.mean([f.get(p, 0.0) for p in L]):.2f}     "
              f"{np.mean(np.array(s['tip_fv'])[wj]):.2f}    {np.mean(np.array(s['tip_core'])[wj]):.2f}    {np.mean(np.array(s['rot_fv'])[wj]):.1f}")
print("\nreplica 2, per coppia (50-300 | 300-500 ns)")
for p in cris:
    print(f"  A{p[0]}{rn['A'][p[0]]}-{p[1]}{p[2]}{rn[p[1]][p[2]]}: {tab[('run2', 50)].get(p, 0.0):.2f} | {tab[('run2', 300)].get(p, 0.0):.2f}")
