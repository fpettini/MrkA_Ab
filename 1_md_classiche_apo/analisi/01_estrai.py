"""Fase 1: estrae dalle traiettorie solo gli atomi della proteina e li salva in .npz.

Uso (dalla cartella del progetto mrkA):  python analisi_convergenza/01_estrai.py . analisi_convergenza/cache
Le traiettorie attese sono traj/<sistema>_run<r>_pr.xtc; per un altro suffisso: SUFFISSO=... python ...
Richiede MDAnalysis >= 2.7.
"""
import sys, os, numpy as np
from multiprocessing import Pool
import MDAnalysis as mda

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "analisi_convergenza", "cache")
SISTEMI = {"xray": "xray", "nmr_ox": "nmr_ox", "nmr_red": "nmr_red"}
REPLICHE = (1, 2, 3)
SUFFISSO = os.environ.get("SUFFISSO", "_pr")   # traj/<sistema>_run<r><SUFFISSO>.xtc

def n_proteina(sistema):
    # conf.gro e' l'uscita di pdb2gmx: contiene solo la proteina, nello stesso
    # ordine dei primi atomi di ions.gro e quindi della traiettoria.
    with open(os.path.join(ROOT, sistema, "conf.gro")) as f:
        f.readline()
        return int(f.readline())

def estrai(arg):
    sistema, r = arg
    n = n_proteina(sistema)
    xtc = os.path.join(ROOT, "traj", f"{sistema}_run{r}{SUFFISSO}.xtc")
    with open(xtc, "rb") as f:                     # l'intestazione XTC dice quanti atomi ci sono
        natoms = int.from_bytes(f.read(8)[4:8], "big")
    # traiettoria di sola proteina (trjconv sul gruppo Protein) o del sistema intero
    top = "conf.gro" if natoms == n else "ions.gro"
    u = mda.Universe(os.path.join(ROOT, sistema, top), xtc)
    prot = u.atoms[:n]
    rif = mda.Universe(os.path.join(ROOT, sistema, "conf.gro")).atoms
    assert list(prot.names) == list(rif.names), "ordine atomi diverso da conf.gro"
    X = np.empty((len(u.trajectory), n, 3), dtype=np.float32)
    t = np.empty(len(u.trajectory))
    for i, ts in enumerate(u.trajectory):
        X[i] = prot.positions
        t[i] = ts.time
    np.savez(os.path.join(OUT, f"{sistema}_run{r}.npz"), X=X, t=t,
             names=np.array(prot.names).astype(str), resids=prot.resids, resnames=np.array(prot.resnames).astype(str),
             ref=rif.positions.astype(np.float32))
    return f"{sistema}_run{r}: {len(t)} frame, {t[0]:.0f}-{t[-1]:.0f} ps, {n} atomi"

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    lavori = [(s, r) for s in SISTEMI for r in REPLICHE]
    with Pool(min(4, os.cpu_count())) as p:
        for msg in p.imap_unordered(estrai, lavori):
            print(msg, flush=True)
