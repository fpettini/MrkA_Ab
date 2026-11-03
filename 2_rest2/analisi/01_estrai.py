"""Estrazione delle traiettorie REST2 di MrkA con la proteina ricomposta (PBC riparate).

Uso (dalla cartella del progetto mrkA):
    python analisi_rest2/01_estrai.py . analisi_rest2/cache [NOME_RUN]
NOME_RUN e' la cartella del run dentro rest2/<sistema>/ (default rest2_16x480; per REST3 per esempio
rest3_16x480_k105, con una cache sua).
Per ogni sistema con traiettoria (xray, nmr) e per r0 e r15 legge xtc + conf.gro, ricompone la catena
(ogni atomo all'immagine minima rispetto all'atomo a cui e' legato) e salva in
CACHE/<sistema>_<replica>.npz solo N, CA, C, O, H di backbone e SG di Cys43/Cys87.
Campi: X (frame, atomi, 3; A), t (ps), names, resids, resnames, ref (da conf.gro).
Il sistema nmr e' ancora in produzione: se manca l'xtc viene saltato.
"""
import sys, os, time, json, warnings, numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from provenienza import prov

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
CACHE = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "analisi_rest2", "cache")
RUN = sys.argv[3] if len(sys.argv) > 3 else "rest2_16x480"
SISTEMI = ["xray", "nmr"]
REPLICHE = ["r0", "r15"]
BB = ["N", "CA", "C", "O", "H"]
SS_RES = (43, 87)
PEP_MAX = 2.0                      # A, legame C(i)-N(i+1) massimo ammesso

def minimo(v, M):
    """Immagine minima dei vettori v (F,3) con le scatole triclinee M (F,3,3, righe = vettori di cella)."""
    f = np.einsum('fi,fij->fj', v, np.linalg.inv(M))
    return v - np.einsum('fi,fij->fj', np.round(f), M)

def ricomponi(X, M, genitore):
    """Percorre la catena: ogni atomo = genitore + vettore a immagine minima (in ordine di indice)."""
    Y = X.astype(np.float64).copy()
    for a, p in enumerate(genitore):
        if p >= 0:
            Y[:, a] = Y[:, p] + minimo(X[:, a].astype(np.float64) - X[:, p].astype(np.float64), M)
    return Y

def estrai(s, r):
    xtc = os.path.join(ROOT, "rest2", s, RUN, r, "rest2.xtc")
    gro = os.path.join(ROOT, "rest2", s, "conf.gro")
    if not os.path.exists(xtc):
        print(f"   {s} {r}: traiettoria assente ({xtc}), salto")
        return
    u = mda.Universe(gro, xtc)                       # mai il .tpr: rinumera i residui
    g0 = mda.Universe(gro)
    if u.trajectory.n_atoms != g0.atoms.n_atoms:
        raise SystemExit(f"{s} {r}: numero di atomi xtc diverso da conf.gro")
    sel = [i for i, (n, k) in enumerate(zip(u.atoms.names, u.atoms.resids)) if n in BB or (n == "SG" and k in SS_RES)]
    ag = u.atoms[sel]
    names, resids = ag.names, ag.resids
    # genitore di ogni atomo scelto (indice nel gruppo), -1 = ancora
    pos = {(n, k): i for i, (n, k) in enumerate(zip(names, resids))}
    gen = []
    for i, (n, k) in enumerate(zip(names, resids)):
        p = {"N": ("C", k - 1), "CA": ("N", k), "C": ("CA", k), "O": ("C", k),
             "H": ("N", k), "SG": ("CA", k)}[n]
        g = pos.get(p, -1)
        if g >= i:
            raise SystemExit("ordine degli atomi inatteso in conf.gro")
        gen.append(g)
    nf = u.trajectory.n_frames
    X = np.empty((nf, len(ag), 3), np.float32); M = np.empty((nf, 3, 3)); t = np.empty(nf)
    for i, ts in enumerate(u.trajectory):
        X[i] = ag.positions; M[i] = ts.triclinic_dimensions; t[i] = ts.time
    X = ricomponi(X, M, gen).astype(np.float32)
    # controllo: legami peptidici C(i)-N(i+1) in tutti i frame
    ks = [k for k in np.unique(resids) if (("C", k) in pos) and (("N", k + 1) in pos)]
    ic = [pos[("C", k)] for k in ks]; inn = [pos[("N", k + 1)] for k in ks]
    pep = np.linalg.norm(X[:, ic].astype(np.float64) - X[:, inn], axis=2).max()
    if not pep < PEP_MAX:
        raise SystemExit(f"{s} {r}: legame C-N massimo {pep:.2f} A dopo la ricomposizione (>= {PEP_MAX})")
    P = prov(__file__, [xtc, gro], {"root": ROOT, "cache": CACHE, "run": RUN, "sistema": s, "replica": r},
             {"ROOT": ROOT, "CACHE": CACHE})
    np.savez(os.path.join(CACHE, f"{s}_{r}.npz"), X=X, t=t, names=np.array(names), resids=np.array(resids),
             resnames=np.array(ag.resnames), ref=g0.atoms[sel].positions.astype(np.float32),
             provenienza=np.array(json.dumps(P)))
    print(f"   {s:5s} {r:3s}: {nf} frame, {t[0] / 1000:.2f}-{t[-1] / 1000:.2f} ns, passo {np.diff(t).min():.1f}-"
          f"{np.diff(t).max():.1f} ps, {len(ag)} atomi, legame C-N max {pep:.2f} A")

if __name__ == "__main__":
    os.makedirs(CACHE, exist_ok=True)
    t0 = time.time()
    for s in SISTEMI:
        for r in REPLICHE:
            estrai(s, r)
    print(f"tempo totale {time.time() - t0:.0f} s")
