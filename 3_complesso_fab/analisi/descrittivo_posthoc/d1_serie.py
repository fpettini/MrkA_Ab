"""DESCRITTIVO (post hoc, 11 nov 2026): serie per fotogramma per caratterizzare lo scopo 2.
Per ogni replica: secondaria (stessa funzione di 13), scostamento per residuo dell'epitopo dopo fit sull'epitopo
e dopo fit sul core (backbone N, CA, C, O), diedri phi/psi di 86-102, distanze scelte.
Non entra in nessun criterio: serve solo a descrivere che cosa cambia quando il libero esce dall'involucro."""
import os, sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # comune_fab.py sta nella cartella sopra
import MDAnalysis as mda
import comune_fab as C
from multiprocessing import Pool
TRJ = "trj_tutti"; PDB = "mrka-fab_completo.pdb"; OUT = "descrittivo"
RES_DIED = list(range(86, 103))

def lavora(arg):
    s, r, xtc = arg
    gro = os.path.join(TRJ, f"proteina_{s}.gro")
    top = C.Topologia(gro, s, PDB)
    piano, gruppi = C.piano_riassemblaggio(top)
    i_epi, r_epi = C.atomi_epitopo(top, PDB)
    xyz, _ = C.leggi_pdb(PDB)
    i_core = np.array([top.atomo("A", k, "CA") for k in C.CORE43]); r_core = np.array([xyz["A"][(k, "CA")] for k in C.CORE43])
    idx = {(k, n): top.atomo("A", k, n) for k in range(85, 104) for n in ("N", "CA", "C")}
    # distanze: Asp40 OD1/OD2 - Ala90 N ; Cys43 SG - Cys87 SG ; Lys94 NZ - Asp96 CG... solo intra MrkA
    def at(k, n):
        return top.atomo("A", k, n)
    coppie = {"D40OD1_A90N": (at(40, "OD1"), at(90, "N")), "D40OD2_A90N": (at(40, "OD2"), at(90, "N")),
              "C43SG_C87SG": (at(43, "SG"), at(87, "SG")), "CA91_CA168": (at(91, "CA"), at(168, "CA")),
              "CA92_CA42": (at(92, "CA"), at(42, "CA")), "CA91_CA96": (at(91, "CA"), at(96, "CA")),
              "CA90_CA94": (at(90, "CA"), at(94, "CA")), "CA92_CA88": (at(92, "CA"), at(88, "CA"))}
    u = mda.Universe(gro, xtc)
    nf = u.trajectory.n_frames
    T = np.empty(nf); sec = np.empty(nf); prim = np.empty(nf)
    dev_self = np.empty((nf, len(C.EPITOPO))); dev_core = np.empty((nf, len(C.EPITOPO)))
    phi = np.empty((nf, len(RES_DIED))); psi = np.empty((nf, len(RES_DIED)))
    dist = {k: np.empty(nf) for k in coppie}
    nb = len(C.BACKBONE)
    for a in range(0, nf, 200):
        b = min(nf, a + 200)
        X = np.empty((b - a, top.n, 3)); M = np.empty((b - a, 3, 3))
        for j, ts in enumerate(u.trajectory[a:b]):
            X[j] = ts.positions; M[j] = ts.triclinic_dimensions; T[a + j] = ts.time
        Y = C.riassembla(X, M, piano, gruppi)
        P = Y[:, i_epi]
        R, cP, cQ = C.kabsch(P, r_epi)
        F = C.applica(P, R, cP, cQ)
        sec[a:b] = C.rmsd(F, r_epi)
        d2 = ((F - r_epi) ** 2).sum(-1)                       # (f, 44)
        dev_self[a:b] = np.sqrt(d2.reshape(b - a, len(C.EPITOPO), nb).mean(-1))
        Rc, cPc, cQc = C.kabsch(Y[:, i_core], r_core)
        Fc = C.applica(P, Rc, cPc, cQc)
        prim[a:b] = C.rmsd(Fc, r_epi)
        d2c = ((Fc - r_epi) ** 2).sum(-1)
        dev_core[a:b] = np.sqrt(d2c.reshape(b - a, len(C.EPITOPO), nb).mean(-1))
        for q, k in enumerate(RES_DIED):
            phi[a:b, q] = C.diedro(Y[:, idx[(k - 1, "C")]], Y[:, idx[(k, "N")]], Y[:, idx[(k, "CA")]], Y[:, idx[(k, "C")]])
            psi[a:b, q] = C.diedro(Y[:, idx[(k, "N")]], Y[:, idx[(k, "CA")]], Y[:, idx[(k, "C")]], Y[:, idx[(k + 1, "N")]])
        for nome, (i, j) in coppie.items():
            dist[nome][a:b] = np.linalg.norm(Y[:, i] - Y[:, j], axis=-1)
    t = np.round((T - T[0]) / 1000.0, 4)
    np.savez_compressed(os.path.join(OUT, f"d1_{s}_{r}.npz"), t_ns=t, sec=sec, prim=prim, dev_self=dev_self, dev_core=dev_core,
                        phi=phi, psi=psi, res_died=np.array(RES_DIED), epitopo=np.array(C.EPITOPO), **dist)
    return f"{s} {r}: {nf} fotogrammi, sec media {sec.mean():.2f}, prim media {prim.mean():.2f}"

if __name__ == "__main__":
    lavori = [("complesso", f"run{i}", os.path.join(TRJ, f"complesso_run{i}_0-500ns.xtc")) for i in range(1, 6)] + \
             [("libero", f"run{i}", os.path.join(TRJ, f"libero_run{i}_0-1000ns.xtc")) for i in range(1, 11)]
    with Pool(4) as p:
        for x in p.imap_unordered(lavora, lavori):
            print(x, flush=True)
