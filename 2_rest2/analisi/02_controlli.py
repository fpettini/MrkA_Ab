"""Controlli preliminari delle traiettorie REST2 di MrkA (r0 e r15, finestra da T0 alla fine della traiettoria).

Uso (dalla cartella del progetto mrkA):
    python analisi_rest2/02_controlli.py . analisi_rest2/cache analisi_rest2/risultati [--t0 NS]
Legge le cache di 01_estrai.py (un sistema senza cache viene saltato). Per ogni sistema, su r0 e r15:
RMSD del CORE dalla propria conf.gro, frazione beta (DSSP, 1 frame/ns, residui 44-199), Rg dei Ca dei
sette loop riscaldati, distanza SG43-SG87. Scrive OUT/controlli.txt e OUT/controlli.json.
"""
import sys, os, json, argparse, warnings, numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP
from provenienza import prov, testo, provenienza_cache

ap = argparse.ArgumentParser()
ap.add_argument("root"); ap.add_argument("cache"); ap.add_argument("out")
ap.add_argument("--t0", type=float, default=20.0, help="scarto iniziale, ns (default 20)")
arg = ap.parse_args()
ROOT, CACHE, OUT, T0 = arg.root, arg.cache, arg.out, arg.t0
SISTEMI = ["xray", "nmr"]
REPLICHE = ["r0", "r15"]
CORE = [44, 45, 46, 47, 78, 79, 80, 81, 82, 83, 84, 85, 88, 94, 102, 103, 104, 105, 107, 108, 122, 134,
        135, 136, 137, 138, 148, 162, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 195, 196, 197, 198, 199]
HOT = {"59-77": range(59, 78), "87-101": range(87, 102),
       "112-118": [k for k in range(112, 119) if k not in (116, 117)],      # solo i 5 Ca senza 116 e 117
       "124-133": range(124, 134), "138-147": range(138, 148), "159-169": range(159, 170),
       "180-192": range(180, 193)}
HOT = {k: list(v) for k, v in HOT.items()}
BETA = (44, 199)                     # residui per la frazione beta
PASSO_DSSP_NS = 1.0
BLOCCO = 20.0                        # ns, blocchi per il SE
SOGLIA_RMSD = 2.0                    # A, flag sulla media in r15

def atomo(names, resids, nome, k):
    w = np.where((names == nome) & (resids == k))[0]
    if len(w) != 1:
        raise SystemExit(f"atomo {nome} del residuo {k} non trovato")
    return w[0]

def fit(P, Q):
    """Sovrappone ogni frame di P (F,N,3) su Q (N,3) (Kabsch con traslazione)."""
    Pm = P.mean(1, keepdims=True); Qm = Q.mean(0)
    H = np.einsum('fni,nj->fij', P - Pm, Q - Qm)
    U, S, Vt = np.linalg.svd(H)
    D = np.tile(np.eye(3), (len(P), 1, 1)); D[:, 2, 2] = np.sign(np.linalg.det(U @ Vt))
    return np.einsum('fni,fij->fnj', P - Pm, U @ D @ Vt) + Qm

def se_blocchi(y, dt, L=BLOCCO):
    """SE dai blocchi di L ns contati dalla fine (solo blocchi interi, l'avanzo si toglie all'inizio)."""
    bl = int(round(L / dt)); n = len(y) // bl
    if n < 2:
        return np.nan
    return y[len(y) - n * bl:].reshape(n, bl).mean(1).std(ddof=1) / np.sqrt(n)

def controlla(s, r):
    d = np.load(os.path.join(CACHE, f"{s}_{r}.npz"), allow_pickle=True)
    X, t, ref = d["X"], d["t"], d["ref"].astype(np.float64)
    names, resids = d["names"].astype(str), d["resids"]
    dt = np.median(np.diff(t)) / 1000
    w = t >= T0 * 1000 - 1e-6
    o = {"frame": int(w.sum())}
    ic = np.array([atomo(names, resids, "CA", k) for k in CORE])
    P = X[w][:, ic].astype(np.float64)
    rm = np.sqrt(((fit(P, ref[ic]) - ref[ic]) ** 2).sum(2).mean(1))
    o["rmsd_core_media"], o["rmsd_core_max"] = float(rm.mean()), float(rm.max())
    # frazione beta: universo con soli N, CA, C, O (stessa numerazione di conf.gro)
    u0 = mda.Universe(os.path.join(ROOT, "rest2", s, "conf.gro"))
    bb = u0.select_atoms("name N CA C O")
    ib = np.where(np.isin(names, ["N", "CA", "C", "O"]))[0]
    if not (np.array_equal(bb.names, names[ib]) and np.array_equal(bb.resids, resids[ib])):
        raise SystemExit("atomi di backbone della cache diversi da conf.gro")
    uu = mda.Merge(bb)
    passo = int(round(PASSO_DSSP_NS / dt))
    uu.load_new(np.ascontiguousarray(X[w][::passo][:, ib]), format=mda.coordinates.memory.MemoryReader)
    e = DSSP(uu.select_atoms("protein")).run().results.dssp == "E"
    o["beta"] = float(e[:, BETA[0] - 42:BETA[1] - 41].mean()); o["frame_dssp"] = int(len(e))
    # Rg dei Ca di ogni loop caldo, con SE a blocchi
    o["rg"] = {}
    for l, rs in HOT.items():
        C = X[w][:, [atomo(names, resids, "CA", k) for k in rs]].astype(np.float64)
        rg = np.sqrt(((C - C.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        o["rg"][l] = {"media": float(rg.mean()), "se": float(se_blocchi(rg, dt))}
    ss = np.linalg.norm(X[w][:, atomo(names, resids, "SG", 43)].astype(np.float64)
                        - X[w][:, atomo(names, resids, "SG", 87)], axis=1)
    o["ss_media"], o["ss_max"] = float(ss.mean()), float(ss.max())
    return o

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    ins = [os.path.join(CACHE, f"{s}_{r}.npz") for s in SISTEMI for r in REPLICHE] + \
          [os.path.join(ROOT, "rest2", s, "conf.gro") for s in SISTEMI]
    P = prov(__file__, ins, vars(arg), {"ROOT": ROOT, "CACHE": CACHE, "OUT": OUT}); P["cache"] = provenienza_cache(CACHE, SISTEMI, REPLICHE)
    R = {"provenienza": P, "T0_ns": T0, "sistemi": {}, "flag": []}
    righe = testo(P) + [f"Controlli preliminari, finestra da {T0:g} ns alla fine della traiettoria", ""]
    for s in SISTEMI:
        mancanti = [r for r in REPLICHE if not os.path.exists(os.path.join(CACHE, f"{s}_{r}.npz"))]
        if mancanti:
            righe.append(f"sistema {s}: cache assente per {', '.join(mancanti)}, salto"); continue
        c = {r: controlla(s, r) for r in REPLICHE}
        R["sistemi"][s] = c
        righe.append(f"=== {s} ===")
        righe.append(f"RMSD dei {len(CORE)} Ca del CORE dalla propria conf.gro (fit sul CORE), A")
        for r in REPLICHE:
            righe.append(f"   {r:4s} media {c[r]['rmsd_core_media']:.2f}  max {c[r]['rmsd_core_max']:.2f}  ({c[r]['frame']} frame)")
        if c["r15"]["rmsd_core_media"] > SOGLIA_RMSD:
            f = f"{s}: RMSD CORE medio in r15 {c['r15']['rmsd_core_media']:.2f} A > {SOGLIA_RMSD:.1f} A"
            R["flag"].append(f); righe.append("   FLAG " + f)
        righe.append(f"Frazione beta (DSSP 'E', residui {BETA[0]}-{BETA[1]}, 1 frame/ns)")
        for r in REPLICHE:
            righe.append(f"   {r:4s} {c[r]['beta']:.3f}  ({c[r]['frame_dssp']} frame)")
        righe.append("Rg dei Ca dei loop caldi, A: media (SE a blocchi da 20 ns)")
        for l in HOT:
            a, b = c["r0"]["rg"][l], c["r15"]["rg"][l]
            soglia = a["media"] - 2 * np.hypot(a["se"], b["se"])
            flag = bool(b["media"] < soglia)
            riga = f"   {l:8s} r0 {a['media']:.2f} ({a['se']:.2f})  r15 {b['media']:.2f} ({b['se']:.2f})"
            if flag:
                f = f"{s}: Rg loop {l} in r15 ({b['media']:.2f}) < r0 - 2 SE ({soglia:.2f})"
                R["flag"].append(f); riga += "  FLAG"
            R["sistemi"][s].setdefault("rg_flag", {})[l] = flag
            righe.append(riga)
        righe.append(f"Distanza SG43-SG87 in r0, A: media {c['r0']['ss_media']:.2f}  max {c['r0']['ss_max']:.2f}")
        righe.append("")
    righe.append("Flag: " + ("; ".join(R["flag"]) if R["flag"] else "nessuno"))
    open(os.path.join(OUT, "controlli.txt"), "w").write("\n".join(righe) + "\n")
    json.dump(R, open(os.path.join(OUT, "controlli.json"), "w"), indent=1)
    print("\n".join(righe))
