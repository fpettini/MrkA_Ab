"""Verifica indipendente dei numeri chiave REST2 di MrkA (scritta da zero).

Uso:  python verifica.py            (dalla cartella mrkA o da qualunque posto)

Legge solo rest2/<sistema>/conf.gro e le xtc indicate; scrive verifica.txt e
verifica.json accanto a questo script.

Scelte esplicite (non dettate dalla consegna):
- finestra [a, b] ns: frame con a*1000 <= t(ps) <= b*1000, estremi inclusi;
- fit di Kabsch non pesato (tutti gli atomi del CORE con peso 1);
- in caso di parita' esatta dX == dN6 il frame va in N (regola "X se dX < dN6, N se no");
- O(i) e' ricomposto all'immagine minima di C(i); la catena N-CA-C e' percorsa
  per residui in ordine, ogni atomo all'immagine minima del precedente.
"""
import hashlib
import json
import os
import platform
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import minimize_vectors

QUI = os.path.dirname(os.path.abspath(__file__))
RADICE = os.path.abspath(os.path.join(QUI, "..", ".."))
REST2 = os.path.join(RADICE, "rest2")

GRO = {"xray": os.path.join(REST2, "xray", "conf.gro"),
       "nmr": os.path.join(REST2, "nmr", "conf.gro")}
TRAIETTORIE = [  # (etichetta, sistema, xtc)
    ("xray_r0", "xray", os.path.join(REST2, "xray", "rest2_16x480", "r0", "rest2.xtc")),
    ("nmr_r0", "nmr", os.path.join(REST2, "nmr", "rest2_16x480", "r0", "rest2.xtc")),
    ("nmr_r15", "nmr", os.path.join(REST2, "nmr", "rest2_16x480", "r15", "rest2.xtc")),
]
FINESTRE = [(20.0, 400.0), (120.0, 400.0)]

RES_PRIMO, RES_ULTIMO = 42, 202


def intervallo(a, b):
    return list(range(a, b + 1))


CORE = (intervallo(44, 47) + intervallo(78, 85) + [88, 94] + intervallo(102, 105)
        + [107, 108, 122] + intervallo(134, 138) + [148, 162]
        + intervallo(170, 179) + intervallo(195, 199))
assert len(CORE) == 43, len(CORE)

LOOP = {
    "59-77": intervallo(59, 77),
    "87-101": intervallo(87, 101),
    "112-118(-116,117)": [r for r in intervallo(112, 118) if r not in (116, 117)],
    "124-133": intervallo(124, 133),
    "138-147": intervallo(138, 147),
    "159-169": intervallo(159, 169),
    "180-192": intervallo(180, 192),
    "48-56": intervallo(48, 56),
}
NOMI_LOOP = ["N", "CA", "C", "O"]


def sha256(percorso, blocco=1 << 22):
    h = hashlib.sha256()
    with open(percorso, "rb") as f:
        while True:
            b = f.read(blocco)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- topologia
def indici(universo):
    """Indici (0-based nell'universo) per catena backbone, O, CORE e loop."""
    pro = universo.atoms
    mappa = {}
    for a in pro:
        chiave = (int(a.resid), a.name)
        if chiave in mappa:
            raise RuntimeError(f"atomo duplicato {chiave}")
        mappa[chiave] = a.index
    resid = sorted({int(r) for r in pro.resids})
    if resid != intervallo(RES_PRIMO, RES_ULTIMO):
        raise RuntimeError(f"residui inattesi: {resid[0]}..{resid[-1]} ({len(resid)})")
    catena = []  # N, CA, C per residuo in ordine
    for r in resid:
        for n in ("N", "CA", "C"):
            catena.append(mappa[(r, n)])
    ossigeni = [mappa[(r, "O")] for r in resid]
    return mappa, np.array(catena), np.array(ossigeni), resid


class Selezioni:
    """Posizioni ordinate di backbone ricomposto e sottoinsiemi CORE/loop."""

    def __init__(self, universo):
        self.mappa, self.catena, self.ossigeni, self.resid = indici(universo)
        self.nres = len(self.resid)
        # array ricomposto: per residuo k -> righe 3k (N), 3k+1 (CA), 3k+2 (C); O in coda
        riga = {}
        for k, r in enumerate(self.resid):
            riga[(r, "N")] = 3 * k
            riga[(r, "CA")] = 3 * k + 1
            riga[(r, "C")] = 3 * k + 2
            riga[(r, "O")] = 3 * self.nres + k
        self.core = np.array([riga[(r, "CA")] for r in CORE])
        self.loop = {nome: np.array([riga[(r, n)] for r in res for n in NOMI_LOOP])
                     for nome, res in LOOP.items()}
        # C(i) -> N(i+1) nella catena
        self.c_righe = np.array([3 * k + 2 for k in range(self.nres - 1)])
        self.n_righe = self.c_righe + 1

    def ricomponi(self, pos, box):
        """Backbone N/CA/C/O integro: ritorna (xyz (4*nres,3), max d C(i)-N(i+1))."""
        grezzo = pos[self.catena].astype(np.float64)
        passi = minimize_vectors((grezzo[1:] - grezzo[:-1]).astype(np.float32), box)
        catena = np.empty_like(grezzo)
        catena[0] = grezzo[0]
        catena[1:] = grezzo[0] + np.cumsum(passi.astype(np.float64), axis=0)
        c = catena[2::3]  # C di ogni residuo
        dO = minimize_vectors((pos[self.ossigeni].astype(np.float64) - c).astype(np.float32), box)
        o = c + dO.astype(np.float64)
        xyz = np.vstack([catena, o])
        dcn = np.linalg.norm(xyz[self.n_righe] - xyz[self.c_righe], axis=1)
        return xyz, float(dcn.max())


# ---------------------------------------------------------------- geometria
def kabsch(mobile, rif):
    """R, cm, cr tali che (mobile - cm) @ R.T + cr sovrappone mobile a rif."""
    cm = mobile.mean(axis=0)
    cr = rif.mean(axis=0)
    H = (mobile - cm).T @ (rif - cr)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, cm, cr


def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def d_loop(xyz, rif, sel):
    """Fit su CA del CORE di rif, poi RMSD di ogni loop senza nuovo fit."""
    R, cm, cr = kabsch(xyz[sel.core], rif[sel.core])
    mosso = (xyz - cm) @ R.T + cr
    return {n: rmsd(mosso[i], rif[i]) for n, i in sel.loop.items()}, rmsd(mosso[sel.core], rif[sel.core])


def carica_riferimento(sistema):
    u = mda.Universe(GRO[sistema])
    sel = Selezioni(u)
    xyz, dcn = sel.ricomponi(u.atoms.positions, u.dimensions)
    nomi_res = [u.select_atoms(f"resid {r}").residues[0].resname for r in sel.resid]
    return u, sel, xyz, dcn, nomi_res


# ---------------------------------------------------------------- traiettorie
def analizza(compito):
    etichetta, sistema, xtc = compito
    _, sel, rif_prop, _, _ = carica_riferimento(sistema)
    _, _, rifX, _, _ = carica_riferimento("xray")
    _, _, rifN, _, _ = carica_riferimento("nmr")
    u = mda.Universe(GRO[sistema], xtc)
    nf = len(u.trajectory)
    t = np.empty(nf)
    dX = {n: np.empty(nf) for n in LOOP}
    dN = {n: np.empty(nf) for n in LOOP}
    core_prop = np.empty(nf)
    dcn = np.empty(nf)
    for k, ts in enumerate(u.trajectory):
        t[k] = ts.time
        xyz, dcn[k] = sel.ricomponi(ts.positions, ts.dimensions)
        a, _ = d_loop(xyz, rifX, sel)
        b, _ = d_loop(xyz, rifN, sel)
        for n in LOOP:
            dX[n][k] = a[n]
            dN[n][k] = b[n]
        R, cm, cr = kabsch(xyz[sel.core], rif_prop[sel.core])
        core_prop[k] = rmsd((xyz[sel.core] - cm) @ R.T + cr, rif_prop[sel.core])
    return {"etichetta": etichetta, "sistema": sistema, "xtc": xtc, "nframe": nf,
            "t": t, "dX": dX, "dN": dN, "core_prop": core_prop, "dcn": dcn,
            "box_triclino": bool(np.any(np.abs(u.trajectory[0].dimensions[3:] - 90.0) > 1e-3))}


def popolazioni(ris, D6, a, b):
    m = (ris["t"] >= a * 1000.0) & (ris["t"] <= b * 1000.0)
    out = {"nframe": int(m.sum())}
    for n in LOOP:
        x = ris["dX"][n][m]
        y = ris["dN"][n][m]
        altro = np.minimum(x, y) > D6[n]
        cx = (~altro) & (x < y)
        cn = (~altro) & ~(x < y)
        out[n] = {"pX": float(cx.mean()), "pN": float(cn.mean()), "pAltro": float(altro.mean()),
                  "dX_medio": float(x.mean()), "dN6_medio": float(y.mean())}
    return out


def main():
    script = os.path.abspath(__file__)
    righe = []
    scrivi = righe.append

    # riferimenti
    uX, selX, X, dcnX, resX = carica_riferimento("xray")
    uN, selN, N6, dcnN, resN = carica_riferimento("nmr")
    if uX.atoms.n_atoms != uN.atoms.n_atoms:
        raise RuntimeError("conf.gro con numero di atomi diverso")
    res_diversi = [(r, a, b) for r, a, b in zip(selX.resid, resX, resN) if a != b]
    D6, core_XN6 = d_loop(X, N6, selX)
    D6_inv, _ = d_loop(N6, X, selX)

    # controllo dell'implementazione di Kabsch contro MDAnalysis (solo come test)
    from MDAnalysis.analysis.rms import rmsd as mda_rmsd
    core_mda = float(mda_rmsd(X[selX.core], N6[selX.core], center=True, superposition=True))

    with ProcessPoolExecutor(max_workers=len(TRAIETTORIE)) as ex:
        risultati = list(ex.map(analizza, TRAIETTORIE))
    ris = {r["etichetta"]: r for r in risultati}

    # --- intestazione / provenienza
    scrivi("VERIFICA INDIPENDENTE REST2 MrkA")
    scrivi("=" * 78)
    scrivi(f"script: {script}")
    scrivi(f"sha256 script: {sha256(script)}")
    scrivi(f"python {platform.python_version()} | numpy {np.__version__} | MDAnalysis {mda.__version__}"
           f" | {platform.platform()}")
    scrivi("")
    scrivi("File letti (sha256):")
    for s in ("xray", "nmr"):
        scrivi(f"  {sha256(GRO[s])}  {os.path.relpath(GRO[s], RADICE)}")
    for r in risultati:
        scrivi(f"  {sha256(r['xtc'])}  {os.path.relpath(r['xtc'], RADICE)}")
    scrivi("")
    scrivi("Traiettorie: frame, tempo iniziale/finale, passo, box, C(i)-N(i+1) massimo")
    for r in risultati:
        dt = np.diff(r["t"])
        scrivi(f"  {r['etichetta']:8s} nframe={r['nframe']}  t0={r['t'][0]/1000:.3f} ns  "
               f"tfin={r['t'][-1]/1000:.3f} ns  dt={np.median(dt):.1f} ps "
               f"(min {dt.min():.1f}, max {dt.max():.1f})  triclino={r['box_triclino']}  "
               f"maxC-N={r['dcn'].max():.3f} A  frame con C-N>=2.0 A: {int((r['dcn'] >= 2.0).sum())}")
    scrivi(f"  conf.gro xray: maxC-N={dcnX:.3f} A ; conf.gro nmr: maxC-N={dcnN:.3f} A")
    scrivi(f"  residui con nome diverso fra X e N6: {res_diversi if res_diversi else 'nessuno'}")
    scrivi("")
    scrivi("Scelte: finestra a<=t<=b (ns, inclusiva); Kabsch non pesato; dX==dN6 -> N;"
           " O all'immagine minima di C dello stesso residuo.")
    scrivi("")

    # --- riferimenti
    scrivi("RIFERIMENTI X (xray/conf.gro) vs N6 (nmr/conf.gro)")
    scrivi("-" * 78)
    scrivi(f"RMSD CORE (43 CA) X vs N6 dopo fit: {core_XN6:.4f} A  (MDAnalysis per controllo: {core_mda:.4f} A)")
    scrivi(f"{'loop':20s} {'atomi':>5s} {'D6=d(X,N6)':>11s} {'d(N6,X)':>9s}")
    for n in LOOP:
        scrivi(f"{n:20s} {len(selX.loop[n]):5d} {D6[n]:11.4f} {D6_inv[n]:9.4f}")
    scrivi("")

    # --- popolazioni
    pop = {}
    for (a, b) in FINESTRE:
        chiave = f"[{a:g}, {b:g}]"
        pop[chiave] = {e: popolazioni(ris[e], D6, a, b) for e in ("xray_r0", "nmr_r0")}
        px, pn = pop[chiave]["xray_r0"], pop[chiave]["nmr_r0"]
        scrivi(f"POPOLAZIONI r0, finestra {chiave} ns  (frame: xray {px['nframe']}, nmr {pn['nframe']})")
        scrivi("-" * 78)
        scrivi(f"{'loop':20s} | {'xray pX':>7s} {'pN':>6s} {'pAltro':>6s} | {'nmr pX':>7s} {'pN':>6s} {'pAltro':>6s}"
               f" | {'Δ pX':>7s} {'Δ pN':>7s}")
        for n in LOOP:
            dpx = px[n]["pX"] - pn[n]["pX"]
            dpn = px[n]["pN"] - pn[n]["pN"]
            pop[chiave].setdefault("delta", {})[n] = {"dpX": dpx, "dpN": dpn}
            scrivi(f"{n:20s} | {px[n]['pX']:7.4f} {px[n]['pN']:6.4f} {px[n]['pAltro']:6.4f} | "
                   f"{pn[n]['pX']:7.4f} {pn[n]['pN']:6.4f} {pn[n]['pAltro']:6.4f} | {dpx:+7.4f} {dpn:+7.4f}")
        scrivi("")

    # --- RMSD del CORE dalla propria conf.gro
    core = {}
    scrivi("RMSD CA del CORE dalla propria conf.gro (fit sul CORE), finestra [20, 400] ns")
    scrivi("-" * 78)
    for e in ("xray_r0", "nmr_r0", "nmr_r15"):
        r = ris[e]
        m = (r["t"] >= 20000.0) & (r["t"] <= 400000.0)
        v = r["core_prop"][m]
        core[e] = {"medio": float(v.mean()), "massimo": float(v.max()), "nframe": int(m.sum()),
                   "t_massimo_ns": float(r["t"][m][np.argmax(v)] / 1000.0)}
        scrivi(f"  {e:8s} medio {v.mean():.4f} A   massimo {v.max():.4f} A "
               f"(t={core[e]['t_massimo_ns']:.1f} ns)   nframe {m.sum()}")
    scrivi("")

    with open(os.path.join(QUI, "verifica.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(righe) + "\n")
    dati = {
        "script": {"path": script, "sha256": sha256(script)},
        "versioni": {"python": platform.python_version(), "numpy": np.__version__,
                     "MDAnalysis": mda.__version__},
        "file": {os.path.relpath(p, RADICE): sha256(p)
                 for p in list(GRO.values()) + [r["xtc"] for r in risultati]},
        "traiettorie": {r["etichetta"]: {"nframe": r["nframe"], "t0_ns": float(r["t"][0] / 1000),
                                         "tfin_ns": float(r["t"][-1] / 1000),
                                         "maxCN_A": float(r["dcn"].max()),
                                         "box_triclino": r["box_triclino"]} for r in risultati},
        "riferimenti": {"core_rmsd_X_N6": core_XN6, "D6": D6, "D6_inverso": D6_inv,
                        "maxCN_conf_xray": dcnX, "maxCN_conf_nmr": dcnN},
        "popolazioni": pop,
        "core_rmsd_propria_conf_20_400": core,
    }
    with open(os.path.join(QUI, "verifica.json"), "w", encoding="utf-8") as f:
        json.dump(dati, f, indent=2, ensure_ascii=False)
    print("\n".join(righe))


if __name__ == "__main__":
    sys.exit(main())
