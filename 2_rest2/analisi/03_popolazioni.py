"""Popolazioni X-ray / NMR dei loop di MrkA nei due run REST2 (r0) e confronto fra i run.

Riferimento NMR principale: il modello 6 di 9HW9, medoide, con la sequenza X-ray (rest2/nmr/conf.gro).
Controllo: il piu' vicino dei 20 modelli di 9HW9.

Uso (dalla cartella del progetto mrkA):
    python analisi_rest2/03_popolazioni.py . analisi_rest2/cache analisi_rest2/risultati [--t0 NS]
Legge le cache di 01_estrai.py (un sistema senza cache viene saltato). Scrive OUT/riferimenti.txt,
OUT/distanze_<sistema>.npz, OUT/popolazioni.txt e OUT/popolazioni.json. Scarto iniziale fisso di 20 ns
(decisione del 24 ottobre 2026); --t0 lo cambia. Sensibilita' con scarto di 120 ns (--t0-sens, 25 ottobre):
un loop e' convergente solo se lo e' con entrambi gli scarti. Sensibilita' alla soglia "altro" (25 ottobre):
assegnazione principale e di controllo anche con 0.75*D e 1.25*D. pX e pN sono frazioni di frame entro la
soglia, non la conformazione dominante; per 87-101 e 112-118 conta solo il confronto fra i run.
Criterio (25 ottobre): test di equivalenza |D| + t*SEc < 0.15 per pX e pN, fra i run e fra le meta' di ogni run;
t = quantile 97.5% di Student a n-1 gradi di liberta' (n = blocchi). Blocchi da 20 ns se la finestra e' sotto
300 ns, da 40 ns sopra (--blocco li fissa); il SE e' il massimo fra il blocco scelto e i blocchi piu' lunghi che
lasciano almeno 4 blocchi, cosi' si prende dove smette di crescere. I t0 di pymbar si stampano solo come diagnostica: con le
fluttuazioni lente di questi loop la stima e' inaffidabile (sull'X-ray da' fino a 181 ns).
"""
import sys, os, json, argparse, warnings, numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from scipy import stats
from pymbar import timeseries
from provenienza import prov, testo, provenienza_cache

ap = argparse.ArgumentParser()
ap.add_argument("root"); ap.add_argument("cache"); ap.add_argument("out")
ap.add_argument("--t0", type=float, default=None, help="scarto iniziale, ns (default 20)")
ap.add_argument("--t0-sens", type=float, default=120.0, help="scarto della sensibilita', ns (default 120; 0 = niente)")
ap.add_argument("--blocco", type=float, default=None, help="blocchi per il SE, ns (default 20 sotto 300 ns di finestra, 40 sopra)")
arg = ap.parse_args()
ROOT, CACHE, OUT = arg.root, arg.cache, arg.out
SISTEMI = ["xray", "nmr"]
CORE = [44, 45, 46, 47, 78, 79, 80, 81, 82, 83, 84, 85, 88, 94, 102, 103, 104, 105, 107, 108, 122, 134,
        135, 136, 137, 138, 148, 162, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 195, 196, 197, 198, 199]
LOOP = {"59-77": range(59, 78), "87-101": range(87, 102),
        "112-118": [k for k in range(112, 119) if k not in (116, 117)],
        "124-133": range(124, 134), "138-147": range(138, 148), "159-169": range(159, 170),
        "180-192": range(180, 193), "48-56": range(48, 57)}
LOOP = {k: list(v) for k, v in LOOP.items()}
BBN = ["N", "CA", "C", "O"]
NMR_PDB = os.path.join(ROOT, "9hw9_20modelli.pdb")
OFFSET = 27                          # wt = costrutto + 27
T_MIN = 20.0                         # ns, scarto iniziale fisso
L_SE = [5, 10, 20, 30, 40, 60, 80]
MARGINE = 0.15                       # margine di equivalenza sulle popolazioni
MIN_BLOCCHI = 4                      # blocchi piu' lunghi del scelto contano solo se ne restano almeno 4
CLASSI = ["X", "N", "altro"]
FATTORI_D = [0.75, 1.25]             # sensibilita' alla soglia "altro" (oltre a 1.0)

def indici(names, resids, lista, nomi):
    """Indici degli atomi (residuo per residuo, nell'ordine di nomi)."""
    out = []
    for k in lista:
        for n in nomi:
            w = np.where((names == n) & (resids == k))[0]
            if len(w) != 1:
                raise SystemExit(f"atomo {n} del residuo {k} non trovato o duplicato")
            out.append(w[0])
    return np.array(out)

def kabsch(P, Q):
    """Rotazione e centri che sovrappongono ogni frame di P (F,N,3) su Q (N,3)."""
    Pm = P.mean(1, keepdims=True); Qm = Q.mean(0)
    H = np.einsum('fni,nj->fij', P - Pm, Q - Qm)
    U, S, Vt = np.linalg.svd(H)
    D = np.tile(np.eye(3), (len(P), 1, 1)); D[:, 2, 2] = np.sign(np.linalg.det(U @ Vt))
    return U @ D @ Vt, Pm, Qm

def dist(Pc, PL, Qc, QL):
    """d(frame, R): fit sui Ca del core (Pc su Qc), poi RMSD degli atomi del loop senza nuovo fit."""
    R, Pm, Qm = kabsch(Pc, Qc)
    A = np.einsum('fni,fij->fnj', PL - Pm, R) + Qm
    return np.sqrt(((A - QL) ** 2).sum(2).mean(1))

def rmsd_core(Pc, Qc):
    R, Pm, Qm = kabsch(Pc[None], Qc)
    A = np.einsum('fni,fij->fnj', Pc[None] - Pm, R) + Qm
    return np.sqrt(((A - Qc) ** 2).sum(2).mean(1))[0]

def geometrie(pos, names, resids):
    """Dizionario {loop: (K,nl,3)} piu' 'core' (K,43,3) da posizioni (K,N,3)."""
    ic = indici(names, resids, CORE, ["CA"])
    g = {"core": pos[:, ic].astype(np.float64)}
    for l, rs in LOOP.items():
        g[l] = pos[:, indici(names, resids, rs, BBN)].astype(np.float64)
    return g

def carica_riferimenti():
    ux = mda.Universe(os.path.join(ROOT, "rest2", "xray", "conf.gro"))
    gx = geometrie(ux.atoms.positions[None], ux.atoms.names, ux.atoms.resids)
    u6 = mda.Universe(os.path.join(ROOT, "rest2", "nmr", "conf.gro"))
    g6 = geometrie(u6.atoms.positions[None], u6.atoms.names, u6.atoms.resids)
    un = mda.Universe(NMR_PDB)
    pos = np.array([un.atoms.positions.copy() for _ in un.trajectory])
    gn = geometrie(pos, un.atoms.names, un.atoms.resids + OFFSET)
    return gx, gn, g6

def se_blocchi(y, dt, L, con_n=False):
    """SE dai blocchi di L ns contati dalla fine (solo blocchi interi, l'avanzo si toglie all'inizio)."""
    bl = int(round(L / dt)); n = len(y) // bl
    if n < 2:
        return (np.nan, n) if con_n else np.nan
    m = y[len(y) - n * bl:].reshape(n, bl).mean(1)
    se = m.std(ddof=1) / np.sqrt(n)
    return (se, n) if con_n else se

def se_usato(y, dt, blocco):
    """SE al blocco scelto, o il piu' grande fra i blocchi piu' lunghi con almeno MIN_BLOCCHI blocchi. -> (SE, L, n)"""
    se, n = se_blocchi(y, dt, blocco, True); best = (se, blocco, n)
    for L in L_SE:
        if L > blocco:
            s2, n2 = se_blocchi(y, dt, L, True)
            if n2 >= MIN_BLOCCHI and s2 > best[0]:
                best = (s2, L, n2)
    return best

def equivalenza(a, b, c):
    """Test |D| + t*SEc < MARGINE fra due stime (a, b con p, se, n della classe c)."""
    d = a["p" + c] - b["p" + c]; sc = float(np.hypot(a["se" + c], b["se" + c]))
    gl = min(a["n" + c], b["n" + c]) - 1
    if np.isnan(sc) or gl < 1:
        return {"delta": float(d), "SEc": sc, "gl": int(gl), "t": np.nan, "limite": np.nan, "equivalenti": None, "diversi": None}
    tq = float(stats.t.ppf(0.975, gl))
    return {"delta": float(d), "SEc": sc, "gl": int(gl), "t": tq, "limite": float(abs(d) + tq * sc),
            "equivalenti": bool(abs(d) + tq * sc < MARGINE), "diversi": bool(abs(d) > tq * sc)}

def esito_di(prove):
    """convergente se tutte le prove sono equivalenti; differenza rilevata se almeno una prova e' significativa E non
    equivalente (25 ott, dalla revisione: una differenza significativa ma entro il margine non conta); altrimenti non dimostrata."""
    if any(p["equivalenti"] is None for p in prove):
        return "indeterminato (SE nan o troppo pochi blocchi)"
    if all(p["equivalenti"] for p in prove):
        return "convergente"
    if any(p["diversi"] and not p["equivalenti"] for p in prove):
        return "differenza rilevata"
    return "equivalenza non dimostrata"

def stima(code, dt, blocco):
    """Frazione di ogni classe e SE a blocchi (se_usato) di una finestra; code: 0 = X, 1 = N, 2 = altro."""
    o = {}
    for j, c in enumerate(CLASSI):
        y = (code == j).astype(float)
        se, L, n = se_usato(y, dt, blocco)
        o["p" + c] = float(y.mean()); o["se" + c] = float(se); o["L" + c] = float(L); o["n" + c] = int(n)
        o["se_L" + c] = {str(L): float(se_blocchi(y, dt, L)) for L in L_SE}
    return o

def scegli_blocco(fine_ns, T0):
    return arg.blocco if arg.blocco is not None else (20.0 if fine_ns - T0 < 300 else 40.0)

def analisi(codici, T0, tabella):
    """Punti 5-7 per un tipo di assegnazione. codici[s][loop] = (t, code)."""
    ris = {"popolazioni": {}, "confronto": {}}
    for s, cs in codici.items():
        ris["popolazioni"][s] = {}
        for l, (t, code) in cs.items():
            dt = np.median(np.diff(t)) / 1000
            w = t >= T0 * 1000 - 1e-6
            tw, cw = t[w], code[w]
            blocco = scegli_blocco(tw[-1] / 1000, T0)
            o = stima(cw, dt, blocco); o["blocco_ns"] = blocco
            mid = T0 + (tw[-1] / 1000 - T0) / 2          # le due meta' uguali, per tempo
            h = []
            for m in (tw / 1000 < mid, tw / 1000 >= mid):
                h.append(stima(cw[m], dt, blocco))
            o["meta"] = [{k: v for k, v in x.items() if k in ("pX", "seX", "nX", "LX", "pN", "seN", "nN", "LN")} for x in h]
            o["meta_dettaglio"] = {c: equivalenza(h[0], h[1], c) for c in ("X", "N")}
            o["esito_meta"] = esito_di(o["meta_dettaglio"].values())
            o["n_frame"] = int(w.sum()); o["t_inizio_ns"] = float(tw[0] / 1000); o["t_fine_ns"] = float(tw[-1] / 1000)
            ris["popolazioni"][s][l] = o
    if all(s in ris["popolazioni"] for s in SISTEMI):
        fx, fn = (list(ris["popolazioni"][s].values())[0]["t_fine_ns"] for s in SISTEMI)
        if abs(fx - fn) > 1:
            print(f"ATTENZIONE: i run hanno lunghezze diverse (xray {fx:g} ns, nmr {fn:g} ns): il confronto vale solo a pari lunghezza")
        for l in LOOP:
            a, b = ris["popolazioni"]["xray"][l], ris["popolazioni"]["nmr"][l]
            c = {k: equivalenza(a, b, k) for k in ("X", "N")}
            prove = [c["X"], c["N"]] + list(a["meta_dettaglio"].values()) + list(b["meta_dettaglio"].values())
            c["esito_run"] = esito_di([c["X"], c["N"]])
            c["esito"] = esito_di(prove)
            c["n_prove"] = len(prove)
            ris["confronto"][l] = c
    return ris

def scrivi(ris, titolo, T0, righe):
    P = ris["popolazioni"]
    righe.append(f"\n{'=' * 100}\n{titolo}\n{'=' * 100}")
    for s in P:
        righe.append(f"\n--- sistema {s}: finestra [{T0:g}, {list(P[s].values())[0]['t_fine_ns']:g}] ns ---")
        righe.append(f"{'loop':8s} {'frame':>6s} | {'pX':>6s} {'SE':>6s} | {'pN':>6s} {'SE':>6s} | {'pAltro':>6s} {'SE':>6s}")
        for l, o in P[s].items():
            righe.append(f"{l:8s} {o['n_frame']:6d} | {o['pX']:6.3f} {o['seX']:6.3f} | {o['pN']:6.3f} {o['seN']:6.3f} | "
                         f"{o['paltro']:6.3f} {o['sealtro']:6.3f}")
        o0 = list(P[s].values())[0]
        righe.append(f"  SE usato = blocchi da {o0['blocco_ns']:g} ns o piu' lunghi se il SE cresce (almeno {MIN_BLOCCHI} blocchi); "
                     f"L/n = blocco usato e numero di blocchi")
        righe.append("  SE con blocchi da " + "/".join(str(L) for L in L_SE) + " ns (pX | pN), nan = meno di 2 blocchi:")
        for l, o in P[s].items():
            f = lambda c: "/".join(f"{o['se_L' + c][str(L)]:.3f}" for L in L_SE)
            righe.append(f"   {l:8s} {f('X')} | {f('N')}   usato X {o['LX']:g}/{o['nX']} N {o['LN']:g}/{o['nN']}")
        righe.append(f"  Meta' (p(SE)); equivalenza |D| + t*SEc < {MARGINE}:")
        for l, o in P[s].items():
            m = o["meta"]; d = o["meta_dettaglio"]
            g = lambda c: (f"p{c} {m[0]['p' + c]:.3f}({m[0]['se' + c]:.3f}) {m[1]['p' + c]:.3f}({m[1]['se' + c]:.3f}) "
                           f"D {d[c]['delta']:+.3f} t {d[c]['t']:.2f} gl {d[c]['gl']} lim {d[c]['limite']:.3f}")
            righe.append(f"   {l:8s} {g('X')} | {g('N')} | {o['esito_meta']}")
    if ris["confronto"]:
        righe.append("\n--- confronto fra i run: Delta = p(xray) - p(nmr) ---")
        righe.append(f"Equivalenza |D| + t*SEc < {MARGINE}; esito = fra i run E meta' di entrambi i run (6 prove per loop)")
        righe.append(f"{'loop':8s} | {'D pX':>7s} {'SEc':>6s} {'lim':>6s} | {'D pN':>7s} {'SEc':>6s} {'lim':>6s} | {'fra i run':26s} | esito")
        for l, c in ris["confronto"].items():
            righe.append(f"{l:8s} | {c['X']['delta']:+7.3f} {c['X']['SEc']:6.3f} {c['X']['limite']:6.3f} | {c['N']['delta']:+7.3f} "
                         f"{c['N']['SEc']:6.3f} {c['N']['limite']:6.3f} | {c['esito_run']:26s} | {c['esito']}")
        nl = len(ris["confronto"])
        righe.append(f"Confronti multipli: {nl} loop x 6 prove x 2 classi = {nl * 12} test, ciascuno al 95% (nessuna correzione).")
        gmin = min(min(c[k]["gl"] for k in ("X", "N")) for c in ris["confronto"].values())
        gm = min(min(o["meta_dettaglio"][k]["gl"] for k in ("X", "N")) for P_ in ris["popolazioni"].values() for o in P_.values())
        if min(gmin, gm) < 3:
            righe.append(f"ATTENZIONE: test poco informativo, gradi di liberta' minimi {min(gmin, gm)} "
                         f"(t = {stats.t.ppf(0.975, max(min(gmin, gm), 1)):.2f}): meta' con troppo pochi blocchi.")
    else:
        righe.append("\nConfronto fra i run non eseguito: manca un sistema.")

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    gx, gn, g6 = carica_riferimenti()
    ins = [os.path.join(CACHE, f"{s}_r0.npz") for s in SISTEMI] + [NMR_PDB] + \
          [os.path.join(ROOT, "rest2", s, "conf.gro") for s in SISTEMI]
    P = prov(__file__, ins, vars(arg), {"ROOT": ROOT, "CACHE": CACHE, "OUT": OUT}); P["cache"] = provenienza_cache(CACHE, SISTEMI, ["r0"])
    # --- 1. tabella dei riferimenti
    NM = len(gn["core"])
    ref = {}
    righe = [f"Riferimenti: {NM} modelli NMR (9HW9, wt = costrutto + {OFFSET}), X-ray, modello 6 mutato (rest2/nmr/conf.gro)",
             f"{'loop':8s} | {'X-NMR piu vicino [mod]':>24s} | {'mediana X-NMR':>13s} | {'dispers. NMR':>12s} | {'X-mod.6 mut.':>12s}"]
    for l in LOOP:
        dX = dist(gn["core"], gn[l], gx["core"][0], gx[l][0])          # (20,): d(X-ray, modello k)
        coppie = [dist(gn["core"][i + 1:], gn[l][i + 1:], gn["core"][i], gn[l][i]) for i in range(NM - 1)]
        disp = np.median(np.concatenate(coppie))
        d6 = dist(g6["core"], g6[l], gx["core"][0], gx[l][0])[0]
        k = int(np.argmin(dX))
        ref[l] = {"D": float(dX[k]), "modello": k + 1, "mediana": float(np.median(dX)), "disp": float(disp), "D6": float(d6)}
        righe.append(f"{l:8s} | {dX[k]:18.2f} [{k + 1:2d}] | {np.median(dX):13.2f} | {disp:12.2f} | {d6:12.2f}")
    assert sum(len(c) for c in coppie) == NM * (NM - 1) // 2          # 190 coppie
    rc = rmsd_core(g6["core"][0], gx["core"][0])
    righe.append(f"\nRMSD dei {len(CORE)} Ca del CORE fra X-ray e rest2/nmr/conf.gro: {rc:.2f} A")
    open(os.path.join(OUT, "riferimenti.txt"), "w").write("\n".join(testo(P) + righe) + "\n")
    print("\n".join(righe))

    # --- 2-3. distanze per frame e assegnazione
    presenti = [s for s in SISTEMI if os.path.exists(os.path.join(CACHE, f"{s}_r0.npz"))]
    for s in SISTEMI:
        if s not in presenti:
            print(f"\nsistema {s}: cache {s}_r0.npz assente, salto")
    if not presenti:
        raise SystemExit("nessun sistema da analizzare")
    cod = {k: {} for k in ["main", "ctrl"] + [f"{n}_{f:g}" for n in ("main", "ctrl") for f in FATTORI_D]}; T = {}; t0tab = []
    for s in presenti:
        d = np.load(os.path.join(CACHE, f"{s}_r0.npz"), allow_pickle=True)
        X, t = d["X"], d["t"]; names, resids = d["names"].astype(str), d["resids"]
        g = geometrie(X, names, resids)
        T[s] = t
        for k in cod:
            cod[k][s] = {}
        sav = {"t": t}
        for l in LOOP:
            dXf = dist(g["core"], g[l], gx["core"][0], gx[l][0])
            dk = np.stack([dist(g["core"], g[l], gn["core"][k], gn[l][k]) for k in range(NM)], 1)
            dN, iN = dk.min(1), dk.argmin(1) + 1
            dN6 = dist(g["core"], g[l], g6["core"][0], g6[l][0])
            sav.update({f"dX_{l}": dXf, f"dN_{l}": dN, f"iN_{l}": iN, f"dN6_{l}": dN6})
            for nome, dn, D in (("main", dN6, ref[l]["D6"]), ("ctrl", dN, ref[l]["D"])):
                code = np.where(np.minimum(dXf, dn) > D, 2, np.where(dXf < dn, 0, 1))
                cod[nome][s][l] = (t, code)
                for f in FATTORI_D:
                    cf = np.where(np.minimum(dXf, dn) > f * D, 2, np.where(dXf < dn, 0, 1))
                    cod[f"{nome}_{f:g}"][s][l] = (t, cf)
        np.savez(os.path.join(OUT, f"distanze_{s}.npz"), provenienza=np.array(json.dumps(P)), **sav)

    # --- 4. scarto iniziale (serie 0/1 a 100 ps, assegnazione principale = modello 6)
    righe = []
    righe.append("Diagnostica (non usata per lo scarto): pymbar.timeseries.detect_equilibration(A, nskip=10) su serie 0/1 a 100 ps")
    righe.append(f"{'sistema':8s} {'loop':8s} {'serie':6s} {'t0 [ns]':>8s}")
    t0max = 0.0
    for s in presenti:
        for l, (t, code) in cod["main"][s].items():
            for j, c in ((0, "X"), (1, "N")):
                A = (code == j).astype(float)[::10]
                if np.ptp(A) == 0:
                    righe.append(f"{s:8s} {l:8s} {c:6s} {'costante':>8s}"); continue
                i0 = timeseries.detect_equilibration(A, nskip=10)[0]
                t0 = i0 * 10 * np.median(np.diff(t)) / 1000
                t0max = max(t0max, t0); t0tab.append({"sistema": s, "loop": l, "serie": c, "t0_ns": float(t0)})
                righe.append(f"{s:8s} {l:8s} {c:6s} {t0:8.1f}")
    T0 = arg.t0 if arg.t0 is not None else T_MIN
    righe.append(f"\nt0 massimo {t0max:.1f} ns (solo diagnostica) -> scarto usato T0 = {T0:g} ns" +
                 (" (da --t0)" if arg.t0 is not None else " (fisso)"))
    print("\n".join(righe))

    # --- 5-8
    R = {"provenienza": P, "T0_ns": float(T0), "t0_serie": t0tab, "riferimenti": ref, "core_rmsd_xray_nmr6": float(rc),
         "sistemi": presenti}
    R["principale_modello6"] = analisi(cod["main"], T0, ref)
    R["controllo_20modelli"] = analisi(cod["ctrl"], T0, ref)
    for f in FATTORI_D:
        R[f"principale_modello6_D{f:g}"] = analisi(cod[f"main_{f:g}"], T0, ref)
        R[f"controllo_20modelli_D{f:g}"] = analisi(cod[f"ctrl_{f:g}"], T0, ref)
    TS = arg.t0_sens
    if TS and TS != T0:
        R["T0_sens_ns"] = float(TS)
        R["sensibilita_modello6"] = analisi(cod["main"], TS, ref)
        if R["principale_modello6"]["confronto"]:
            R["esito_finale"] = {}
            for l in LOOP:
                e1 = R["principale_modello6"]["confronto"][l]["esito"]; e2 = R["sensibilita_modello6"]["confronto"][l]["esito"]
                R["esito_finale"][l] = "convergente" if e1 == e2 == "convergente" else f"non convergente ({T0:g} ns: {e1}; {TS:g} ns: {e2})"
    out = testo(P) + list(righe)
    out.append("\nRiferimenti usati (D6 = d(X-ray, modello 6 mutato), D = d(X-ray, modello NMR piu' vicino)):")
    for l, r_ in ref.items():
        out.append(f"   {l:8s} D6 {r_['D6']:.2f}  D {r_['D']:.2f} [modello {r_['modello']}]")
    out.append("Regola: 'altro' se min(dX, dN) > D; altrimenti X se dX < dN, N se no (principale: dN6 e D6).")
    scrivi(R["principale_modello6"], "ASSEGNAZIONE PRINCIPALE (dN6 = modello 6 mutato, D6)", T0, out)
    scrivi(R["controllo_20modelli"], "ASSEGNAZIONE DI CONTROLLO (dN = minimo sui 20 modelli NMR, D)", T0, out)
    for f in FATTORI_D:
        scrivi(R[f"principale_modello6_D{f:g}"], f"SENSIBILITA' ALLA SOGLIA: principale con {f:g}*D6", T0, out)
        scrivi(R[f"controllo_20modelli_D{f:g}"], f"SENSIBILITA' ALLA SOGLIA: controllo con {f:g}*D", T0, out)
    if R["principale_modello6"]["confronto"]:
        out.append(f"\n--- RIEPILOGO ESITI (scarto {T0:g} ns) per soglia ---")
        cols = [("princ D6", "principale_modello6")] + [(f"princ {f:g}D6", f"principale_modello6_D{f:g}") for f in FATTORI_D] + \
               [("contr D", "controllo_20modelli")] + [(f"contr {f:g}D", f"controllo_20modelli_D{f:g}") for f in FATTORI_D]
        out.append(f"{'loop':8s} " + " ".join(f"{c:>26s}" for c, _ in cols))
        for l in LOOP:
            out.append(f"{l:8s} " + " ".join(f"{R[k]['confronto'][l]['esito']:>26s}" for _, k in cols))
        out.append("Per 87-101 e 112-118 conta solo il confronto fra i run (separazione X/NMR debole).")
    if "sensibilita_modello6" in R:
        scrivi(R["sensibilita_modello6"], f"SENSIBILITA' ALLO SCARTO: assegnazione principale con T0 = {TS:g} ns", TS, out)
    if "esito_finale" in R:
        out.append(f"\n--- ESITO FINALE (assegnazione principale): convergente solo se lo e' con scarto {T0:g} e {TS:g} ns ---")
        for l, e in R["esito_finale"].items():
            out.append(f"{l:8s} {e}")
    open(os.path.join(OUT, "popolazioni.txt"), "w").write("\n".join(out) + "\n")
    json.dump(R, open(os.path.join(OUT, "popolazioni.json"), "w"), indent=1)
    print("\n".join(out[len(righe):]))
