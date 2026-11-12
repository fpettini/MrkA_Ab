"""DESCRITTIVO PREVISTO DAI CRITERI (sezione "Analisi descrittive": interfaccia, cioe' superficie sepolta, legami H e ponti
salini, con occupanza per replica), eseguito il 12 nov 2026 dopo i risultati finali. Nessun esito, nessun test.
Le definizioni geometriche qui sotto NON erano nei criteri: sono fissate in questo script prima di guardare i numeri.

Legame H MrkA-Fab: donatore N o O con un H legato, accettore O oppure N senza H; distanza donatore-accettore <= 3.5 A e
  angolo H-donatore-accettore <= 30 gradi (criterio geometrico di gmx hbond). Occupanza per coppia di residui.
Ponte salino: N di Lys (NZ) o di Arg (NE, NH1, NH2) entro 4.0 A da un O carbossilico di Asp (OD1, OD2) o Glu (OE1, OE2).
Superficie sepolta: Shrake-Rupley, sonda 1.4 A, 200 punti per atomo, soli atomi pesanti con raggi di Bondi (C 1.70, N 1.55,
  O 1.52, S 1.80 A); BSA = SASA(MrkA) + SASA(Fab) - SASA(complesso), somma sui due lati; un fotogramma per ns.
Finestre: 50-500 ns (e 100-500 ns) come lo scopo 1. Riferimento: mrka-fab_completo.pdb (cristallo completato), dove per i
  legami H e i ponti salini si puo' dare solo la distanza fra gli atomi pesanti (il modello cristallografico non ha H).

Uso: python3 d9_interfaccia.py <BASE con trj_tutti/ e mrka-fab_completo.pdb> <cartella degli script di analisi> [processi]
Scrive BASE/descrittivo/posthoc12/d9_interfaccia.json e .txt"""
import os, sys, json, warnings, numpy as np
warnings.filterwarnings("ignore")
BASE = sys.argv[1]; sys.path.insert(0, sys.argv[2]); NPROC = int(sys.argv[3]) if len(sys.argv) > 3 else 4
import comune_fab as C
from scipy.spatial import cKDTree
PDB = os.path.join(BASE, "mrka-fab_completo.pdb"); TRJ = os.path.join(BASE, "trj_tutti")
D_HB, ANG_HB, D_SB, SONDA, NPT = 3.5, 30.0, 4.0, 1.4, 200
RAGGI = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
CARB = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
BASI = {"LYS": ("NZ",), "ARG": ("NE", "NH1", "NH2")}

def base3(rn):
    """Nome a tre lettere senza prefissi N/C dei terminali e con CYX/HIE ecc. lasciati come sono."""
    return rn[1:] if len(rn) == 4 and rn[0] in "NC" else rn

def sfera(n):
    i = np.arange(n) + 0.5
    ph = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)]
SF = sfera(NPT)

def sasa_atomi(X, rad, quali, albero=None):
    """SASA (A^2) degli atomi 'quali' (indici in X) nell'insieme X (N,3) con raggi rad (N)."""
    t = albero if albero is not None else cKDTree(X)
    out = np.empty(len(quali)); rmax = rad.max()
    for q, i in enumerate(quali):
        ri = rad[i] + SONDA
        vic = [j for j in t.query_ball_point(X[i], ri + rmax + SONDA) if j != i]
        P = X[i] + ri * SF
        if vic:
            vic = np.array(vic)
            d2 = ((P[:, None, :] - X[vic][None, :, :]) ** 2).sum(-1)
            lib = (d2 >= ((rad[vic] + SONDA) ** 2)[None, :]).all(1)
        else:
            lib = np.ones(NPT, bool)
        out[q] = 4 * np.pi * ri * ri * lib.mean()
    return out

def bsa(X, rad, ia, ifab):
    """BSA = somma sugli atomi dei due lati di (SASA da solo - SASA nel complesso); cambia solo vicino all'interfaccia."""
    ta = cKDTree(X[ia]); tf = cKDTree(X[ifab]); lim = 2 * (rad.max() + SONDA)
    va = np.array(sorted({i for l in tf.query_ball_tree(ta, lim) for i in l}), dtype=int)     # atomi di A vicini al Fab
    vf = np.array(sorted({i for l in ta.query_ball_tree(tf, lim) for i in l}), dtype=int)
    if len(va) == 0:
        return 0.0, 0.0, 0.0
    tutti = np.r_[ia, ifab]; Xt = X[tutti]; rt = rad[tutti]; tt = cKDTree(Xt)
    sa_soli = sasa_atomi(X[ia], rad[ia], va, ta); sa_cpx = sasa_atomi(Xt, rt, va, tt)
    sf_soli = sasa_atomi(X[ifab], rad[ifab], vf, tf); sf_cpx = sasa_atomi(Xt, rt, len(ia) + vf, tt)
    a = float((sa_soli - sa_cpx).sum()); f = float((sf_soli - sf_cpx).sum())
    return a + f, a, f

def prepara(top):
    """Atomi polari, H legati (dalla geometria del gro), cariche, pesanti; per catena."""
    import MDAnalysis as mda
    u = mda.Universe(os.path.join(TRJ, "proteina_complesso.gro")); X0 = u.atoms.positions.astype(float)
    el = np.array([n.lstrip("0123456789")[0] for n in top.nomi])
    pes = np.where(el != "H")[0]; idr = np.array([f"{c}{k}" for c, k in zip(top.catena, top.resid)])
    tp = cKDTree(X0[pes]); legato = {}
    for h in np.where(el == "H")[0]:
        d, j = tp.query(X0[h]); p = pes[j]
        if d > 1.3 or idr[p] != idr[h]:
            raise SystemExit(f"H {h} senza atomo pesante legato (d {d:.2f})")
        legato.setdefault(int(p), []).append(int(h))
    don = np.array(sorted(p for p in legato if el[p] in "NO"), dtype=int)
    acc = np.array([i for i in pes if el[i] == "O" or (el[i] == "N" and i not in legato)], dtype=int)
    return el, pes, legato, don, acc

def lavora(arg):
    r, = arg
    import MDAnalysis as mda
    gro = os.path.join(TRJ, "proteina_complesso.gro"); xtc = os.path.join(TRJ, f"complesso_{r}_0-500ns.xtc")
    top = C.Topologia(gro, "complesso", PDB); piano, gruppi = C.piano_riassemblaggio(top)
    el, pes, legato, don, acc = prepara(top)
    inA = top.catena == "A"
    donA, donF = don[inA[don]], don[~inA[don]]; accA, accF = acc[inA[acc]], acc[~inA[acc]]
    rb = np.array([base3(x) for x in top.resnomi])
    carb = np.array([i for i in pes if rb[i] in CARB and top.nomi[i] in CARB[rb[i]]], dtype=int)
    basi = np.array([i for i in pes if rb[i] in BASI and top.nomi[i] in BASI[rb[i]]], dtype=int)
    rad = np.array([RAGGI[e] for e in el[pes]]); pa = np.where(inA[pes])[0]; pf = np.where(~inA[pes])[0]
    rid = lambda i: (str(top.catena[i]), int(top.resid[i]))
    u = mda.Universe(gro, xtc); nf = u.trajectory.n_frames
    T = np.empty(nf); hb = []; sb = []; B = []
    for a in range(0, nf, 100):
        b = min(nf, a + 100)
        X = np.empty((b - a, top.n, 3)); M = np.empty((b - a, 3, 3))
        for j, ts in enumerate(u.trajectory[a:b]):
            X[j] = ts.positions; M[j] = ts.triclinic_dimensions; T[a + j] = ts.time
        Y = C.riassembla(X, M, piano, gruppi)
        for j in range(b - a):
            y = Y[j]; f = a + j; visti = set()
            for D, A in ((donA, accF), (donF, accA)):
                t = cKDTree(y[A])
                for qi, js in enumerate(t.query_ball_point(y[D], D_HB)):
                    d = D[qi]
                    for jj in js:
                        ac = A[jj]; v = y[ac] - y[d]; nv = np.linalg.norm(v)
                        for h in legato[int(d)]:
                            w = y[h] - y[d]
                            if np.degrees(np.arccos(np.clip(v @ w / (nv * np.linalg.norm(w)), -1, 1))) <= ANG_HB:
                                ra, rf = (rid(d), rid(ac)) if inA[d] else (rid(ac), rid(d))
                                visti.add((ra[1], rf[0], rf[1], str(top.nomi[d]), str(top.nomi[ac]), "A" if inA[d] else "F"))
                                break
            hb += [(f,) + v for v in visti]
            vs = set()
            for C1, C2 in ((carb[inA[carb]], basi[~inA[basi]]), (basi[inA[basi]], carb[~inA[carb]])):
                if len(C1) and len(C2):
                    t = cKDTree(y[C2])
                    for qi, js in enumerate(t.query_ball_point(y[C1], D_SB)):
                        for jj in js:
                            vs.add((rid(C1[qi])[1],) + rid(C2[jj]))
            sb += [(f,) + v for v in vs]
            if f % 10 == 0:
                B.append((f,) + bsa(y[pes], rad, pa, pf))
    t = np.round((T - T[0]) / 1000.0, 4)
    return r, t.tolist(), hb, sb, B

def cristallo():
    """Dal pdb completato (senza H): distanze N/O-N/O minime per coppia di residui, ponti salini, BSA."""
    xyz, rn = C.leggi_pdb(PDB)
    at = [(c, k, n) for c in ("A", "L", "H") for (k, n) in xyz[c]]
    X = np.array([xyz[c][(k, n)] for c, k, n in at]); el = np.array([n.lstrip("0123456789")[0] for _, _, n in at])
    rad = np.array([RAGGI[e] for e in el]); ia = np.array([i for i, a in enumerate(at) if a[0] == "A"]); iff = np.array([i for i, a in enumerate(at) if a[0] != "A"])
    tot, a, f = bsa(X, rad, ia, iff)
    pol = {}; sal = {}
    for i in ia:
        for j in iff:
            d = np.linalg.norm(X[i] - X[j])
            if d <= 4.5 and el[i] in "NO" and el[j] in "NO":
                k = (at[i][1], at[j][0], at[j][1]); pol[k] = min(pol.get(k, 9.0), float(d))
            ri, rj = rn["A"][at[i][1]], rn[at[j][0]][at[j][1]]
            ci = ri in CARB and at[i][2] in CARB[ri]; bi = ri in BASI and at[i][2] in BASI[ri]
            cj = rj in CARB and at[j][2] in CARB[rj]; bj = rj in BASI and at[j][2] in BASI[rj]
            if (ci and bj) or (bi and cj):
                k = (at[i][1], at[j][0], at[j][1]); sal[k] = min(sal.get(k, 99.0), float(d))
    return {"bsa": tot, "bsa_A": a, "bsa_Fab": f, "polari": pol, "salini": sal, "rn": rn}

if __name__ == "__main__":
    from multiprocessing import Pool
    cr = cristallo(); rn = cr["rn"]
    with Pool(NPROC) as p:
        R = {r: (t, hb, sb, B) for r, t, hb, sb, B in p.imap_unordered(lavora, [(f"run{i}",) for i in range(1, 6)])}
    reps = [f"run{i}" for i in range(1, 6)]
    nome = lambda k: f"{rn['A'][k[0]].capitalize()}{k[0]}-{k[1]} {rn[k[1]][k[2]].capitalize()}{k[2]}"
    ris = {"definizioni": {"D_HB": D_HB, "ANG_HB": ANG_HB, "D_SB": D_SB, "sonda": SONDA, "punti": NPT, "raggi": RAGGI},
           "cristallo": {"bsa": cr["bsa"], "bsa_MrkA": cr["bsa_A"], "bsa_Fab": cr["bsa_Fab"],
                         "polari": {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in cr["polari"].items()},
                         "salini": {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in cr["salini"].items()}}, "scarti": {}}
    L = [f"Interfaccia MrkA-Fab, descrittivo previsto dai criteri (12 nov 2026). Legame H: D-A <= {D_HB} A e H-D-A <= {ANG_HB:.0f} gradi; "
         f"ponte salino: N(Lys/Arg)-O(Asp/Glu) <= {D_SB} A; BSA: Shrake-Rupley, sonda {SONDA} A, atomi pesanti, raggi di Bondi, somma dei due lati.",
         f"Cristallo completato: BSA {cr['bsa']:.0f} A^2 (MrkA {cr['bsa_A']:.0f}, Fab {cr['bsa_Fab']:.0f})"]
    for sc in (50.0, 100.0):
        S = {"hb": {}, "hb_atomi": {}, "sb": {}, "bsa": {}}
        for r in reps:
            t, hb, sb, B = R[r]; t = np.array(t); w = (t >= sc - 1e-6) & (t <= 500.0 + 1e-6); n = int(w.sum()); ok = set(np.where(w)[0].tolist())
            ch = {}; ca = {}
            for f, ka, cf, kf, nd, na, lato in hb:
                if f in ok:
                    ch.setdefault((ka, cf, kf), set()).add(f); ca.setdefault((ka, cf, kf, nd, na, lato), set()).add(f)
            for k, v in ch.items(): S["hb"].setdefault(k, {})[r] = len(v) / n
            for k, v in ca.items(): S["hb_atomi"].setdefault(k, {})[r] = len(v) / n
            cs = {}
            for f, ka, cf, kf in sb:
                if f in ok: cs.setdefault((ka, cf, kf), set()).add(f)
            for k, v in cs.items(): S["sb"].setdefault(k, {})[r] = len(v) / n
            bb = np.array([x for x in B if x[0] in ok]); S["bsa"][r] = {"media": float(bb[:, 1].mean()), "ds": float(bb[:, 1].std(ddof=1)), "min": float(bb[:, 1].min()),
                                                                     "max": float(bb[:, 1].max()), "MrkA": float(bb[:, 2].mean()), "Fab": float(bb[:, 3].mean()), "n": int(len(bb))}
        ris["scarti"][f"{sc:.0f}"] = {"hb": {f"{k[0]}|{k[1]}|{k[2]}": [v.get(r, 0.0) for r in reps] for k, v in S["hb"].items()},
                                      "hb_atomi": {"|".join(map(str, k)): [v.get(r, 0.0) for r in reps] for k, v in S["hb_atomi"].items()},
                                      "sb": {f"{k[0]}|{k[1]}|{k[2]}": [v.get(r, 0.0) for r in reps] for k, v in S["sb"].items()}, "bsa": S["bsa"]}
        L += ["", f"=== scarto {sc:.0f} ns (finestra {sc:.0f}-500 ns)",
              "BSA (A^2) per replica, media (DS sui fotogrammi, min-max; un fotogramma per ns): " + "; ".join(f"{r} {S['bsa'][r]['media']:.0f} ({S['bsa'][r]['ds']:.0f}, {S['bsa'][r]['min']:.0f}-{S['bsa'][r]['max']:.0f})" for r in reps),
              "legami H per coppia di residui con occupanza >= 0.10 in almeno una replica: run1..run5 | distanza N/O-N/O minima nel cristallo (A)"]
        for k, v in sorted(S["hb"].items(), key=lambda x: -np.mean([x[1].get(r, 0.0) for r in reps])):
            f = [v.get(r, 0.0) for r in reps]
            if max(f) >= 0.10:
                at = sorted(((kk, vv) for kk, vv in S["hb_atomi"].items() if kk[:3] == k), key=lambda x: -np.mean([x[1].get(r, 0.0) for r in reps]))[0][0]
                dc = cr["polari"].get(k)
                L.append(f"  {nome(k):28s} " + " ".join(f"{x:.2f}" for x in f) + f" | media {np.mean(f):.2f} | {'%.2f' % dc if dc is not None else '> 4.5'} | piu' frequente: "
                         + (f"{at[3]} (MrkA) -> {at[4]} (Fab)" if at[5] == "A" else f"{at[3]} (Fab) -> {at[4]} (MrkA)"))
        L.append("ponti salini con occupanza >= 0.10 in almeno una replica: run1..run5 | distanza N-O minima nel cristallo (A)")
        for k, v in sorted(S["sb"].items(), key=lambda x: -np.mean([x[1].get(r, 0.0) for r in reps])):
            f = [v.get(r, 0.0) for r in reps]
            if max(f) >= 0.10:
                dc = cr["salini"].get(k)
                L.append(f"  {nome(k):28s} " + " ".join(f"{x:.2f}" for x in f) + f" | media {np.mean(f):.2f} | {'%.2f' % dc if dc is not None else '-'}")
    out = os.path.join(BASE, "descrittivo", "posthoc12")
    json.dump(ris, open(os.path.join(out, "d9_interfaccia.json"), "w"), indent=1)
    open(os.path.join(out, "d9_interfaccia.txt"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))
