"""Scopo 3: effetto del legame sulla mobilita' (criteri, sezione "Scopo 3: effetto del legame sulla mobilita'").

Uso (dalla cartella complesso_fab, dopo 11_controlli.py):
    python3 analisi/14_scopo3_rmsf.py [--cache analisi/cache] [--out analisi/risultati]

RMSF dei Ca della catena A dopo fit sui 43 Ca di strand (riferimento = cristallo), per replica, finestra dallo
scarto (50 ns; 100 di sensibilita') a T3; il libero e' tagliato alla stessa lunghezza del complesso.
Unita' statistica: media dell'RMSF sui residui del segmento, una per replica. 87-101: legato < libero, Welch
unilaterale; 138-147: Welch bilaterale; Holm sui due test, alfa 0.05; sensibilita' Mann-Whitney esatto.
Esito positivo (differenza) solo se regge sia con scarto 50 sia con 100 ns, altrimenti "non robusto".
I test si fanno SOLO nell'analisi finale dello scopo 3 (tutte le repliche del disegno alla durata confermativa, o
--finale in 11_controlli.py): nell'analisi intermedia lo script riporta solo medie e RMSF per residuo (descrittivi).
Scrive OUT/scopo3_rmsf.txt, .json e scopo3_rmsf_per_residuo.csv.
"""
import os, json, csv, argparse
import numpy as np
import comune_fab as C

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="analisi/cache")
ap.add_argument("--out", default="analisi/risultati")
arg = ap.parse_args()

ctrl = C.leggi_controlli(arg.out)
serie = C.leggi_cache(arg.cache)
T3 = ctrl["T3_ns"]; stato = ctrl["stato"]["scopo3"]
ALFA = 0.05


def rmsf(d, w):
    X = d["ca_fit"][w].astype(np.float64)
    return np.sqrt(((X - X.mean(0)) ** 2).sum(-1).mean(0))


def holm(p):
    o = np.argsort(p); m = len(p); adj = np.empty(m); cur = 0.0
    for i, j in enumerate(o):
        cur = max(cur, min(1.0, (m - i) * p[j])); adj[j] = cur
    return adj


ris = {}
for sc in C.SCARTI:
    vc = C.valide(ctrl, "complesso", "complesso", sc)["valide"]
    vl = C.valide(ctrl, "libero", "scopo3", sc)["valide"]
    R = {}
    for s, reps in (("complesso", vc), ("libero", vl)):
        for r in reps:
            d = serie[(s, r)]
            R[(s, r)] = rmsf(d, C.finestra(d["t_ns"], sc, T3))
    resid = serie[next(iter(R))]["ca_resid"] if R else np.zeros(0, int)
    seg = {}
    for nome, res in C.SEGMENTI.items():
        m = np.isin(resid, res)
        seg[nome] = {s: {r: float(R[(s, r)][m].mean()) for r in reps} for s, reps in (("complesso", vc), ("libero", vl))}
    test = None
    if stato == "FINALE" and len(vc) >= 2 and len(vl) >= 2:
        from scipy import stats
        test = {}
        for nome, alt in (("87-101", "less"), ("138-147", "two-sided")):
            a = np.array(list(seg[nome]["complesso"].values())); b = np.array(list(seg[nome]["libero"].values()))
            w = stats.ttest_ind(a, b, equal_var=False, alternative=alt)
            mw = stats.mannwhitneyu(a, b, alternative=alt, method="exact")
            test[nome] = {"alternativa": alt, "welch_t": float(w.statistic), "welch_p": float(w.pvalue), "mw_U": float(mw.statistic), "mw_p": float(mw.pvalue)}
        ph = holm(np.array([test[k]["welch_p"] for k in test])); pm = holm(np.array([test[k]["mw_p"] for k in test]))
        for i, k in enumerate(test):
            test[k]["welch_p_holm"] = float(ph[i]); test[k]["mw_p_holm"] = float(pm[i])
            test[k]["differenza"] = bool(ph[i] < ALFA); test[k]["differenza_mw"] = bool(pm[i] < ALFA)
    ris[sc] = {"complesso_valide": vc, "libero_valide": vl, "segmenti": seg, "test": test,
               "rmsf": {f"{s}_{r}": v.tolist() for (s, r), v in R.items()}, "resid": [int(x) for x in resid]}

fin = None
if all(ris[sc]["test"] for sc in C.SCARTI):
    fin = {}
    for nome in C.SEGMENTI:
        e = [ris[sc]["test"][nome]["differenza"] for sc in C.SCARTI]
        fin[nome] = ("differenza" if e[0] else "nessuna differenza") if e[0] == e[1] else "non robusto"

P = C.prov(__file__, [v["_file"] for v in serie.values()] + [os.path.join(arg.out, "controlli.json")],
           {"cache": arg.cache, "out": arg.out}, {"CACHE": arg.cache, "OUT": arg.out})
L = C.testo_prov(P) + C.avviso_intermedia(stato)
L += ["", f"Scopo 3 ({stato}): RMSF dei Ca (fit sui 43 Ca di strand), finestra dallo scarto a T3 = {T3:.1f} ns per entrambi i sistemi."]
for sc in C.SCARTI:
    r = ris[sc]
    L += ["", f"=== Scarto {sc:.0f} ns: complesso {r['complesso_valide']}, libero {r['libero_valide']}"]
    for nome, v in r["segmenti"].items():
        for s in ("complesso", "libero"):
            x = np.array(list(v[s].values()))
            L.append(f"  {nome:8s} {s:10s} media di segmento per replica (A): " + ", ".join(f"{k} {y:.2f}" for k, y in v[s].items())
                     + (f"  | media {x.mean():.2f}, DS {x.std(ddof=1):.2f}" if len(x) > 1 else ""))
        a = np.array(list(v["complesso"].values())); b = np.array(list(v["libero"].values()))
        if len(a) and len(b):
            L.append(f"  {nome:8s} differenza delle medie (complesso - libero): {a.mean() - b.mean():+.2f} A (descrittiva)")
    if r["test"]:
        for nome, t in r["test"].items():
            L.append(f"  test {nome} ({t['alternativa']}): Welch t {t['welch_t']:.2f}, p {t['welch_p']:.4f}, Holm {t['welch_p_holm']:.4f}; "
                     f"Mann-Whitney U {t['mw_U']:.0f}, p {t['mw_p']:.4f}, Holm {t['mw_p_holm']:.4f}")
    else:
        L.append("  test non eseguiti: analisi intermedia (i test si fanno una volta sola, sui dati finali)")
if fin:
    L += ["", "Esito (Welch + Holm, stesso a 50 e 100 ns): " + "; ".join(f"{k}: {v}" for k, v in fin.items())]
open(os.path.join(arg.out, "scopo3_rmsf.txt"), "w").write("\n".join(L) + "\n")
with open(os.path.join(arg.out, "scopo3_rmsf_per_residuo.csv"), "w", newline="") as f:
    w = csv.writer(f)
    for sc in C.SCARTI:
        r = ris[sc]
        cols = list(r["rmsf"])
        w.writerow([f"scarto_{int(sc)}ns", "resid_wt"] + cols)
        for i, k in enumerate(r["resid"]):
            w.writerow(["", k] + [f"{r['rmsf'][c][i]:.3f}" for c in cols])
json.dump({"provenienza": P, "stato": stato, "T3_ns": T3, "per_scarto": {str(int(k)): v for k, v in ris.items()}, "esito": fin},
          open(os.path.join(arg.out, "scopo3_rmsf.json"), "w"), indent=1)
print("\n".join(L))
