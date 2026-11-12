"""POST HOC, 12 nov 2026: numeri di contorno chiesti dalla lettura critica dei testi finali. Nessun esito, nessun test nuovo.
(1) RMSF del libero oltre la durata confermativa (descrittivo PREVISTO dai criteri: "il libero da 500 ns a 1 us per lo scopo 3");
(2) differenza rilevabile per 138-147 (come la Nota 1, bilaterale); (3) RMSF dei residui 36-41; (4) frazione dei fotogrammi con
phi95 > -90 e phi97 < -110 (seconda disposizione); (5) correlazione di rango con lo scarto a 100 ns e dentro i due gruppi;
(6) legame Asp40-Ala90 e involucro; (7) involucro per replica piu' stretto: quante repliche del complesso ci stanno dentro;
(8) Welch bilaterale per 87-101.
Uso: python3 d11_verifiche.py <BASE con cache/, risultati/, descrittivo/> <cartella degli script di analisi>"""
import sys, json, numpy as np
from scipy import stats
W = sys.argv[1]; sys.path.insert(0, sys.argv[2])
import comune_fab as C
P = []; JS = {}
S3 = json.load(open(f"{W}/risultati/scopo3_rmsf.json")); S2 = json.load(open(f"{W}/risultati/scopo2_epitopo.json"))
k50 = [k for k in S3["per_scarto"] if float(k) == 50.0][0] if "per_scarto" in S3 else None
R3 = S3["per_scarto"][k50] if k50 else S3[[k for k in S3 if k.startswith("50")][0]]
resid = np.array(R3["resid"])
def rmsf(d, w):
    X = d["ca_fit"][w].astype(np.float64)
    return np.sqrt(((X - X.mean(0)) ** 2).sum(-1).mean(0))
cache = {}
for s, n in (("complesso", 5), ("libero", 10)):
    for r in range(1, n + 1):
        z = np.load(f"{W}/cache/serie_{s}_run{r}.npz", allow_pickle=True)
        cache[(s, r)] = {"t_ns": z["t_ns"], "ca_fit": z["ca_fit"], "ca_resid": z["ca_resid"]}
P.append("=== (1) RMSF medio di segmento (A) del libero per finestra; controllo 50-500 contro il json di 14")
seg = {"87-101": np.isin(resid, C.SEGMENTI["87-101"]), "138-147": np.isin(resid, C.SEGMENTI["138-147"])}
fin = {"50-500": (50, 500), "550-1000": (550, 1000), "500-1000": (500, 1000), "50-1000": (50, 1000)}
tab = {}
for nome, m in seg.items():
    for f, (a, b) in fin.items():
        v = []
        for r in range(1, 11):
            d = cache[("libero", r)]; w = (d["t_ns"] >= a - 1e-6) & (d["t_ns"] <= b + 1e-6)
            v.append(float(rmsf(d, w)[m].mean()))
        tab[(nome, f)] = np.array(v); JS.setdefault('rmsf_libero', {}).setdefault(nome, {})[f] = v
        if f == "50-500":
            rif = np.array([R3["segmenti"][nome]["libero"][f"run{r}"] for r in range(1, 11)])
            assert np.abs(rif - tab[(nome, f)]).max() < 1e-6, np.abs(rif - tab[(nome, f)]).max()
        P.append(f"  {nome:8s} {f:9s}: " + " ".join(f"{x:.2f}" for x in v) + f" | media {np.mean(v):.3f}, DS {np.std(v, ddof=1):.3f}, min {min(v):.2f}, max {max(v):.2f}")
    c = np.array([R3["segmenti"][nome]["complesso"][f"run{r}"] for r in range(1, 6)])
    P.append(f"  {nome:8s} complesso 50-500: " + " ".join(f"{x:.2f}" for x in c) + f" | media {c.mean():.3f}, DS {c.std(ddof=1):.3f}, min {c.min():.2f}, max {c.max():.2f}")
P.append("  fotogrammi per finestra: " + ", ".join(f"{f} {int(((cache[('libero', 1)]['t_ns'] >= a - 1e-6) & (cache[('libero', 1)]['t_ns'] <= b + 1e-6)).sum())}" for f, (a, b) in fin.items()))

P.append("\n=== (2) differenza rilevabile con potenza 80% (t non centrale, gradi di liberta' di Welch-Satterthwaite)")
def rilevabile(s1, n1, s2, n2, alfa, lati):
    se = np.sqrt(s1 ** 2 / n1 + s2 ** 2 / n2); df = se ** 4 / ((s1 ** 2 / n1) ** 2 / (n1 - 1) + (s2 ** 2 / n2) ** 2 / (n2 - 1))
    def pot(d):
        nc = d / se
        if lati == 1:
            return 1 - stats.nct.cdf(stats.t.ppf(1 - alfa, df), df, nc)
        tc = stats.t.ppf(1 - alfa / 2, df)
        return 1 - stats.nct.cdf(tc, df, nc) + stats.nct.cdf(-tc, df, nc)
    lo, hi = 0.0, 5.0
    for _ in range(80):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if pot(mid) < 0.8 else (lo, mid)
    return hi, se, df
for nome, lati in (("87-101", 1), ("138-147", 2)):
    c = np.array([R3["segmenti"][nome]["complesso"][f"run{r}"] for r in range(1, 6)]); l = tab[(nome, "50-500")]
    for alfa in (0.025, 0.05):
        d, se, df = rilevabile(c.std(ddof=1), 5, l.std(ddof=1), 10, alfa, lati); JS.setdefault('rilevabile', {})[f'{nome}_{alfa}'] = d
        P.append(f"  {nome} ({'unilaterale' if lati == 1 else 'bilaterale'}), alfa {alfa}: {d:.3f} A (ES {se:.4f}, gl {df:.2f}; DS {c.std(ddof=1):.3f} e {l.std(ddof=1):.3f})")
c = np.array([R3["segmenti"]["87-101"]["complesso"][f"run{r}"] for r in range(1, 6)]); l = tab[("87-101", "50-500")]
P.append(f"  scala: media del libero {l.mean():.2f}, replica meno mobile libero {l.min():.2f}, complesso {c.min():.2f}; 0.54/media = {0.536 / l.mean():.2f}")
w2 = stats.ttest_ind(c, l, equal_var=False); c2 = np.array([R3["segmenti"]["138-147"]["complesso"][f"run{r}"] for r in range(1, 6)]); l2 = tab[("138-147", "50-500")]
w3 = stats.ttest_ind(c2, l2, equal_var=False)
P.append(f"\n=== (8) 87-101 con Welch bilaterale: p {w2.pvalue:.4f}; 138-147 bilaterale p {w3.pvalue:.4f}; Holm sui due bilaterali: {min(1, 2 * min(w2.pvalue, w3.pvalue)):.3f} e {max(min(1, 2 * min(w2.pvalue, w3.pvalue)), max(w2.pvalue, w3.pvalue)):.3f}")

P.append("\n=== (3) RMSF dei Ca (A, 50-500 ns) dei residui 36-41: media sulle 15 repliche, min-max; rango del residuo 40 per replica")
allr = {k: np.array(v) for k, v in R3["rmsf"].items()}
for k in range(36, 42):
    v = np.array([x[resid == k][0] for x in allr.values()]); P.append(f"  residuo {k}: media {v.mean():.2f}, min {v.min():.2f}, max {v.max():.2f}")
rg = {n: int((x > x[resid == 40][0]).sum()) + 1 for n, x in allr.items()}
P.append("  rango del residuo 40: " + ", ".join(f"{n.replace('complesso_', 'c').replace('libero_', 'l')} {v}" for n, v in rg.items()))
P.append(f"  repliche in cui 36-39 sono i quattro piu' mobili: {sum(set(resid[np.argsort(-x)[:4]]) == {36, 37, 38, 39} for x in allr.values())}/15; "
         f"in cui 36 e 37 sono i due piu' mobili: {sum(set(resid[np.argsort(-x)[:2]]) == {36, 37} for x in allr.values())}/15; mediana di catena min-max: "
         f"{min(np.median(x) for x in allr.values()):.2f}-{max(np.median(x) for x in allr.values()):.2f}")
med = np.mean([x for x in allr.values()], axis=0); o = np.argsort(-med)[:10]
P.append("  dieci residui piu' mobili per RMSF medio sulle 15 repliche: " + ", ".join(f"{resid[i]} ({med[i]:.2f})" for i in o))

P.append("\n=== (4) seconda disposizione (phi95 > -90 e phi97 < -110): frazione dei fotogrammi per replica")
D1 = {(s, r): np.load(f"{W}/descrittivo/d1_{s}_run{r}.npz") for s, n in (("complesso", 5), ("libero", 10)) for r in range(1, n + 1)}
for (s, r), d in D1.items():
    t = d["t_ns"]; fine = 500 if s == "complesso" else 1000; w = (t >= 50 - 1e-6) & (t <= fine + 1e-6)
    res = list(d["res_died"]); a = d["phi"][w, res.index(95)]; b = d["phi"][w, res.index(97)]
    B = (a > -90) & (b < -110); A = (a < -110) & (b > -90)
    P.append(f"  {s} run{r}: seconda {B.mean():.3f}, prima (phi95 < -110 e phi97 > -90) {A.mean():.3f}" + (f"; primi 50 ns seconda {((d['phi'][t < 50, res.index(95)] > -90) & (d['phi'][t < 50, res.index(97)] < -110)).mean():.2f}" if s == "libero" else ""))

P.append("\n=== (5) correlazione di rango fra frazione dentro l'involucro a fit proprio (libero, scarto-1000 ns) e RMSF di 87-101 (scarto-500 ns)")
for sc in ("50", "100"):
    Rj = S2["secondaria"]["per_scarto"][sc]; fr = np.array([Rj["libero"][f"run{r}"]["P95"]["frazione"] for r in range(1, 11)])
    kk = [k for k in S3["per_scarto"] if float(k) == float(sc)][0]; rm = np.array([S3["per_scarto"][kk]["segmenti"]["87-101"]["libero"][f"run{r}"] for r in range(1, 11)])
    rho = stats.spearmanr(fr, rm).statistic
    g1 = [0, 1, 2, 7]; g2 = [3, 4, 5, 6, 8, 9]
    P.append(f"  scarto {sc}: rho {rho:.4f}; dentro il gruppo 1, 2, 3, 8: {stats.spearmanr(fr[g1], rm[g1]).statistic:.2f}; dentro gli altri sei: {stats.spearmanr(fr[g2], rm[g2]).statistic:.2f}; "
             f"RMSF run4 {rm[3]:.6f}, run9 {rm[8]:.6f}")

P.append("\n=== (6) legame Asp40-Ala90 (O carbossilico entro 3.5 A da N) e involucro a fit sul core (3.880 A), fotogrammi uniti")
s95 = S2["per_scarto"]["50"]["soglie"]["P95"]
for s, n, fine in (("complesso", 5, 500), ("libero", 10, 1000)):
    lg = []; pr = []
    for r in range(1, n + 1):
        d = D1[(s, r)]; t = d["t_ns"]; w = (t >= 50 - 1e-6) & (t <= fine + 1e-6)
        lg.append(np.minimum(d["D40OD1_A90N"][w], d["D40OD2_A90N"][w]) <= 3.5); pr.append(d["prim"][w])
    lg = np.concatenate(lg); pr = np.concatenate(pr)
    P.append(f"  {s}: legato {lg.mean():.3f}; dentro l'involucro: legati {(pr[lg] <= s95).mean():.3f}, non legati {(pr[~lg] <= s95).mean():.3f}; RMSD a fit sul core medio: legati {pr[lg].mean():.2f}, non legati {pr[~lg].mean():.2f}")

P.append("\n=== (7) involucri registrati (fit sul core): soglia, repliche del libero dentro (>= 0.5) e con cicli; frazione delle repliche del COMPLESSO dentro il piu' stretto")
for sc in ("50", "100"):
    Rj = S2["per_scarto"][sc]
    for nome, s in Rj["soglie"].items():
        fr = [Rj["libero"][f"run{r}"][nome]["frazione"] for r in range(1, 11)]; cy = [Rj["libero"][f"run{r}"][nome]["cicli"] for r in range(1, 11)]
        P.append(f"  scarto {sc} {nome:16s} {s:.3f} A: dentro {sum(f >= 0.5 for f in fr)}/10, cicli in {sum(c >= 1 for c in cy)}/10, esito {Rj['esiti'][nome]}")
    smin = Rj["soglie"]["P95_replica_min"]; v = []
    for r in range(1, 6):
        d = D1[("complesso", r)]; t = d["t_ns"]; w = (t >= float(sc) - 1e-6) & (t <= 500 + 1e-6); v.append(float((d["prim"][w] <= smin).mean()))
    P.append(f"  scarto {sc}: frazione dei fotogrammi di ogni replica del complesso entro {smin:.3f} A: " + " ".join(f"{x:.2f}" for x in v) + f" (>= 0.5 in {sum(x >= 0.5 for x in v)}/5)")
JS["rmsf_36_41_media"] = {str(k): float(np.mean([x[resid == k][0] for x in allr.values()])) for k in range(36, 42)}
JS["rmsf_40_minmax"] = [float(min(x[resid == 40][0] for x in allr.values())), float(max(x[resid == 40][0] for x in allr.values()))]
json.dump(JS, open(f"{W}/descrittivo/posthoc12/d11_verifiche.json", "w"), indent=1)
print("\n".join(P))
