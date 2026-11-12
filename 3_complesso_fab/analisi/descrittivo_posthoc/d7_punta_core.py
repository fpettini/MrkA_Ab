#!/usr/bin/env python3
"""Descrittivo POST HOC (12 nov 2026): scostamento dal cristallo dei Ca dell'epitopo dopo fit sul core, per residuo, nel
complesso (50-500 ns) e nel libero (50-1000 ns). Stessa definizione di 17_cerniera.py (per_residuo_finestra: media sui
fotogrammi della distanza fra il Ca dopo fit sul core, dalla cache di 10, e il Ca del cristallo), estesa al libero.
In piu': RMSF dei Ca dei residui 36-40 per replica (da scopo3_rmsf.json, scarto 50, 50-500 ns).
Non e' un criterio, nessun test. Uso: python3 d7_punta_core.py <cache> <risultati> <pdb> <cartella degli script di analisi>"""
import sys, json, numpy as np
cache, ris, pdb, scr = sys.argv[1:5]
sys.path.insert(0, scr); sys.dont_write_bytecode = True
import comune_fab as C
xyz, rn = C.leggi_pdb(pdb)
EPI = list(C.EPITOPO)
rif = np.array([xyz["A"][(k, "CA")] for k in EPI], float)
print("POST HOC 12 nov 2026, descrittivo, nessun test. Scostamento medio (A) dei Ca dal cristallo dopo fit sul core.")
print("replica          finestra   " + "".join(f"{k:>7d}" for k in EPI))
out = {}
for sis, n, t1 in (("complesso", 5, 500.0), ("libero", 10, 1000.0)):
    for r in range(1, n + 1):
        d = np.load(f"{cache}/serie_{sis}_run{r}.npz", allow_pickle=True)
        t = d["t_ns"]; w = (t >= 50.0) & (t <= t1)
        idx = [list(d["ca_resid"]).index(k) for k in EPI]
        dev = np.linalg.norm(d["ca_fit"][:, idx].astype(float) - rif, axis=-1)[w].mean(0)
        out[(sis, r)] = dev
        print(f"{sis+' run'+str(r):<16} 50-{t1:<6.0f} " + "".join(f"{v:7.2f}" for v in dev))
for sis, n in (("complesso", 5), ("libero", 10)):
    M = np.array([out[(sis, r)] for r in range(1, n + 1)])
    punta = [EPI.index(k) for k in (90, 91, 92, 93)]; altri = [i for i in range(len(EPI)) if i not in punta]
    print(f"{sis}: Ca 90-93 min {M[:, punta].min():.2f} max {M[:, punta].max():.2f}; altri residui dell'epitopo min {M[:, altri].min():.2f} max {M[:, altri].max():.2f}; "
          f"residuo 91 per replica " + ", ".join(f"{v:.1f}" for v in M[:, EPI.index(91)]))
s3 = json.load(open(f"{ris}/scopo3_rmsf.json"))
P = s3["per_scarto"]["50"]; resid = P["resid"] if "resid" in P else s3["resid"]
print()
print("RMSF dei Ca (A), scarto 50, 50-500 ns, residui 36-40 e, per confronto, mediana di tutti i residui")
for k, v in P["rmsf"].items():
    v = np.array(v); i = [resid.index(x) for x in (36, 37, 38, 39, 40)]
    print(f"{k:<16} " + " ".join(f"{v[j]:6.2f}" for j in i) + f"   mediana {np.median(v):.2f}")
tutti = np.array([np.array(v)[[resid.index(x) for x in (36, 37, 38, 39, 40)]] for v in P["rmsf"].values()])
print("residuo 40: min %.2f max %.2f su 15 repliche; residui 36-40: min %.2f max %.2f" % (tutti[:, 4].min(), tutti[:, 4].max(), tutti.min(), tutti.max()))
