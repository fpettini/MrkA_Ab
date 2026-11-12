"""DESCRITTIVO (post hoc, 12 nov 2026): diedri phi/psi di 86-102 nel libero, dentro e fuori dall'involucro secondario,
confrontati NELLA STESSA REPLICA, e confronto fra repliche. Corregge d3, che confrontava i fotogrammi fuori di una
replica con i fotogrammi dentro di TUTTE le repliche messe insieme (dominati dalle repliche che non escono mai) e quindi
mescolava la differenza fuori/dentro con la differenza fra repliche.
Nessun test, nessun esito.
Uso: python3 d8_diedri_per_replica.py <cartella con descrittivo/ e risultati/> <mrka-fab_completo.pdb> <cartella degli script di analisi>"""
import json, os, sys
import numpy as np
BASE, PDB, ANALISI = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, ANALISI)
import comune_fab as C
J = json.load(open(os.path.join(BASE, "risultati", "scopo2_epitopo.json")))
P95s = J["secondaria"]["per_scarto"]["50"]["soglie"]["P95"]
L = {i: np.load(os.path.join(BASE, "descrittivo", f"d1_libero_run{i}.npz")) for i in range(1, 11)}
Cx = {i: np.load(os.path.join(BASE, "descrittivo", f"d1_complesso_run{i}.npz")) for i in range(1, 6)}
res = [int(x) for x in L[1]["res_died"]]

def cmean(a):
    a = np.radians(a); return np.degrees(np.arctan2(np.sin(a).mean(0), np.cos(a).mean(0)))
def clen(a):
    a = np.radians(a); return np.hypot(np.sin(a).mean(0), np.cos(a).mean(0))
def cdiff(a, b):
    return (a - b + 180) % 360 - 180
riga = lambda v, f="{:6.0f}": "".join(f.format(x) for x in v)

print(f"involucro secondario (P95, scarto 50): {P95s:.3f} A; finestre: complesso 50-500 ns, libero 50-1000 ns")
print("medie circolari in gradi; R = lunghezza del risultante (1 = una sola popolazione stretta)")
MINF = 500   # fotogrammi minimi (50 ns) per riportare una media dentro o fuori
tab = {}
for nome in ("phi", "psi"):
    print(f"\n===== {nome} =====")
    print(f"{'':26s}" + riga(res, "{:6d}"))
    X = C.leggi_pdb(PDB)[0]["A"]
    if nome == "phi":
        cri = [float(C.diedro(X[(k - 1, "C")], X[(k, "N")], X[(k, "CA")], X[(k, "C")])) for k in res]
    else:
        cri = [float(C.diedro(X[(k, "N")], X[(k, "CA")], X[(k, "C")], X[(k + 1, "N")])) for k in res]
    print(f"{'cristallo (modello)':26s}" + riga(cri))
    for i, d in Cx.items():
        w = d["t_ns"] >= 50 - 1e-6
        print(f"{'complesso run%d' % i:26s}" + riga(cmean(d[nome][w])))
    for i, d in L.items():
        w = d["t_ns"] >= 50 - 1e-6
        den = w & (d["sec"] <= P95s); fuo = w & (d["sec"] > P95s)
        fr = den.sum() / w.sum()
        tutto = cmean(d[nome][w]); tab[(nome, i, "tutto")] = tutto
        print(f"{'libero run%d tutto (%.2f)' % (i, fr):26s}" + riga(tutto))
        print(f"{'   R':26s}" + riga(clen(d[nome][w]), "{:6.2f}"))
        if den.sum() >= MINF and fuo.sum() >= MINF:
            a, b = cmean(d[nome][den]), cmean(d[nome][fuo]); tab[(nome, i, "diff")] = cdiff(b, a)
            print(f"{'   dentro (n %d)' % den.sum():26s}" + riga(a))
            print(f"{'   fuori  (n %d)' % fuo.sum():26s}" + riga(b))
            print(f"{'   fuori - dentro':26s}" + riga(cdiff(b, a)))

print("\n===== riassunto =====")
escono = [i for i in L if ("phi", i, "diff") in tab]
print("repliche con almeno 50 ns dentro e 50 ns fuori:", escono)
for i in escono:
    for nome in ("phi", "psi"):
        dd = tab[(nome, i, "diff")]
        grandi = [(res[k], int(round(dd[k]))) for k in np.argsort(-np.abs(dd)) if abs(dd[k]) >= 30]
        print(f"run{i} {nome}: differenze fuori - dentro di almeno 30 gradi nella stessa replica: {grandi if grandi else 'nessuna'} (massimo {np.abs(dd).max():.0f})")
gA = [1, 2, 3, 8]; gB = [4, 5, 6, 7, 9, 10]
for nome in ("phi", "psi"):
    A = np.array([tab[(nome, i, "tutto")] for i in gA]); B = np.array([tab[(nome, i, "tutto")] for i in gB])
    print(f"\n{nome}: intervallo delle medie per replica (tutti i fotogrammi), repliche {gA} e repliche {gB}")
    print(f"{'':10s}" + riga(res, "{:6d}"))
    print(f"{'A min':10s}" + riga(A.min(0))); print(f"{'A max':10s}" + riga(A.max(0)))
    print(f"{'B min':10s}" + riga(B.min(0))); print(f"{'B max':10s}" + riga(B.max(0)))
    sep = [res[k] for k in range(len(res)) if A.max(0)[k] < B.min(0)[k] - 30 or A.min(0)[k] > B.max(0)[k] + 30]
    print(f"residui in cui i due gruppi sono separati di almeno 30 gradi (senza attraversare +-180): {sep}")

# ---- decorso: phi di 95 e 97 a blocchi di 100 ns (dall'inizio della produzione), per vedere se una replica cambia stato
print("\n===== decorso di phi95 / phi97 (medie circolari a blocchi di 100 ns, da 0 ns) =====")
k95, k97 = res.index(95), res.index(97)
def blocchi(d, fine):
    out = []
    for a in range(0, fine, 100):
        w = (d["t_ns"] >= a - 1e-6) & (d["t_ns"] < a + 100 - 1e-6)
        out.append(f"{cmean(d['phi'][w][:, k95]):5.0f}/{cmean(d['phi'][w][:, k97]):4.0f}")
    return " ".join(out)
for i, d in Cx.items():
    print(f"complesso run{i:<3d}" + blocchi(d, 500))
for i, d in L.items():
    print(f"libero run{i:<6d}" + blocchi(d, 1000))
print("\nprimi 10 ns di produzione (media circolare) phi95 / phi97; e frazione dei fotogrammi 50-fine con phi95 > -105 (stato delle repliche 1, 2, 3, 8)")
for nome, G in (("complesso", Cx), ("libero", L)):
    for i, d in G.items():
        w0 = d["t_ns"] < 10 - 1e-6; w = d["t_ns"] >= 50 - 1e-6
        p95 = d["phi"][w][:, k95]
        fr = np.mean((p95 > -105) & (p95 < 0))
        print(f"{nome} run{i:<3d} primi 10 ns: {cmean(d['phi'][w0][:, k95]):5.0f}/{cmean(d['phi'][w0][:, k97]):4.0f}   frazione con -105 < phi95 < 0: {fr:.2f}")
