"""DESCRITTIVO (post hoc): dentro/fuori dall'involucro secondario nel libero, residuo per residuo."""
import json, numpy as np
np.set_printoptions(linewidth=220, suppress=True, precision=2)
J = json.load(open("risultati/scopo2_epitopo.json")); S3 = json.load(open("risultati/scopo3_rmsf.json"))
P95s = J["secondaria"]["per_scarto"]["50"]["soglie"]["P95"]; P95p = J["per_scarto"]["50"]["soglie"]["P95"]
EPI = [42, 86, 88, 90, 91, 92, 93, 94, 96, 167, 168]
L = {f"run{i}": np.load(f"descrittivo/d1_libero_run{i}.npz") for i in range(1, 11)}
Cx = {f"run{i}": np.load(f"descrittivo/d1_complesso_run{i}.npz") for i in range(1, 6)}
rm = S3["per_scarto"]["50"]["segmenti"]["87-101"]["libero"]

def tratti(dentro, t):
    """tratto continuo piu' lungo fuori (ns), numero di fotogrammi fuori"""
    best = 0; cur = 0
    for v in dentro:
        cur = 0 if v else cur + 1
        best = max(best, cur)
    return best * (t[1] - t[0])

print("== libero 50-1000 ns ==")
print(f"{'rep':6s}{'fr.sec':>7s}{'fr.prim':>8s}{'sec med':>8s}{'prim med':>9s}{'fuori max ns':>13s}{'RMSF87-101':>11s}{'corr(sec,prim)':>15s}")
fs, rr = [], []
for r, d in L.items():
    w = d["t_ns"] >= 50 - 1e-6; t = d["t_ns"][w]
    ds = d["sec"][w] <= P95s; dp = d["prim"][w] <= P95p
    fs.append(ds.mean()); rr.append(rm[r])
    print(f"{r:6s}{ds.mean():7.2f}{dp.mean():8.2f}{d['sec'][w].mean():8.2f}{d['prim'][w].mean():9.2f}{tratti(ds, t):13.1f}{rm[r]:11.2f}{np.corrcoef(d['sec'][w], d['prim'][w])[0,1]:15.2f}")
from scipy.stats import spearmanr
print(f"Spearman fra frazione dentro (sec, 50-1000) e RMSF 87-101 (50-500), 10 repliche: rho {spearmanr(fs, rr)[0]:.2f} (descrittivo, nessun test)")

print("\n== scostamento per residuo (backbone, A) dopo fit sull'epitopo: complesso 50-500 (media su repliche) e libero dentro / fuori ==")
print("residuo      ", "".join(f"{k:>7d}" for k in EPI))
mc = np.mean([d["dev_self"][d["t_ns"] >= 50 - 1e-6].mean(0) for d in Cx.values()], axis=0)
print("complesso    ", "".join(f"{x:7.2f}" for x in mc))
for r, d in L.items():
    w = d["t_ns"] >= 50 - 1e-6; ds = d["sec"][w] <= P95s; dv = d["dev_self"][w]
    a = dv[ds].mean(0) if ds.any() else np.full(11, np.nan); b = dv[~ds].mean(0) if (~ds).sum() > 50 else np.full(11, np.nan)
    print(f"{r:6s} dentro", "".join(f"{x:7.2f}" for x in a))
    if (~ds).sum() > 50:
        print(f"{r:6s} fuori ", "".join(f"{x:7.2f}" for x in b), f"  (n fuori {(~ds).sum()})")

print("\n== scostamento per residuo dopo fit sul CORE (backbone, A): complesso per replica e libero per replica, 50-fine ==")
print("residuo      ", "".join(f"{k:>7d}" for k in EPI))
for r, d in Cx.items():
    print(f"C {r:5s}     ", "".join(f"{x:7.2f}" for x in d["dev_core"][d["t_ns"] >= 50 - 1e-6].mean(0)))
for r, d in L.items():
    print(f"L {r:5s}     ", "".join(f"{x:7.2f}" for x in d["dev_core"][d["t_ns"] >= 50 - 1e-6].mean(0)))

print("\n== diedri phi/psi 86-102: differenza circolare media fuori - dentro (gradi) nelle repliche con molto 'fuori' ==")
def cmean(a):
    a = np.radians(a); return np.degrees(np.arctan2(np.sin(a).mean(0), np.cos(a).mean(0)))
def cdiff(a, b):
    return (a - b + 180) % 360 - 180
res = L["run1"]["res_died"]
# riferimento "dentro": tutte le repliche del libero, fotogrammi dentro
phi_in = np.concatenate([d["phi"][(d["t_ns"] >= 50 - 1e-6) & (d["sec"] <= P95s)] for d in L.values()])
psi_in = np.concatenate([d["psi"][(d["t_ns"] >= 50 - 1e-6) & (d["sec"] <= P95s)] for d in L.values()])
print("residuo        ", "".join(f"{k:>6d}" for k in res))
print("phi dentro     ", "".join(f"{x:6.0f}" for x in cmean(phi_in)))
print("psi dentro     ", "".join(f"{x:6.0f}" for x in cmean(psi_in)))
phc = np.concatenate([d["phi"][d["t_ns"] >= 50 - 1e-6] for d in Cx.values()]); psc = np.concatenate([d["psi"][d["t_ns"] >= 50 - 1e-6] for d in Cx.values()])
print("phi complesso  ", "".join(f"{x:6.0f}" for x in cmean(phc)))
print("psi complesso  ", "".join(f"{x:6.0f}" for x in cmean(psc)))
for r in ("run1", "run2", "run3", "run7", "run8"):
    d = L[r]; w = (d["t_ns"] >= 50 - 1e-6) & (d["sec"] > P95s)
    print(f"{r} fuori dphi", "".join(f"{x:6.0f}" for x in cdiff(cmean(d['phi'][w]), cmean(phi_in))))
    print(f"{r} fuori dpsi", "".join(f"{x:6.0f}" for x in cdiff(cmean(d['psi'][w]), cmean(psi_in))))

print("\n== distanze (A): media complesso 50-500; libero dentro / fuori ==")
nomi = ["D40OD1_A90N", "D40OD2_A90N", "C43SG_C87SG", "CA91_CA168", "CA92_CA42", "CA91_CA96", "CA90_CA94", "CA92_CA88"]
print(f"{'':14s}" + "".join(f"{n:>13s}" for n in nomi))
print(f"{'complesso':14s}" + "".join(f"{np.mean([d[n][d['t_ns'] >= 50 - 1e-6].mean() for d in Cx.values()]):13.2f}" for n in nomi))
for r, d in L.items():
    w = d["t_ns"] >= 50 - 1e-6; ds = d["sec"][w] <= P95s
    print(f"{r+' dentro':14s}" + "".join(f"{d[n][w][ds].mean():13.2f}" for n in nomi))
    if (~ds).sum() > 50:
        print(f"{r+' fuori':14s}" + "".join(f"{d[n][w][~ds].mean():13.2f}" for n in nomi))
