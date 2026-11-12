"""DESCRITTIVO (post hoc): riassunto di che cosa cambia quando il libero esce dall'involucro secondario."""
import json, numpy as np
np.set_printoptions(linewidth=200, suppress=True)
J = json.load(open("risultati/scopo2_epitopo.json"))
S3 = json.load(open("risultati/scopo3_rmsf.json"))
P95s = J["secondaria"]["per_scarto"]["50"]["soglie"]["P95"]; P95p = J["per_scarto"]["50"]["soglie"]["P95"]
print(f"soglie scarto 50: primaria {P95p:.3f}, secondaria {P95s:.3f}")
EPI = [42, 86, 88, 90, 91, 92, 93, 94, 96, 167, 168]
D = {}
for s, n in (("complesso", 5), ("libero", 10)):
    for i in range(1, n + 1):
        D[(s, f"run{i}")] = np.load(f"descrittivo/d1_{s}_run{i}.npz")
# controllo: prim = cache rmsd_epi; sec libero = medie di 13
mx = 0
for (s, r), d in D.items():
    c = np.load(f"cache/serie_{s}_{r}.npz", allow_pickle=True)
    mx = max(mx, np.abs(c["rmsd_epi"] - d["prim"]).max())
    if s == "complesso": mx = max(mx, np.abs(c["rmsd_epi_proprio"] - d["sec"]).max())
print(f"controllo: |prim - cache rmsd_epi| e |sec - cache| (complesso) max {mx:.1e}")
for r in [f"run{i}" for i in range(1, 11)]:
    d = D[("libero", r)]; w = d["t_ns"] >= 50 - 1e-6
    m13 = J["secondaria"]["per_scarto"]["50"]["libero"][r]["rmsd_media"]
    assert abs(d["sec"][w].mean() - m13) < 1e-9, (r, d["sec"][w].mean(), m13)
print("controllo: medie della secondaria nel libero = uscita di 13 (1e-9)")

print("\n== libero, finestra 50-1000 ns: frazione dentro (sec P95), tratto fuori piu' lungo, RMSF 87-101 (50-500, scopo 3) ==")
rmsf3 = S3["per_scarto"]["50"]
def cerca(o, chiave, prof=0):
    if isinstance(o, dict):
        for k, v in o.items():
            if k == chiave: return v
            x = cerca(v, chiave, prof + 1)
            if x is not None: return x
    return None
print("chiavi scopo3 per_scarto/50:", list(rmsf3.keys())[:12])
seg = rmsf3.get("segmenti") or rmsf3
print(json.dumps(seg, indent=0)[:600])
