"""Diagnostica descrittiva della rotazione Fv-core nel complesso (6 nov 2026, dopo l'analisi intermedia a 250 ns).
NON e' un criterio del doc: non decide validita' ne' esiti e non cambia nulla degli scopi 1-3.
Solo complesso: nel libero l'RMSD dell'epitopo a fit proprio non si calcola prima dei dati finali (revisore, 6 nov).

Uso (dalla cartella complesso_fab, dopo 09 e 10):
    python3 analisi/17_cerniera.py [--cache analisi/cache] [--trj analisi/traiettorie] [--out analisi/risultati]
                                   [--pdb mrka-fab_completo.pdb] [--passo 10]

Per ogni replica del complesso in CACHE si legge un fotogramma ogni PASSO dall'xtc di 09 usato per la serie (sha256
controllato contro la provenienza della cache), si riassembla come in 10 e si calcola:
  epi|core   RMSD del backbone dell'epitopo dopo fit sui 43 Ca del core (= rmsd_epi della cache: controllo < 1e-4 A)
  epi|Fv     lo stesso dopo fit sui Ca dell'Fv (VL 1-107 + VH 1-117, ricostruiti esclusi)
  epi|epi    lo stesso dopo fit sull'epitopo stesso (forma dell'epitopo)
  tip|core, tip|Fv   RMSD dei Ca 90-93 dopo fit sul core e dopo fit sull'Fv
  rot Fv     angolo della rotazione che porta l'Fv, nel riferimento del core, sull'Fv del cristallo (gradi)
  trasl Fv   spostamento del centro dei Ca dell'Fv nel riferimento del core (A)
Le stesse misure su proteina_complesso.gro (fine di equil_s4, partenza comune delle repliche), e per residuo lo
scostamento dal cristallo dei Ca dell'epitopo dopo fit sul core (dalla cache, ca_fit), al primo ns e da SCARTO in poi.
Scrive OUT/diagnostica_cerniera.txt, OUT/diagnostica_cerniera.json e OUT/fig_diagnostica_cerniera.png (con matplotlib).
"""
import os, json, argparse, warnings
import numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
import comune_fab as C

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="analisi/cache")
ap.add_argument("--trj", default="analisi/traiettorie")
ap.add_argument("--out", default="analisi/risultati")
ap.add_argument("--pdb", default="mrka-fab_completo.pdb")
ap.add_argument("--passo", type=int, default=10, help="un fotogramma ogni PASSO (10 = 1 ns con fotogrammi da 100 ps)")
ap.add_argument("--scarto", type=float, default=C.SCARTI[0])
arg = ap.parse_args()
os.makedirs(arg.out, exist_ok=True)
S = "complesso"
PUNTA = [90, 91, 92, 93]
TOLL_CACHE = 1e-4


def xtc_della_cache(d):
    P = json.loads(str(d["provenienza"]))
    x = [k for k in P["input"] if k.endswith(".xtc")]
    if len(x) != 1 or not x[0].startswith("<TRJ>/"):
        raise SystemExit(f"{d['_file']}: xtc della serie non riconosciuto nella provenienza")
    f = os.path.join(arg.trj, x[0][len("<TRJ>/"):])
    if not os.path.exists(f) or C.sha256(f) != P["input"][x[0]]:
        raise SystemExit(f"{f}: assente o diverso da quello usato per la cache (rilanciare 10)")
    return f


def angolo(R):
    return np.degrees(np.arccos(np.clip((np.trace(R, axis1=1, axis2=2) - 1.0) / 2.0, -1.0, 1.0)))


gro = os.path.join(arg.trj, f"proteina_{S}.gro")
xyz, rn = C.leggi_pdb(arg.pdb)
rif = lambda c, k, n: xyz[c][(k, n)]
top = C.Topologia(gro, S, arg.pdb)
piano, gruppi = C.piano_riassemblaggio(top)
i_core = np.array([top.atomo("A", k, "CA") for k in C.CORE43]); r_core = np.array([rif("A", k, "CA") for k in C.CORE43])
epi = [(k, n) for k in C.EPITOPO for n in C.BACKBONE]
i_epi = np.array([top.atomo("A", k, n) for k, n in epi]); r_epi = np.array([rif("A", k, n) for k, n in epi])
fv = [(c, k) for c in ("L", "H") for k, _, _ in top.res[c] if C.FV[c][0] <= k <= C.FV[c][1] and k not in C.RICOSTRUITI[c]]
i_fv = np.array([top.atomo(c, k, "CA") for c, k in fv]); r_fv = np.array([rif(c, k, "CA") for c, k in fv])
i_tip = np.array([top.atomo("A", k, "CA") for k in PUNTA]); r_tip = np.array([rif("A", k, "CA") for k in PUNTA])
r_epi_ca = np.array([rif("A", k, "CA") for k in C.EPITOPO])


def misure(Y):
    """Y (F, n, 3) gia' riassemblato -> dizionario di serie (F,)"""
    R, cP, cQ = C.kabsch(Y[:, i_core], r_core)
    Fc = C.applica(Y[:, i_fv], R, cP, cQ)                      # Fv nel riferimento del core
    R2, _, _ = C.kabsch(Fc, r_fv)
    R3, c3, q3 = C.kabsch(Y[:, i_fv], r_fv)
    R4, c4, q4 = C.kabsch(Y[:, i_epi], r_epi)
    return {"epi_core": C.rmsd(C.applica(Y[:, i_epi], R, cP, cQ), r_epi),
            "epi_fv": C.rmsd(C.applica(Y[:, i_epi], R3, c3, q3), r_epi),
            "epi_epi": C.rmsd(C.applica(Y[:, i_epi], R4, c4, q4), r_epi),
            "tip_core": C.rmsd(C.applica(Y[:, i_tip], R, cP, cQ), r_tip),
            "tip_fv": C.rmsd(C.applica(Y[:, i_tip], R3, c3, q3), r_tip),
            "rot_fv": angolo(R2),
            "trasl_fv": np.linalg.norm(Fc.mean(1) - r_fv.mean(0), axis=-1)}


NOMI = ["epi_core", "epi_fv", "epi_epi", "tip_core", "tip_fv", "rot_fv", "trasl_fv"]
TITOLI = ["epi|core", "epi|Fv", "epi|epi", "tip|core", "tip|Fv", "rot Fv", "trasl Fv"]

# partenza comune: fine di equil_s4
u0 = mda.Universe(gro)
Y0 = C.riassembla(u0.atoms.positions.copy()[None].astype(float), u0.trajectory.ts.triclinic_dimensions.copy()[None].astype(float),
                  piano, gruppi)
s4 = {k: float(v[0]) for k, v in misure(Y0).items()}

serie = C.leggi_cache(arg.cache)
ris, ingressi, max_diff = {}, [], 0.0
for (s, r), d in sorted(serie.items(), key=lambda x: (x[0][0], int(x[0][1][3:]))):
    if s != S:
        continue                                            # il libero non si tocca (condizione del revisore)
    xtc = xtc_della_cache(d)
    u = mda.Universe(gro, xtc)
    if u.trajectory.n_frames != len(d["t_ns"]):
        raise SystemExit(f"{s} {r}: fotogrammi nell'xtc {u.trajectory.n_frames}, nella cache {len(d['t_ns'])}")
    fr = np.arange(0, u.trajectory.n_frames, arg.passo)
    X = np.empty((len(fr), top.n, 3)); M = np.empty((len(fr), 3, 3)); T = np.empty(len(fr))
    for j, f in enumerate(fr):
        ts = u.trajectory[int(f)]
        X[j] = ts.positions; M[j] = ts.triclinic_dimensions; T[j] = ts.time
    t_ns = d["t_ns"][fr]
    if np.abs((T - float(d["tinit_ps"])) / 1000.0 - t_ns).max() > 1e-3:
        raise SystemExit(f"{s} {r}: tempi dell'xtc diversi da quelli della cache")
    m = misure(C.riassembla(X, M, piano, gruppi))
    diff = float(np.abs(m["epi_core"] - d["rmsd_epi"][fr]).max()); max_diff = max(max_diff, diff)
    # per residuo, dalla cache: Ca dell'epitopo dopo fit sul core
    idx = [list(d["ca_resid"]).index(k) for k in C.EPITOPO]
    dev = np.linalg.norm(d["ca_fit"][:, idx].astype(float) - r_epi_ca, axis=-1)
    w0 = d["t_ns"] < 1.0 - 1e-6; w = C.finestra(d["t_ns"], arg.scarto, float(d["t_ns"][-1]))
    wf = C.finestra(t_ns, arg.scarto, float(d["t_ns"][-1]))
    ris[r] = {"t_ns": t_ns.tolist(), "serie": {k: v.tolist() for k, v in m.items()}, "diff_cache_A": diff,
              "primo": {k: float(v[0]) for k, v in m.items()},
              "finestra": {k: {"media": float(v[wf].mean()), "min": float(v[wf].min()), "max": float(v[wf].max())} for k, v in m.items()},
              "corr_epi_core_rot": float(np.corrcoef(m["epi_core"][wf], m["rot_fv"][wf])[0, 1]),
              "per_residuo_primo_ns": dict(zip(map(str, C.EPITOPO), dev[w0].mean(0).tolist())),
              "per_residuo_finestra": dict(zip(map(str, C.EPITOPO), dev[w].mean(0).tolist()))}
    ingressi += [d["_file"], xtc]
    print(f"   {s} {r}: {len(fr)} fotogrammi, scarto con la cache {diff:.2e} A", flush=True)

if not ris:
    raise SystemExit(f"nessuna replica del {S} in {arg.cache}")
if max_diff > TOLL_CACHE:
    raise SystemExit(f"epi|core diverso da rmsd_epi della cache ({max_diff:.2e} A > {TOLL_CACHE}): riassemblaggio o fit non coincidono con 10")

P = C.prov(__file__, ingressi + [gro, arg.pdb], {"cache": arg.cache, "trj": arg.trj, "out": arg.out, "pdb": arg.pdb,
           "passo": arg.passo, "scarto": arg.scarto}, {"CACHE": arg.cache, "TRJ": arg.trj, "OUT": arg.out})
L = C.testo_prov(P)
L += ["# DIAGNOSTICA DESCRITTIVA, non e' un criterio del doc: non decide validita' ne' esiti. Solo complesso (nel libero",
      "# l'RMSD dell'epitopo a fit proprio non si calcola prima dei dati finali, condizione del revisore del 6 nov 2026).",
      *C.avviso_intermedia("INTERMEDIA"), "",
      f"Un fotogramma ogni {arg.passo} della serie; epi|core coincide con rmsd_epi della cache (scarto massimo {max_diff:.1e} A).",
      "epi = backbone dei residui " + ", ".join(map(str, C.EPITOPO)) + "; tip = Ca " + "-".join(map(str, (PUNTA[0], PUNTA[-1])))
      + f"; Fv = {len(fv)} Ca (ricostruiti esclusi); core = 43 Ca. RMSD e traslazione in A, rotazione in gradi.", "",
      f"{'':22s}" + "".join(f"{x:>10s}" for x in TITOLI),
      f"{'fine equil_s4 (gro)':22s}" + "".join(f"{s4[k]:10.2f}" for k in NOMI)]
for r, v in ris.items():
    L.append(f"{S + ' ' + r + ' primo':22s}" + "".join(f"{v['primo'][k]:10.2f}" for k in NOMI))
    L.append(f"{S + ' ' + r + ' ' + str(int(arg.scarto)) + '-fine':22s}" + "".join(f"{v['finestra'][k]['media']:10.2f}" for k in NOMI))
L += ["", f"Da {int(arg.scarto)} ns in poi: minimo-massimo e correlazione fra epi|core e rot Fv"]
for r, v in ris.items():
    f = v["finestra"]
    L.append(f"{S} {r}: rot Fv {f['rot_fv']['min']:.0f}-{f['rot_fv']['max']:.0f} gradi, trasl Fv {f['trasl_fv']['min']:.1f}-"
             f"{f['trasl_fv']['max']:.1f} A, tip|Fv max {f['tip_fv']['max']:.2f} A, epi|epi max {f['epi_epi']['max']:.2f} A, "
             f"corr(epi|core, rot Fv) {v['corr_epi_core_rot']:.2f}")
L += ["", "Ca dell'epitopo: scostamento dal cristallo dopo fit sul core (A, dalla cache; primo = primo ns di produzione)",
      f"{'':22s}" + "".join(f"{rn['A'][k] + str(k):>8s}" for k in C.EPITOPO)]
for r, v in ris.items():
    for lab, key in (("primo", "per_residuo_primo_ns"), (f"{int(arg.scarto)}-fine", "per_residuo_finestra")):
        L.append(f"{S + ' ' + r + ' ' + lab:22s}" + "".join(f"{v[key][str(k)]:8.1f}" for k in C.EPITOPO))
open(os.path.join(arg.out, "diagnostica_cerniera.txt"), "w").write("\n".join(L) + "\n")
json.dump({"provenienza": P, "equil_s4": s4, "repliche": ris, "diff_cache_max_A": max_diff},
          open(os.path.join(arg.out, "diagnostica_cerniera.json"), "w"), indent=1)
print("\n".join(L))

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib non installato: figura saltata"); raise SystemExit(0)
COL = {"core": "#2a78d6", "fv": "#eb6834", "epi": "#52514e"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e5e0", "grid.linewidth": 0.5, "lines.linewidth": 0.8, "figure.dpi": 150})
n = len(ris)
fig, ax = plt.subplots(2, n, figsize=(3.2 * n + 0.6, 5.0), sharex=True, sharey="row", squeeze=False)
for j, (r, v) in enumerate(ris.items()):
    t = np.array(v["t_ns"]); sv = {k: np.array(x) for k, x in v["serie"].items()}
    for a in ax[:, j]:
        a.axvspan(0, arg.scarto, color="#efeee9", zorder=0)
    ax[0, j].plot(t, sv["rot_fv"], color=COL["core"]); ax[0, j].set_title(f"{S} {r}", loc="left", fontsize=9)
    ax[1, j].plot(t, sv["epi_core"], color=COL["core"], label="fit sul core MrkA (misura dello scopo 2)")
    ax[1, j].plot(t, sv["epi_fv"], color=COL["fv"], label="fit sull'Fv")
    ax[1, j].plot(t, sv["epi_epi"], color=COL["epi"], label="fit sull'epitopo stesso")
    ax[1, j].set_xlabel("tempo di produzione (ns)")
ax[0, 0].set_ylabel("rotazione dell'Fv rispetto\nal core MrkA (gradi)"); ax[1, 0].set_ylabel("RMSD backbone epitopo (A)")
h, l = ax[1, 0].get_legend_handles_labels(); fig.legend(h, l, frameon=False, loc="lower center", ncol=3)
fig.suptitle(f"Diagnostica descrittiva, non criterio (un fotogramma ogni {arg.passo}; grigio = scarto {int(arg.scarto)} ns) - "
             "ANALISI INTERMEDIA", x=0.01, ha="left", fontsize=9)
fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(os.path.join(arg.out, "fig_diagnostica_cerniera.png"), dpi=150)
print("figura: fig_diagnostica_cerniera.png")
