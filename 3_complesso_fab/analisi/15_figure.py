"""Figure descrittive delle analisi del complesso MrkA-Fab e del libero (dopo 11-14).

Uso (dalla cartella complesso_fab):
    python3 analisi/15_figure.py [--cache analisi/cache] [--out analisi/risultati]

Scrive in OUT/:
  fig_controlli.png      RMSD dei 43 Ca del core, contatti del cristallo presenti (Q), crescita della distanza dei centri,
                         distanza minima dalle immagini, per replica, con le soglie di validita'
  fig_scopo1_contatti.png  frazione di fotogrammi in contatto per ogni contatto del cristallo (4.0 A, scarto 50 ns)
  fig_scopo2_epitopo.png   RMSD del backbone dell'epitopo per replica, con la soglia P95 dal complesso
  fig_scopo3_rmsf.png      RMSF dei Ca per residuo (media e intervallo fra repliche) e medie di segmento per replica
Le figure sono solo descrittive; gli esiti stanno nei txt. Colori: complesso blu, libero arancio.
"""
import os, json, argparse
import numpy as np
import comune_fab as C

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit("matplotlib non installato: le figure si saltano (i risultati sono nei txt)")

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="analisi/cache")
ap.add_argument("--out", default="analisi/risultati")
arg = ap.parse_args()

COL = {"complesso": "#2a78d6", "libero": "#eb6834"}
INK, INK2, GRIGIO = "#0b0b0b", "#52514e", "#b5b3ad"
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e5e0",
                     "grid.linewidth": 0.5, "lines.linewidth": 0.8, "figure.dpi": 150, "savefig.bbox": "tight"})

ctrl = C.leggi_controlli(arg.out)
serie = C.leggi_cache(arg.cache)
T3, T2 = ctrl["T3_ns"], ctrl["T2_ns"]
s1 = json.load(open(os.path.join(arg.out, "scopo1_contatti.json")))
s2 = json.load(open(os.path.join(arg.out, "scopo2_epitopo.json")))
s3 = json.load(open(os.path.join(arg.out, "scopo3_rmsf.json")))
interm = any(v != "FINALE" for v in ctrl["stato"].values())
sott = " - ANALISI INTERMEDIA, descrittiva" if interm else ""
reps = {s: C.ordina_repliche([r for (x, r) in serie if x == s]) for s in ("complesso", "libero")}


def media_mobile(y, n=10):
    """media mobile su n fotogrammi (1 ns), normalizzata ai bordi"""
    if len(y) < n:
        return y
    k = np.ones(n)
    return np.convolve(y, k, mode="same") / np.convolve(np.ones(len(y)), k, mode="same")


def scarto(ax):
    ax.axvspan(0, C.SCARTI[0], color="#f1f0ec", zorder=0, lw=0)


# --- controlli
fig, axs = plt.subplots(2, 2, figsize=(9, 5.5), sharex=True)
pannelli = [("rmsd_core", "RMSD 43 Ca del core (A)", C.RMSD_CORE_MAX), ("Q", "contatti del cristallo presenti (frazione)", 0.5),
            ("dcentri", "crescita distanza dei centri (A)", C.CENTRI_MAX), ("mindist_pi", "distanza minima dalle immagini (A)", C.PI_MIN)]
cr = [tuple(p) for p in ctrl["contatti_cristallo_4A"]]
for ax, (k, tit, soglia) in zip(axs.ravel(), pannelli):
    for s in ("complesso", "libero"):
        for r in reps[s]:
            d = serie[(s, r)]; t = d["t_ns"]
            if k == "Q":
                if s != "complesso":
                    continue
                sel = d["c_d"] <= C.CONTATTO_PRINCIPALE
                cc = np.where(d["c_c"][sel] == 0, "L", "H")
                pres = {}
                for f, a, c, kk in zip(d["c_f"][sel], d["c_a"][sel], cc, d["c_k"][sel]):
                    if (int(a), c, int(kk)) in cr:
                        pres[f] = pres.get(f, 0) + 1
                y = np.array([pres.get(f, 0) for f in range(len(t))]) / len(cr)
            elif k == "dcentri":
                if s != "complesso":
                    continue
                y = d["dcentri"] - float(d["dcentri_cristallo"])
            elif k == "mindist_pi":
                y = np.where(np.isnan(d["mindist_pi"]), C.PI_CUT, d["mindist_pi"])
            else:
                y = d[k]
            ax.plot(t, media_mobile(y), color=COL[s], alpha=0.75, label=s if r == reps[s][0] else None)
    ax.axhline(soglia, color=INK2, ls="--", lw=0.7)
    scarto(ax)
    ax.set_title(tit, loc="left", color=INK, fontsize=8)
    if k == "mindist_pi":
        ax.set_ylim(0, C.PI_CUT + 1)
        ax.text(ax.get_xlim()[1], C.PI_CUT, "> 25 A disegnato a 25", ha="right", va="bottom", color=INK2, fontsize=6)
for ax in axs[1]:
    ax.set_xlabel("tempo di produzione (ns)")
axs[0, 0].legend(frameon=False, loc="upper left")
fig.suptitle(f"Controlli per replica (media mobile 1 ns; tratteggio = soglia; grigio = scarto 50 ns){sott}", x=0.01, ha="left", fontsize=9)
fig.savefig(os.path.join(arg.out, "fig_controlli.png")); plt.close(fig)

# --- scopo 1
r50 = s1["risultati"][str(C.CONTATTO_PRINCIPALE)]["50"]
crist = [(k, v) for k, v in r50["coppie"].items() if v["cristallo"]]
crist.sort(key=lambda kv: (int(kv[0].split("|")[0]), kv[0].split("|")[1], int(kv[0].split("|")[2])))
if crist and r50["repliche"]:
    _, rn = C.leggi_pdb("mrka-fab_completo.pdb") if os.path.exists("mrka-fab_completo.pdb") else (None, None)
    nomi = []
    for k, _ in crist:
        a, c, f = k.split("|")
        nomi.append(f"A{a}{rn['A'][int(a)] if rn else ''} - {c}{f}{rn[c][int(f)] if rn else ''}")
    Z = np.array([v["frazioni"] for _, v in crist])
    fig, ax = plt.subplots(figsize=(1.6 + 0.55 * Z.shape[1], 0.6 + 0.2 * len(nomi)))
    im = ax.imshow(Z, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(Z.shape[1]), r50["repliche"]); ax.set_yticks(range(len(nomi)), nomi, fontsize=6.5)
    ax.grid(False)
    ax.set_xticks(np.arange(-0.5, Z.shape[1]), minor=True); ax.set_yticks(np.arange(-0.5, len(nomi)), minor=True)
    ax.grid(which="minor", color="white", lw=2); ax.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.05); cb.set_label("frazione di fotogrammi in contatto (4.0 A)")
    ax.set_title(f"Contatti del cristallo nel complesso, scarto 50 ns, fino a {T3:.0f} ns{sott}", loc="left", fontsize=8)
    fig.savefig(os.path.join(arg.out, "fig_scopo1_contatti.png")); plt.close(fig)

# --- scopo 2
p50 = s2["per_scarto"]["50"]["soglie"]["P95"] if "50" in s2["per_scarto"] else None
tutte = [("complesso", r) for r in reps["complesso"]] + [("libero", r) for r in reps["libero"]]
nc = 5
nr = int(np.ceil(len(tutte) / nc))
fig, axs = plt.subplots(nr, nc, figsize=(11, 1.7 * nr), sharex=True, sharey=True, squeeze=False)
for ax, (s, r) in zip(axs.ravel(), tutte):
    d = serie[(s, r)]
    ax.plot(d["t_ns"], d["rmsd_epi"], color=COL[s], lw=0.4, alpha=0.9)
    if p50 is not None:
        ax.axhline(p50, color=INK, ls="--", lw=0.7)
    scarto(ax)
    ax.axvline(T3 if s == "complesso" else T2, color=GRIGIO, lw=0.6)
    ax.set_title(f"{s} {r}", loc="left", fontsize=7, color=COL[s])
for ax in axs.ravel()[len(tutte):]:
    ax.axis("off")
for ax in axs[-1]:
    ax.set_xlabel("ns")
for ax in axs[:, 0]:
    ax.set_ylabel("RMSD epitopo (A)")
fig.suptitle(f"Scopo 2: RMSD del backbone dell'epitopo dal cristallo (fit 43 Ca); tratteggio = soglia P95 dal complesso "
             f"({p50:.2f} A, scarto 50 ns){sott}" if p50 is not None else "Scopo 2", x=0.01, ha="left", fontsize=9)
fig.savefig(os.path.join(arg.out, "fig_scopo2_epitopo.png")); plt.close(fig)

# --- scopo 3
r = s3["per_scarto"]["50"]
resid = np.array(r["resid"])
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 3.2), gridspec_kw={"width_ratios": [4, 1.3]})
for nome, res in C.SEGMENTI.items():
    ax.axvspan(min(res) - 0.5, max(res) + 0.5, color="#f1f0ec", lw=0, zorder=0)
    ax.text((min(res) + max(res)) / 2, 1.0, nome, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7, color=INK2)
for s in ("complesso", "libero"):
    Y = np.array([v for k, v in r["rmsf"].items() if k.startswith(s + "_")])
    if len(Y) == 0:
        continue
    ax.fill_between(resid, Y.min(0), Y.max(0), color=COL[s], alpha=0.18, lw=0)
    ax.plot(resid, Y.mean(0), color=COL[s], lw=1.2, label=f"{s} (n = {len(Y)}, media e intervallo)")
ax.plot(C.EPITOPO, np.zeros(len(C.EPITOPO)), "|", color=INK, ms=8, transform=ax.get_xaxis_transform(), clip_on=False)
ax.set_xlabel("residuo (numerazione wt; | = epitopo)"); ax.set_ylabel("RMSF Ca (A)")
ax.legend(frameon=False, loc="upper right"); ax.set_xlim(resid.min() - 1, resid.max() + 1)
for i, nome in enumerate(C.SEGMENTI):
    for j, s in enumerate(("complesso", "libero")):
        v = list(r["segmenti"][nome][s].values())
        x = i * 3 + j + np.linspace(-0.15, 0.15, len(v)) if v else []
        ax2.plot(x, v, "o", ms=4, color=COL[s], mec="white", mew=0.8)
        if v:
            ax2.plot([i * 3 + j - 0.3, i * 3 + j + 0.3], [np.mean(v)] * 2, color=INK, lw=1)
ax2.set_xticks([0.5, 3.5], list(C.SEGMENTI)); ax2.set_ylabel("media di segmento per replica (A)")
ax2.set_xlim(-0.8, 4.8)
fig.suptitle(f"Scopo 3: RMSF dei Ca, scarto 50 ns, fino a {T3:.0f} ns per entrambi i sistemi{sott}", x=0.01, ha="left", fontsize=9)
fig.savefig(os.path.join(arg.out, "fig_scopo3_rmsf.png")); plt.close(fig)
print("figure scritte in", arg.out)
