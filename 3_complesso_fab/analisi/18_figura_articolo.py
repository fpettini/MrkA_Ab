"""Figure per l'articolo (in inglese) dalle uscite FINALI di 12, 13, 14 e 17: SOLO DISEGNO.
Condizioni del revisore (11 nov 2026): non ricalcola e non riaggrega niente (nessuna media, nessun test, nessun esito
ricavato qui); ogni numero disegnato e' una foglia dei json letti, le cui sha256 stanno in OUT/fig_articolo_md.txt.

Uso (dalla cartella complesso_fab, dopo 11-14 e 17):
    python3 analisi/18_figura_articolo.py [--out analisi/risultati] [--pdb mrka-fab_completo.pdb]

Scrive in OUT (png e pdf):
  fig_articolo_md      A  contatti del cristallo nel complesso: frazione dei fotogrammi in contatto per replica
                          (4.0 A, scarto 50 ns) ed esito di 12 (si ferma se l'esito a scarto 100 ns e' diverso)
                       B  scopo 2: frazione dei fotogrammi del libero dentro l'involucro del legato (P95, scarto 50 ns)
                          per replica, misura primaria e secondaria, con i cicli uscita-rientro di 13
                       C  scopo 3: RMSF dei Ca per residuo, una curva per replica (scarto 50 ns), tratti 87-101 e 138-147
                          (asse tagliato a YMAX_C: la coda N-terminale esce di scala e la figura lo scrive)
                       D  scopo 3: medie di segmento per replica e p di Welch corretto con Holm, come scritti da 14
  fig_articolo_md_S1   scopo 2 nel tempo: frazione dentro l'involucro per blocchi di 50 ns ("blocchi" di 13), due misure
  fig_articolo_md_S2   rotazione dell'Fv rispetto al core di MrkA nel tempo, per replica del complesso (serie di 17);
                       saltata se diagnostica_cerniera.json manca o non dichiara lo stato FINALE
  fig_articolo_md.txt  provenienza: sha256 dello script e dei json letti, e che cosa mostra ogni pannello
Si ferma se uno dei tre scopi non e' FINALE o se manca la misura secondaria del libero.
Colori: complesso blu, libero arancio (come 15_figure.py).
"""
import os, json, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import comune_fab as C

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="analisi/risultati")
ap.add_argument("--pdb", default="mrka-fab_completo.pdb")
arg = ap.parse_args()

BLU, ARANCIO = "#2a78d6", "#eb6834"
INK, INK2, GRIGIO, FASCIA, RETE = "#0b0b0b", "#52514e", "#b5b3ad", "#f1f0ec", "#e6e5e0"
RAMPA_BLU = LinearSegmentedColormap.from_list("blu", ["#ffffff", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
RAMPA_ARANCIO = LinearSegmentedColormap.from_list("arancio", ["#ffffff", "#fde0d2", "#f8b294", "#ee8053", "#cf5521", "#8a3510"])
plt.rcParams.update({"font.size": 7, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "lines.linewidth": 0.9, "figure.dpi": 150,
                     "savefig.dpi": 300, "pdf.fonttype": 42})

F = {k: os.path.join(arg.out, f"{k}.json") for k in ("scopo1_contatti", "scopo2_epitopo", "scopo3_rmsf", "diagnostica_cerniera")}
s1, s2, s3 = (json.load(open(F[k])) for k in ("scopo1_contatti", "scopo2_epitopo", "scopo3_rmsf"))
for nome, s in (("scopo 1", s1), ("scopo 2", s2), ("scopo 3", s3)):
    if s["stato"] != "FINALE":
        raise SystemExit(f"{nome} non e' FINALE: le figure dell'articolo si fanno solo sui dati finali")
if not s2.get("finale_per_lunghezze") or not s2["secondaria"]["libero_calcolato"]:
    raise SystemExit("scopo 2: misura secondaria del libero assente (rilanciare 13 con gli script attuali)")
s17 = json.load(open(F["diagnostica_cerniera"])) if os.path.exists(F["diagnostica_cerniera"]) else None
if s17 is not None and s17.get("stato") != "FINALE":
    print("diagnostica_cerniera.json senza stato FINALE (17 vecchio o dati non finali): figura S2 saltata"); s17 = None
letti = [F[k] for k in ("scopo1_contatti", "scopo2_epitopo", "scopo3_rmsf")] + ([F["diagnostica_cerniera"]] if s17 else [])
_, rn = C.leggi_pdb(arg.pdb)
ESITO = {"persistente": "persistent", "perso": "lost", "senza esito": "no outcome", "nuovo": "new"}
YMAX_C = 4.6                              # A: limite dell'asse del pannello C (la coda N-terminale esce di scala, dichiarato)
SC = str(int(C.SCARTI[0])); SC2 = str(int(C.SCARTI[1])); CUT = f"{C.CONTATTO_PRINCIPALE:.1f}"
rep = lambda d: sorted(d, key=lambda r: int(r[3:]))


def salva(fig, nome):
    for est in ("png", "pdf"):
        fig.savefig(os.path.join(arg.out, f"{nome}.{est}"))
    plt.close(fig); print("scritta", os.path.join(arg.out, nome + ".png"))


# =================================================================== figura principale
fig = plt.figure(figsize=(7.2, 8.1))
gs = fig.add_gridspec(2, 2, height_ratios=[1.6, 1.0], width_ratios=[1.0, 1.3], hspace=0.50, wspace=0.50,
                      left=0.175, right=0.985, top=0.925, bottom=0.065)

# ---------------- A: contatti del cristallo nel complesso
axA = fig.add_subplot(gs[0, 0])
r50 = s1["risultati"][CUT][SC]; r100 = s1["risultati"][CUT][SC2]
crist = [(k, v) for k, v in r50["coppie"].items() if v["cristallo"]]
crist.sort(key=lambda kv: (int(kv[0].split("|")[0]), kv[0].split("|")[1], int(kv[0].split("|")[2])))
Z = np.array([v["frazioni"] for _, v in crist])
im = axA.imshow(Z, cmap=RAMPA_BLU, vmin=0, vmax=1, aspect="auto")
nomi, esiti = [], []
for k, v in crist:
    a, c, f = k.split("|")
    if r100["coppie"][k]["esito"] != v["esito"]:
        raise SystemExit(f"{k}: esito diverso con scarto {SC} e {SC2} ns ('non robusto'): la figura mostra solo esiti uguali, "
                         "va deciso come disegnarla")
    nomi.append(f"{rn['A'][int(a)].capitalize()}{a} – {c} {rn[c][int(f)].capitalize()}{f}"); esiti.append(ESITO[v["esito"]])
axA.set_xticks(range(Z.shape[1]), [r[3:] for r in r50["repliche"]])
axA.set_yticks(range(len(nomi)), nomi, fontsize=5.8)
axA.set_xticks(np.arange(-0.5, Z.shape[1]), minor=True); axA.set_yticks(np.arange(-0.5, len(nomi)), minor=True)
axA.grid(which="minor", color="white", lw=1.6); axA.tick_params(which="both", length=0)
for sp in axA.spines.values():
    sp.set_visible(False)
axA.set_xlabel("complex replica", labelpad=2)
for i, e in enumerate(esiti):
    axA.text(Z.shape[1] - 0.38, i, e, va="center", ha="left", fontsize=5.6, color=INK2 if e == "persistent" else INK,
             fontweight="normal" if e == "persistent" else "bold", clip_on=False)
pa = axA.get_position()
cax = fig.add_axes([pa.x0, pa.y0 - 0.058, pa.width, 0.008])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label(f"fraction of frames in contact ({CUT} Å)", fontsize=6.5, labelpad=2)
cb.outline.set_linewidth(0.4); cb.ax.tick_params(length=2, labelsize=6)
fig.text(0.012, pa.y1 + 0.014, "A", fontweight="bold", fontsize=10, va="bottom")

# ---------------- B: scopo 2, frazione dentro l'involucro per replica del libero
axB = fig.add_subplot(gs[0, 1])
p50 = s2["per_scarto"][SC]; q50 = s2["secondaria"]["per_scarto"][SC]
reps = rep(p50["libero"]); x = np.arange(1, len(reps) + 1)
fp = [p50["libero"][r]["P95"]["frazione"] for r in reps]; fs = [q50["libero"][r]["P95"]["frazione"] for r in reps]
cp = [p50["libero"][r]["P95"]["cicli"] for r in reps]; cs = [q50["libero"][r]["P95"]["cicli"] for r in reps]
axB.axhline(0.5, color=INK2, ls="--", lw=0.6); axB.axhline(0.05, color=INK2, ls=":", lw=0.6)
DX_B = 0.16                               # i due simboli di una replica stanno affiancati: sovrapposti, quello vuoto copriva il pieno
for xi, a, b in zip(x, fp, fs):
    axB.plot([xi - DX_B, xi + DX_B], [a, b], color=GRIGIO, lw=0.8, zorder=1)
axB.plot(x - DX_B, fp, "o", ms=5.5, color=ARANCIO, mec="white", mew=0.8, zorder=3,
         label=f"core-fitted (registered), envelope {p50['soglie']['P95']:.2f} Å")
axB.plot(x + DX_B, fs, "o", ms=4.8, color="white", mec=ARANCIO, mew=1.3, zorder=4,
         label=f"self-fitted (post hoc), envelope {q50['soglie']['P95']:.2f} Å")
axB.set_xlim(0.4, len(reps) + 0.6); axB.set_ylim(-0.02, 1.03)
axB.set_xticks(x, [r[3:] for r in reps]); axB.set_ylabel("fraction of frames inside the bound envelope")
axB.text(len(reps) + 0.58, 0.51, "0.5", ha="right", va="bottom", fontsize=6, color=INK2)
axB.text(len(reps) + 0.58, 0.06, "0.05", ha="right", va="bottom", fontsize=6, color=INK2)
axB.legend(frameon=False, loc="lower left", bbox_to_anchor=(-0.02, 1.005), fontsize=6.3, handletextpad=0.3, borderaxespad=0)
tr = axB.get_xaxis_transform()
axB.text(0.3, -0.036, "free replica", transform=tr, ha="right", va="center", fontsize=6.5, color=INK, clip_on=False)
axB.text(0.3, -0.085, "cycles, core-fitted", transform=tr, ha="right", va="center", fontsize=6, color=INK2, clip_on=False)
axB.text(0.3, -0.125, "cycles, self-fitted", transform=tr, ha="right", va="center", fontsize=6, color=INK2, clip_on=False)
for xi, a, b in zip(x, cp, cs):
    axB.text(xi, -0.085, str(a), transform=tr, ha="center", va="center", fontsize=6, color=INK if a else GRIGIO, clip_on=False)
    axB.text(xi, -0.125, str(b), transform=tr, ha="center", va="center", fontsize=6, color=INK if b else GRIGIO, clip_on=False)
axB.tick_params(axis="x", pad=1.5, length=2)
pb = axB.get_position()
fig.text(pb.x0 - 0.115, pa.y1 + 0.014, "B", fontweight="bold", fontsize=10, va="bottom")

# ---------------- C, D: scopo 3 (una curva e un punto per replica: nessuna media calcolata qui)
sub = gs[1, :].subgridspec(1, 2, width_ratios=[4, 1.2], wspace=0.24)
axC = fig.add_subplot(sub[0, 0]); axD = fig.add_subplot(sub[0, 1])
r = s3["per_scarto"][SC]; resid = np.array(r["resid"])
for nome, v in C.SEGMENTI.items():
    axC.axvspan(v[0] - 0.5, v[-1] + 0.5, color=FASCIA, lw=0, zorder=0)
    axC.text((v[0] + v[-1]) / 2, 1.0, nome.replace("-", "–"), transform=axC.get_xaxis_transform(), ha="center", va="bottom",
             fontsize=6.3, color=INK2)
for s, col, et in (("libero", ARANCIO, "free MrkA"), ("complesso", BLU, "complex")):
    chiavi = rep([k[len(s) + 1:] for k in r["rmsf"] if k.startswith(s + "_")])
    for i, k in enumerate(chiavi):
        axC.plot(resid, r["rmsf"][f"{s}_{k}"], color=col, lw=0.55, alpha=0.6,
                 label=f"{et} ({len(chiavi)} replicas, one line each)" if i == 0 else None)
axC.plot(C.EPITOPO, np.zeros(len(C.EPITOPO)), "|", color=INK, ms=6, mew=0.9, transform=axC.get_xaxis_transform(), clip_on=False)
axC.set_xlabel("MrkA residue (ticks on the axis: epitope)"); axC.set_ylabel("Cα RMSF (Å)")
h, l = axC.get_legend_handles_labels(); leg = axC.legend(h[::-1], l[::-1], frameon=False, loc="upper right", fontsize=6.5)
for lh in leg.legend_handles:
    lh.set_alpha(1.0); lh.set_linewidth(1.2)
fuori = sorted({int(resid[i]) for v in r["rmsf"].values() for i, y in enumerate(v) if y > YMAX_C})
if fuori and (fuori[0] != int(resid.min()) or fuori != list(range(fuori[0], fuori[-1] + 1))):
    raise SystemExit(f"pannello C: RMSF oltre {YMAX_C} A fuori dalla coda N-terminale ({fuori}): cambiare YMAX_C")
axC.set_xlim(resid.min() - 1, resid.max() + 1); axC.set_ylim(0, YMAX_C)
if fuori:
    axC.text(fuori[-1] + 1.5, YMAX_C * 0.985, f"← {fuori[0]}–{fuori[-1]} off scale", ha="left", va="top", fontsize=6, color=INK2)
axC.yaxis.grid(True, color=RETE, lw=0.4); axC.set_axisbelow(True)
pc = axC.get_position()
fig.text(0.012, pc.y1 + 0.014, "C", fontweight="bold", fontsize=10, va="bottom")
for i, nome in enumerate(C.SEGMENTI):
    for j, (s, col) in enumerate((("complesso", BLU), ("libero", ARANCIO))):
        v = [r["segmenti"][nome][s][k] for k in rep(r["segmenti"][nome][s])]
        axD.plot(i * 3 + j + np.linspace(-0.2, 0.2, len(v)), v, "o", ms=3.8, color=col, mec="white", mew=0.7)
    axD.text(i * 3 + 0.5, 1.0, f"p = {r['test'][nome]['welch_p_holm']:.2f}", transform=axD.get_xaxis_transform(), ha="center",
             va="bottom", fontsize=6.3, color=INK2)
axD.set_xticks([0.5, 3.5], [n.replace("-", "–") for n in C.SEGMENTI]); axD.set_ylabel("segment mean RMSF per replica (Å)")
axD.set_xlim(-0.8, 4.8); axD.set_ylim(bottom=0); axD.yaxis.grid(True, color=RETE, lw=0.4); axD.set_axisbelow(True)
pd_ = axD.get_position()
fig.text(pd_.x0 - 0.075, pc.y1 + 0.014, "D", fontweight="bold", fontsize=10, va="bottom")
salva(fig, "fig_articolo_md")

# =================================================================== S1: scopo 2 nel tempo, blocchi di 50 ns
fig, ax = plt.subplots(2, 1, figsize=(7.2, 4.9), sharex=True)
fig.subplots_adjust(left=0.085, right=0.79, top=0.93, bottom=0.115, hspace=0.34)
for a, d, tit in ((ax[0], p50, "core-fitted (registered)"), (ax[1], q50, "self-fitted (post hoc)")):
    B = [d["libero"][k]["P95"]["blocchi"] for k in reps]
    t0 = [b[0] for b in B[0]]; t1 = [b[1] for b in B[0]]
    if any([b[0] for b in bb] != t0 or [b[1] for b in bb] != t1 for bb in B):
        raise SystemExit("scopo 2: i blocchi non sono gli stessi per tutte le repliche del libero")
    Z = np.array([[b[2] for b in bb] for bb in B])
    im = a.imshow(Z, cmap=RAMPA_ARANCIO, vmin=0, vmax=1, aspect="auto", extent=(t0[0], t1[-1], len(reps) + 0.5, 0.5))
    a.set_yticks(range(1, len(reps) + 1), [k[3:] for k in reps]); a.set_ylabel("free replica")
    a.set_xticks(np.array(t0[1:]), minor=True); a.set_yticks(np.arange(1.5, len(reps)), minor=True)
    a.grid(which="minor", color="white", lw=1.0); a.tick_params(which="both", length=0)
    for sp in a.spines.values():
        sp.set_visible(False)
    a.set_title(f"{tit}, envelope {d['soglie']['P95']:.2f} Å", loc="left", fontsize=7.5, color=INK, pad=4)
    tr = a.get_yaxis_transform()
    a.text(1.035, 0.2, "fraction", transform=tr, ha="center", va="center", fontsize=6, color=INK2, clip_on=False)
    a.text(1.115, 0.2, "cycles", transform=tr, ha="center", va="center", fontsize=6, color=INK2, clip_on=False)
    for i, k in enumerate(reps):
        a.text(1.035, i + 1, f"{d['libero'][k]['P95']['frazione']:.2f}", transform=tr, ha="center", va="center", fontsize=6,
               color=INK, clip_on=False)
        a.text(1.115, i + 1, str(d["libero"][k]["P95"]["cicli"]), transform=tr, ha="center", va="center", fontsize=6,
               color=INK, clip_on=False)
ax[1].set_xlabel(f"production time (ns); {int(t1[0] - t0[0])}-ns blocks from {int(t0[0])} ns")
cax = fig.add_axes([0.925, 0.30, 0.012, 0.42])
cb = fig.colorbar(im, cax=cax); cb.set_label("fraction of frames inside the bound envelope", fontsize=6.5)
cb.outline.set_linewidth(0.4); cb.ax.tick_params(length=2, labelsize=6)
salva(fig, "fig_articolo_md_S1")

# =================================================================== S2: rotazione dell'Fv (serie di 17)
if s17:
    rr = rep(s17["repliche"]); sc17 = s17["provenienza"]["argomenti"]["scarto"]
    fig, ax = plt.subplots(1, len(rr), figsize=(7.2, 2.1), sharex=True, sharey=True, squeeze=False)
    fig.subplots_adjust(left=0.075, right=0.99, top=0.86, bottom=0.2, wspace=0.1)
    for a, k in zip(ax[0], rr):
        a.axvspan(0, sc17, color=FASCIA, lw=0, zorder=0)
        a.plot(s17["repliche"][k]["t_ns"], s17["repliche"][k]["serie"]["rot_fv"], color=BLU, lw=0.6)
        a.set_title(f"complex replica {k[3:]}", loc="left", fontsize=7, color=INK, pad=3)
        a.yaxis.grid(True, color=RETE, lw=0.4); a.set_axisbelow(True); a.set_xlabel("production time (ns)")
    ax[0, 0].set_ylabel("Fv rotation relative to\nthe MrkA core (degrees)"); ax[0, 0].set_ylim(bottom=0)
    salva(fig, "fig_articolo_md_S2")

P = C.prov(__file__, letti + [arg.pdb], {"out": arg.out, "pdb": arg.pdb}, {"OUT": arg.out})
L = C.testo_prov(P)
L += ["# SOLO DISEGNO: nessun ricalcolo e nessuna aggregazione; ogni numero disegnato e' una foglia dei json qui sopra.", "",
      "fig_articolo_md",
      f"  A  scopo1_contatti.json: risultati[{CUT}][{SC}].coppie con cristallo = vero: frazioni per replica ed esito",
      f"     (controllato uguale all'esito con scarto {SC2} ns); nomi dei residui dal pdb.",
      f"  B  scopo2_epitopo.json: per_scarto[{SC}] (primaria) e secondaria.per_scarto[{SC}]: libero.<replica>.P95.frazione e",
      "     .cicli, soglie.P95; linee a 0.5 e 0.05 = soglie dei criteri.",
      f"  C  scopo3_rmsf.json: per_scarto[{SC}].rmsf, una curva per replica; tratti ed epitopo da comune_fab; asse tagliato a",
      f"     {YMAX_C} A, fuori scala i residui " + (f"{fuori[0]}-{fuori[-1]}" if fuori else "nessuno") + ".",
      f"  D  scopo3_rmsf.json: per_scarto[{SC}].segmenti (un punto per replica) e test.<tratto>.welch_p_holm.",
      "fig_articolo_md_S1",
      f"  scopo2_epitopo.json: libero.<replica>.P95.blocchi (frazione dentro l'involucro per blocco), .frazione e .cicli,",
      f"  per_scarto[{SC}] (primaria, sopra) e secondaria.per_scarto[{SC}] (sotto).",
      "fig_articolo_md_S2" + ("" if s17 else " (NON scritta: diagnostica_cerniera.json assente o senza stato FINALE)"),
      "  diagnostica_cerniera.json: repliche.<replica>.t_ns e .serie.rot_fv; fascia = scarto usato da 17."]
open(os.path.join(arg.out, "fig_articolo_md.txt"), "w").write("\n".join(L) + "\n")
print("scritta", os.path.join(arg.out, "fig_articolo_md.txt"))
