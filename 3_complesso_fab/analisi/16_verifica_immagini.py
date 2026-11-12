"""Controllo del revisore (5 nov 2026, punto 6b): distanza dalle immagini periodiche di 10_serie.py contro
gmx mindist -pi, su alcuni fotogrammi per replica.

Uso (dalla cartella complesso_fab, con GROMACS caricato, dopo 10_serie.py):
    python3 analisi/16_verifica_immagini.py [--n 8] [--sistemi "complesso libero"] [--gmx gmx]

Per ogni replica in CACHE si prendono N fotogrammi: quello con la distanza minima dalle immagini e N-1 equidistanti.
 1. gmx convert-tpr: tpr della sola proteina (gruppo Protein_Other di <sistema>/index.ndx);
 2. gmx trjconv -pbc cluster sul gruppo Protein_Other (A+L+H nel complesso) dalla traiettoria di produzione
    <sistema>/run<r>/run<r>.xtc: il complesso esce intero e contiguo (-pbc mol non basta: rende intere le molecole
    ma puo' lasciare A e il Fab, che non sono legati, in immagini diverse);
 3. gmx mindist -pi sul risultato (tpr della sola proteina, gruppo System).
Gli stessi fotogrammi si rileggono dall'xtc di 09 usato per la cache, si riassemblano come in 10 e se ne ricalcola la
distanza dalle immagini senza il tetto di 25 A (la cache registra "> 25" come nan). Accordo atteso entro 0.005 A
(gmx scrive la distanza in nm con 3 decimali); dove la cache ha un valore, deve coincidere col ricalcolo.
Il tempo di ogni fotogramma di gmx deve essere quello della cache (stesso fotogramma).
Scrive OUT/verifica_immagini.txt e .json; i file di gmx restano in LAVORO/<sistema>_run<r>/.
"""
import os, json, shlex, argparse, subprocess, shutil, warnings
import numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
import comune_fab as C

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="analisi/cache")
ap.add_argument("--out", default="analisi/risultati")
ap.add_argument("--trj", default="analisi/traiettorie")
ap.add_argument("--lavoro", default="analisi/verifica_immagini")
ap.add_argument("--radice", default=".", help="cartella complesso_fab (per <sistema>/run<r>/run<r>.xtc e .tpr)")
ap.add_argument("--pdb", default="mrka-fab_completo.pdb")
ap.add_argument("--n", type=int, default=8, help="fotogrammi per replica (5-10)")
ap.add_argument("--sistemi", default="complesso libero")
ap.add_argument("--gmx", default=os.environ.get("GMX", "gmx"))
arg = ap.parse_args()
GMX = shlex.split(arg.gmx)
if shutil.which(GMX[0]) is None:
    raise SystemExit(f"{GMX[0]} non trovato: caricare GROMACS (o --gmx)")
os.makedirs(arg.out, exist_ok=True)
TAGLIO = 100.0                              # A: ricalcolo senza il tetto di 25 A della cache


def gmx(args, stdin, log):
    with open(log, "w") as f:
        p = subprocess.run(GMX + args, input=stdin, text=True, stdout=f, stderr=subprocess.STDOUT,
                           env={**os.environ, "GMX_MAXBACKUP": "-1"})          # niente copie #file#
    if p.returncode:
        raise SystemExit(f"gmx {args[0]} fallito (vedi {log})")


def scegli(d, n):
    nf = len(d["t_ns"]); pi = d["mindist_pi"]
    sel = set(np.round(np.linspace(0, nf - 1, max(n - 1, 2))).astype(int).tolist())
    if np.any(~np.isnan(pi)):
        sel.add(int(np.nanargmin(pi)))
    return sorted(sel)


def xtc_della_cache(d):
    P = json.loads(str(d["provenienza"]))
    x = [k for k in P["input"] if k.endswith(".xtc")]
    if len(x) != 1 or not x[0].startswith("<TRJ>/"):
        raise SystemExit(f"{d['_file']}: xtc della serie non riconosciuto nella provenienza")
    f = os.path.join(arg.trj, x[0][len("<TRJ>/"):])
    if not os.path.exists(f) or C.sha256(f) != P["input"][x[0]]:
        raise SystemExit(f"{f}: assente o diverso da quello usato per la cache (rilanciare 10)")
    return f


versione = subprocess.run(GMX + ["--version"], capture_output=True, text=True).stdout
versione = "; ".join(l.strip() for l in versione.splitlines() if l.strip().startswith(("GROMACS version", "Precision")))
serie = C.leggi_cache(arg.cache)
ris, ingressi, prod = {}, [], {}
for (s, r), d in sorted(serie.items(), key=lambda x: (x[0][0], int(x[0][1][3:]))):
    if s not in arg.sistemi.split():
        continue
    D = os.path.join(arg.radice, s); xtc = os.path.join(D, r, f"{r}.xtc"); tpr = os.path.join(D, r, f"{r}.tpr")
    ndx = os.path.join(D, "index.ndx")
    if not (os.path.exists(xtc) and os.path.exists(tpr) and os.path.exists(ndx)):
        print(f"   {s} {r}: xtc, tpr o index.ndx di produzione assenti, salto"); continue
    W = os.path.join(arg.lavoro, f"{s}_{r}"); os.makedirs(W, exist_ok=True)
    fr = scegli(d, arg.n)
    t_ps = float(d["tinit_ps"]) + d["t_ns"] * 1000.0
    ptpr = os.path.join(W, "proteina.tpr"); fndx = os.path.join(W, "fotogrammi.ndx")
    # trr, non xtc: con -pbc cluster gli atomi spostati di un vettore di cella non si riarrotondano a 0.001 nm
    ctrr = os.path.join(W, "cluster.trr"); xvg = os.path.join(W, "mindist_pi.xvg")
    gmx(["convert-tpr", "-s", tpr, "-n", ndx, "-o", ptpr], "Protein_Other\n", os.path.join(W, "convert-tpr.log"))
    open(fndx, "w").write("[ fotogrammi ]\n" + " ".join(str(f + 1) for f in fr) + "\n")     # gmx conta da 1
    fine = t_ps[fr[-1]] + 0.5 * float(d["dt_ps"])
    print(f"   {s} {r}: trjconv -pbc cluster su {len(fr)} fotogrammi (legge {os.path.basename(xtc)} fino a {fine / 1000:.1f} ns) ...", flush=True)
    gmx(["trjconv", "-f", xtc, "-s", tpr, "-n", ndx, "-fr", fndx, "-pbc", "cluster", "-e", f"{fine:.3f}", "-o", ctrr],
        "Protein_Other\nProtein_Other\n", os.path.join(W, "trjconv.log"))
    gmx(["mindist", "-f", ctrr, "-s", ptpr, "-pi", "-od", xvg], "0\n", os.path.join(W, "mindist.log"))
    righe = [[float(v) for v in l.split()] for l in open(xvg) if l.strip() and l[0] not in "#@"]
    if len(righe) != len(fr) or any(abs(x[0] - t_ps[f]) > 0.5 for x, f in zip(righe, fr)):
        raise SystemExit(f"{s} {r}: fotogrammi di gmx ({[x[0] for x in righe]} ps) diversi da quelli scelti "
                         f"({[round(float(t_ps[f]), 1) for f in fr]} ps)")
    # ricalcolo come in 10, senza tetto
    top = C.Topologia(os.path.join(arg.trj, f"proteina_{s}.gro"), s, arg.pdb)
    piano, gruppi = C.piano_riassemblaggio(top)
    pxtc = xtc_della_cache(d)
    u = mda.Universe(os.path.join(arg.trj, f"proteina_{s}.gro"), pxtc)
    righe_r = []
    for (tg, dg, *_), f in zip(righe, fr):
        ts = u.trajectory[f]
        M = ts.triclinic_dimensions[None].astype(float)
        Y = C.riassembla(ts.positions[None].astype(float), M, piano, gruppi)[0]
        py = C.mindist_immagini(Y, M[0], cut=TAGLIO)
        ca = float(d["mindist_pi"][f])
        righe_r.append({"fotogramma": int(f), "t_ns": float(d["t_ns"][f]), "gmx_A": 10.0 * dg, "python_A": float(py),
                        "diff_A": float(py - 10.0 * dg), "cache_A": None if np.isnan(ca) else ca,
                        "cache_uguale": None if np.isnan(ca) else bool(abs(ca - py) < 1e-6)})
    ris[f"{s}_{r}"] = righe_r
    ingressi += [d["_file"], pxtc, tpr, ndx]
    prod[f"{s}_{r}"] = {"xtc": os.path.relpath(xtc, arg.radice), "byte": os.path.getsize(xtc)}

if not ris:
    raise SystemExit("nessuna replica confrontata")
tutte = [x for v in ris.values() for x in v]
dmax = max(abs(x["diff_A"]) for x in tutte)
cache_ok = all(x["cache_uguale"] for x in tutte if x["cache_uguale"] is not None)
esito = "OK" if dmax <= C.PI_TOLL + 1e-4 and cache_ok else "DIFFERENZA"
P = C.prov(__file__, ingressi + [arg.pdb], {"n": arg.n, "sistemi": arg.sistemi, "gmx": arg.gmx, "cache": arg.cache,
           "trj": arg.trj, "radice": arg.radice}, {"CACHE": arg.cache, "TRJ": arg.trj, "RADICE": arg.radice})
P["gromacs"] = versione; P["xtc_produzione"] = prod
L = C.testo_prov(P) + [f"# GROMACS: {versione}", ""]
L += [f"Distanza dalle immagini periodiche: 10_serie.py (riassemblaggio per residuo, Python) contro gmx mindist -pi dopo",
      f"trjconv -pbc cluster sul gruppo Protein_Other; tolleranza {C.PI_TOLL} A (gmx scrive nm con 3 decimali).", "",
      f"{'replica':16s}{'fotog':>6s}{'t (ns)':>8s}{'gmx (A)':>9s}{'Python (A)':>11s}{'diff':>8s}{'cache (A)':>10s}"]
for x, v in ris.items():
    for y in v:
        L.append(f"{x:16s}{y['fotogramma']:6d}{y['t_ns']:8.1f}{y['gmx_A']:9.2f}{y['python_A']:11.3f}{y['diff_A']:8.3f}"
                 + (f"{y['cache_A']:10.3f}" if y["cache_A"] is not None else f"{'> 25':>10s}"))
L += ["", f"Esito: {esito}; differenza massima {dmax:.4f} A su {len(tutte)} fotogrammi di {len(ris)} repliche; "
      f"cache = ricalcolo dove la cache ha un valore: {'si' if cache_ok else 'NO'}."]
open(os.path.join(arg.out, "verifica_immagini.txt"), "w").write("\n".join(L) + "\n")
json.dump({"provenienza": P, "tolleranza_A": C.PI_TOLL, "esito": esito, "diff_max_A": dmax, "cache_uguale": cache_ok,
           "repliche": ris}, open(os.path.join(arg.out, "verifica_immagini.json"), "w"), indent=1)
print("\n".join(L))
