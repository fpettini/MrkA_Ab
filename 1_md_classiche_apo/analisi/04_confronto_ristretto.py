"""Fase 4: confronto ristretto fra i sistemi sulle grandezze convergenti (core beta, S-S, Rg).

Uso (dalla cartella del progetto mrkA):
    python analisi_convergenza/04_confronto_ristretto.py . analisi_convergenza/cache analisi_convergenza/500ns
Legge le cache .npz di 01_estrai.py e 02_analisi.py. L'unita' statistica e' la replica:
ogni grandezza e' mediata su 50-500 ns dentro la replica, poi media e SD fra le 3 repliche.
Scrive confronto_ristretto.txt e confronto_ristretto.json (righe per i grafici del report).
"""
import sys, os, json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
CACHE = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "analisi_convergenza", "cache")
OUT = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "analisi_convergenza", "500ns")
SISTEMI = ["xray", "nmr_ox", "nmr_red"]
REPLICHE = [1, 2, 3]
RES = np.arange(42, 203)
DT = 0.1                    # ns fra due frame
T0 = 50.0                   # ns scartati come equilibrazione
PASSO_DSSP = 10             # DSSP ogni 1 ns
STRAND = [44, 45, 46, 47, 78, 79, 80, 81, 82, 83, 84, 85, 88, 94, 102, 103, 104, 105, 107, 108, 122, 134,
          135, 136, 137, 138, 148, 162, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 195, 196, 197, 198, 199]
STRAND2 = [k for k in STRAND if k >= 57]          # core beta senza lo strand 44-47 (sta sul solco del donatore)
SEGMENTI = [(44, 47), (78, 85), (102, 108), (134, 138), (170, 179), (195, 199)]
BINS = np.arange(1.5, 9.01, 0.25)                 # istogramma SG-SG (A)

def carica(s, r):
    d = np.load(os.path.join(CACHE, f"{s}_run{r}.npz"), allow_pickle=True)
    d = {k: d[k] for k in d.files}
    d["names"] = d["names"].astype(str)
    return d

def atomo(d, nome, resid):
    return np.where((d["names"] == nome) & (d["resids"] == resid))[0][0]

def fit(P, Q, sel):
    Ps = P[:, sel]; Pm = Ps.mean(1, keepdims=True); Qm = Q[sel].mean(0)
    H = np.einsum('fni,nj->fij', Ps - Pm, Q[sel] - Qm)
    U, S, Vt = np.linalg.svd(H)
    D = np.tile(np.eye(3), (len(P), 1, 1)); D[:, 2, 2] = np.sign(np.linalg.det(U @ Vt))
    return np.einsum('fni,fij->fnj', P - Pm, U @ D @ Vt) + Qm

def rmsd(P, Q, sel):
    Pa = fit(P, Q, sel)
    return np.sqrt(((Pa[:, sel] - Q[sel]) ** 2).sum(2).mean(1))

def ca_core(d, lista):
    return np.array([atomo(d, "CA", k) for k in lista])

# strutture di partenza: core beta di ciascun sistema, per l'RMSD incrociato
REF = {}
for s in SISTEMI:
    d = carica(s, 1)
    REF[s] = d["ref"][ca_core(d, STRAND)].astype(np.float64)
i2 = np.array([STRAND.index(k) for k in STRAND2])

def una_replica(arg):
    s, r = arg
    d = carica(s, r)
    serie = np.load(os.path.join(CACHE, f"{s}_run{r}_serie.npz"))
    w = serie["t"] >= T0
    P = d["X"][w][:, ca_core(d, STRAND)].astype(np.float64)
    tutti = np.arange(len(STRAND))
    out = {"rb": rmsd(P, REF[s], tutti).mean(),
           "rb2": rmsd(P, REF[s], i2).mean(),
           "rx": rmsd(P, REF["xray"], tutti).mean(),        # core verso la struttura X-ray di partenza
           "rn": rmsd(P, REF["nmr_ox"], tutti).mean(),      # core verso la struttura NMR di partenza
           "rg": serie["rg"][w].mean(),
           "ss": serie["ss"][w].mean(), "ssmax": serie["ss"][w].max(),
           "ss35": 100 * (serie["ss"][w] < 3.5).mean(),
           "hist": np.histogram(serie["ss"][w], BINS)[0] / w.sum()}
    u = mda.Universe(os.path.join(ROOT, s, "conf.gro"))
    u.load_new(d["X"][w][::PASSO_DSSP].copy(), format=mda.coordinates.memory.MemoryReader)
    e = DSSP(u.select_atoms("protein")).run().results.dssp == "E"
    out["beta_tot"] = e[:, 44 - 42:199 - 41].mean()
    out["beta_seg"] = [e[:, a - 42:b - 41].mean() for a, b in SEGMENTI]
    return s, r, out

if __name__ == "__main__":
    with Pool(4) as pool:
        ris = pool.map(una_replica, [(s, r) for s in SISTEMI for r in REPLICHE])
    R = {s: [o for ss, r, o in ris if ss == s] for s in SISTEMI}
    righe = []
    def stampa(t=""):
        print(t); righe.append(t)
    stampa(f"Confronto ristretto, finestra {T0:.0f}-500 ns, media fra 3 repliche +/- SD (valori delle repliche)")
    xs = ca_core(carica("xray", 1), STRAND)
    stampa(f"Core beta fra le strutture di partenza: X-ray vs NMR {rmsd(REF['nmr_ox'][None], REF['xray'], np.arange(len(STRAND)))[0]:.2f} A "
           f"({len(STRAND)} Ca); NMR ox vs NMR red {rmsd(REF['nmr_red'][None], REF['nmr_ox'], np.arange(len(STRAND)))[0]:.2f} A")
    voci = [("rb", "RMSD core beta dalla propria partenza (A)", "{:.2f}"),
            ("rb2", "RMSD core beta senza 44-47 dalla propria partenza (A)", "{:.2f}"),
            ("rx", "RMSD core beta dalla struttura X-ray (A)", "{:.2f}"),
            ("rn", "RMSD core beta dalla struttura NMR (A)", "{:.2f}"),
            ("beta_tot", "frazione beta DSSP 44-199", "{:.3f}"),
            ("rg", "Rg atomi pesanti (A)", "{:.2f}"),
            ("ss", "SG43-SG87 media (A)", "{:.3f}"),
            ("ssmax", "SG43-SG87 max (A)", "{:.2f}"),
            ("ss35", "frame con SG-SG < 3.5 A (%)", "{:.1f}")]
    J = {"sistemi": {}, "segmenti": [f"{a}-{b}" for a, b in SEGMENTI], "bins": BINS.tolist()}
    for k, nome, f in voci:
        stampa(f"\n{nome}")
        for s in SISTEMI:
            v = np.array([o[k] for o in R[s]])
            stampa(f"   {s:8s} {f.format(v.mean())} +/- {f.format(v.std(ddof=1))}   ({', '.join(f.format(x) for x in v)})")
    stampa("\nFrazione beta DSSP per segmento: media fra repliche (min-max)")
    for s in SISTEMI:
        B = np.array([o["beta_seg"] for o in R[s]])
        stampa(f"   {s:8s} " + " | ".join(f"{a}-{b} {B[:, j].mean():.2f} ({B[:, j].min():.2f}-{B[:, j].max():.2f})"
                                          for j, (a, b) in enumerate(SEGMENTI)))
    stampa("\nDifferenze fra sistemi (media delle repliche; 'separati' = intervalli min-max delle repliche non sovrapposti)")
    for a, b, cosa in (("xray", "nmr_ox", "struttura + 8 posizioni di sequenza"), ("nmr_ox", "nmr_red", "solo stato redox")):
        stampa(f"   {a} vs {b} ({cosa})")
        for k, nome, f in voci[:6]:
            va = np.array([o[k] for o in R[a]]); vb = np.array([o[k] for o in R[b]])
            sep = va.max() < vb.min() or vb.max() < va.min()
            stampa(f"      {nome:55s} {f.format(vb.mean() - va.mean()):>7s}   {'separati' if sep else 'sovrapposti'}")
    for s in SISTEMI:
        J["sistemi"][s] = [{k: (np.asarray(o[k]).tolist()) for k in o} for o in R[s]]
    os.makedirs(OUT, exist_ok=True)
    open(os.path.join(OUT, "confronto_ristretto.txt"), "w").write("\n".join(righe) + "\n")
    json.dump(J, open(os.path.join(OUT, "confronto_ristretto.json"), "w"), indent=1)
