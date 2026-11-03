"""Fase 2: verifica di convergenza sulle traiettorie estratte da 01_estrai.py.

Uso (dalla cartella del progetto mrkA):
    python analisi_convergenza/02_analisi.py . analisi_convergenza/cache
Scrive analisi_convergenza/cache/*.npz (usati da 03_figure.py) e stampa tutti i numeri
riportati in RIEPILOGO_convergenza.md. Richiede MDAnalysis >= 2.8 e pymbar >= 4.
"""
import sys, os, warnings, numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP
from pymbar import timeseries

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
CACHE = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "analisi_convergenza", "cache")
SISTEMI = ["xray", "nmr_ox", "nmr_red"]
REPLICHE = [1, 2, 3]
RES = np.arange(42, 203)
CORE = (44, 199)           # esclusi 42-43 e 200-202 (cariche artificiali ai terminali)
DT = 0.1                   # ns fra due frame

def carica(s, r):
    d = np.load(os.path.join(CACHE, f"{s}_run{r}.npz"), allow_pickle=True)
    d = {k: d[k] for k in d.files}
    d["names"] = d["names"].astype(str)
    return d

def atomo(d, nome, resid):
    return np.where((d["names"] == nome) & (d["resids"] == resid))[0][0]

def fit(P, Q, sel):
    """Sovrappone ogni frame di P (F,N,3) su Q (N,3) usando gli atomi sel (Kabsch)."""
    Ps = P[:, sel]; Pm = Ps.mean(1, keepdims=True); Qm = Q[sel].mean(0)
    H = np.einsum('fni,nj->fij', Ps - Pm, Q[sel] - Qm)
    U, S, Vt = np.linalg.svd(H)
    D = np.tile(np.eye(3), (len(P), 1, 1)); D[:, 2, 2] = np.sign(np.linalg.det(U @ Vt))
    return np.einsum('fni,fij->fnj', P - Pm, U @ D @ Vt) + Qm

def rmsd(P, Q, sel):
    Pa = fit(P, Q, sel)
    return np.sqrt(((Pa[:, sel] - Q[sel]) ** 2).sum(2).mean(1))

def rmsip(A, B, k=10):           # Amadei et al. 1999
    return np.sqrt(((A[:, :k].T @ B[:, :k]) ** 2).sum() / k)

def cosine_content(p, i):        # Hess 2002
    N = len(p); tt = (np.arange(N) + 0.5) / N
    return 2.0 / N * (np.cos(i * np.pi * tt) @ p) ** 2 / (p @ p)

NBLOCCHI = 10

def blocchi(y, n=NBLOCCHI):
    return [y[i * len(y) // n:(i + 1) * len(y) // n].mean() for i in range(n)]

# ---------------------------------------------------------------- core beta comune
print("1. Core beta: residui in strand (DSSP) nella struttura di partenza di tutti e tre i sistemi")
E = {}
for s in SISTEMI:
    u = mda.Universe(os.path.join(ROOT, s, "conf.gro"))
    code = "".join(DSSP(u.select_atoms("protein")).run().results.dssp[0])
    E[s] = {42 + i for i, c in enumerate(code) if c == "E"}
    print(f"   {s:8s} {code}")
STRAND = sorted(k for k in E["xray"] & E["nmr_ox"] & E["nmr_red"] if CORE[0] <= k <= CORE[1])
print(f"   {len(STRAND)} Ca di strand comuni dentro {CORE[0]}-{CORE[1]}: {STRAND}\n")

# ---------------------------------------------------------------- serie per frame
print("2. Controlli di integrita' e serie temporali (riferimento = conf.gro, struttura di partenza)")
serie = {}
for s in SISTEMI:
    for r in REPLICHE:
        d = carica(s, r); X = d["X"].astype(np.float64); ref = d["ref"].astype(np.float64); t = d["t"] / 1000
        pep = np.linalg.norm(X[:, [atomo(d, "C", k) for k in range(42, 202)]]
                             - X[:, [atomo(d, "N", k) for k in range(43, 203)]], axis=2)
        ss = np.linalg.norm(X[:, atomo(d, "SG", 43)] - X[:, atomo(d, "SG", 87)], axis=1)
        ca_core = np.array([atomo(d, "CA", k) for k in range(CORE[0], CORE[1] + 1)])
        ca_beta = np.array([atomo(d, "CA", k) for k in STRAND])
        heavy = np.array([not n.startswith("H") for n in d["names"]])
        Xh = X[:, heavy]
        rg = np.sqrt(((Xh - Xh.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
        serie[s, r] = dict(t=t, ss=ss, rg=rg,
                           rmsd=rmsd(X[:, ca_core], ref[ca_core], np.arange(len(ca_core))),
                           rb=rmsd(X[:, ca_beta], ref[ca_beta], np.arange(len(ca_beta))))
        np.savez(os.path.join(CACHE, f"{s}_run{r}_serie.npz"), **serie[s, r])
        print(f"   {s:8s} r{r}: {len(t)} frame, dt {np.diff(t).min():.3f}-{np.diff(t).max():.3f} ns, "
              f"legame peptidico max {pep.max():.2f} A (proteina intera se < 2)")

print("\n3. Distanza SG43-SG87 (A)")
for s in SISTEMI:
    for r in REPLICHE:
        ss = serie[s, r]["ss"]
        print(f"   {s:8s} r{r}: media {ss.mean():.3f}  min {ss.min():.2f}  95perc {np.percentile(ss, 95):.2f}  "
              f"max {ss.max():.2f}  frame < 3.5 A: {np.mean(ss < 3.5) * 100:.1f}%")

for chiave, nome in (("rmsd", "RMSD Ca 44-199"), ("rb", "RMSD core beta"), ("rg", "Rg atomi pesanti")):
    tmax = serie[SISTEMI[0], 1]["t"][-1]; lb = tmax / NBLOCCHI
    print(f"\n4. {nome} (A), medie su blocchi da {lb:.0f} ns (inizio blocco in ns)")
    print("   " + " " * 11 + "".join(f"{a * lb:>6.0f}" for a in range(NBLOCCHI)))
    for s in SISTEMI:
        for r in REPLICHE:
            print(f"   {s + ' r' + str(r):11s}" + "".join(f"{v:6.2f}" for v in blocchi(serie[s, r][chiave])))

print("\n5. Equilibrazione automatica (Chodera 2016) e deriva lineare dopo t0")
print("   t0 vicino a 200 ns = la serie non e' stazionaria: ha fluttuazioni lente quanto la traiettoria")
for chiave in ("rmsd", "rb", "rg"):
    for s in SISTEMI:
        for r in REPLICHE:
            A = serie[s, r][chiave]; t = serie[s, r]["t"]
            t0, g, neff = timeseries.detect_equilibration(A, nskip=10)
            p, cov = np.polyfit(t[t0:], A[t0:], 1, cov=True)
            ci = 1.96 * np.sqrt(cov[0, 0] * g) * 100
            print(f"   {chiave:4s} {s:8s} r{r}: t0 {t[t0]:6.1f} ns  g {g:5.1f}  Neff {neff:5.0f}  "
                  f"deriva {p[0] * 100:+.2f} +/- {ci:.2f} A/100 ns")

# ---------------------------------------------------------------- confronto fra repliche
print("\n6. Confronto fra repliche (Ca 42-202, sovrapposizione sul core beta)")
m = (RES >= CORE[0]) & (RES <= CORE[1])
bsel = np.array([np.where(RES == k)[0][0] for k in STRAND])
for s in SISTEMI:
    Ys = []
    for r in REPLICHE:
        d = carica(s, r)
        ca = np.array([atomo(d, "CA", k) for k in RES])
        Ys.append(d["X"][:, ca].astype(np.float64)); M = d["ref"][ca].astype(np.float64)
    for _ in range(3):                                  # media iterativa comune alle tre repliche
        Ys = [fit(Y, M, bsel) for Y in Ys]; M = np.concatenate(Ys).mean(0)
    rmsf = np.array([np.sqrt(((Y - Y.mean(0)) ** 2).sum(2).mean(0)) for Y in Ys])
    NF = len(Ys[0]); H2 = NF // 2
    half = np.array([[np.sqrt(((Y[a:b] - Y[a:b].mean(0)) ** 2).sum(2).mean(0))
                      for a, b in ((0, H2), (H2, NF))] for Y in Ys])
    print(f"   == {s} ==")
    for r in range(3):
        print(f"   r{r + 1}: RMSF prima vs seconda meta'  r={np.corrcoef(half[r, 0, m], half[r, 1, m])[0, 1]:.2f}  "
              f"|diff| medio {np.abs(half[r, 0, m] - half[r, 1, m]).mean():.2f} A")
    for a, b in ((0, 1), (0, 2), (1, 2)):
        print(f"   r{a + 1} vs r{b + 1}: RMSF  r={np.corrcoef(rmsf[a, m], rmsf[b, m])[0, 1]:.2f}  "
              f"|diff| medio {np.abs(rmsf[a, m] - rmsf[b, m]).mean():.2f} A")
    sd = rmsf[:, m].std(0, ddof=1)
    top = np.argsort(sd)[::-1][:6]
    print(f"   SD dell'RMSF fra repliche: mediana {np.median(sd):.2f} A; residui peggiori: "
          + ", ".join(f"{RES[m][i]} ({sd[i]:.2f})" for i in top))
    Z = [Y[:, m].reshape(len(Y), -1) for Y in Ys]
    vecs, cc = [], []
    for Zr in Z:
        Zc = Zr - Zr.mean(0); w, v = np.linalg.eigh(Zc.T @ Zc); v = v[:, ::-1]
        vecs.append(v); p = Zc @ v[:, :3]
        cc.append([cosine_content(p[:, i], i + 1) for i in range(3)])
    for r in range(3):
        print(f"   r{r + 1}: cosine content PC1-PC3 = " + ", ".join(f"{c:.2f}" for c in cc[r]))
    for a, b in ((0, 1), (0, 2), (1, 2)):
        print(f"   RMSIP (10 PC) r{a + 1}-r{b + 1} = {rmsip(vecs[a], vecs[b]):.2f}")
    passo = int(round(1.0 / DT))                         # un frame ogni ns
    F = np.concatenate([Y[::passo][:, m] for Y in Ys])
    idx = np.arange(m.sum())
    D = np.array([np.sqrt(((fit(F, F[i], idx) - F[i]) ** 2).sum(2).mean(1)) for i in range(len(F))])
    lab = np.repeat([0, 1, 2], len(F) // 3)
    n3 = len(F) // 3
    intra = np.mean([D[np.ix_(lab == k, lab == k)][np.triu_indices(n3, 1)].mean() for k in range(3)])
    inter = np.mean([D[np.ix_(lab == a, lab == b)].mean() for a, b in ((0, 1), (0, 2), (1, 2))])
    print(f"   RMSD incrociato: dentro la replica {intra:.2f} A, fra repliche {inter:.2f} A, rapporto {inter / intra:.2f}")
    # --- come migliorano gli indicatori allungando la simulazione: finestre 0-T ---
    finestre = np.arange(100, int(round(n3 / 1.0)) + 1, 50)      # ns
    evo = []
    for T in finestre:
        k = int(T / DT)
        vv = []
        for Zr in Z:
            Zc = Zr[:k] - Zr[:k].mean(0); w, v = np.linalg.eigh(Zc.T @ Zc); vv.append(v[:, ::-1])
        rip = np.mean([rmsip(vv[a], vv[b]) for a, b in ((0, 1), (0, 2), (1, 2))])
        sel = np.concatenate([np.arange(T) + q * n3 for q in range(3)]).astype(int)
        Dw = D[np.ix_(sel, sel)]; lw = np.repeat([0, 1, 2], T)
        ia = np.mean([Dw[np.ix_(lw == q, lw == q)][np.triu_indices(T, 1)].mean() for q in range(3)])
        ie = np.mean([Dw[np.ix_(lw == a, lw == b)].mean() for a, b in ((0, 1), (0, 2), (1, 2))])
        rr = np.mean([np.corrcoef(np.sqrt(((Ys[a][:k] - Ys[a][:k].mean(0)) ** 2).sum(2).mean(0))[m],
                                  np.sqrt(((Ys[b][:k] - Ys[b][:k].mean(0)) ** 2).sum(2).mean(0))[m])[0, 1]
                      for a, b in ((0, 1), (0, 2), (1, 2))])
        evo.append((T, rip, ie / ia, rr))
    print("   finestra 0-T ns:  RMSIP medio | RMSD fra/dentro | r RMSF medio fra repliche")
    for T, rip, rap, rr in evo:
        print(f"     0-{T:3d}: {rip:.2f} | {rap:.2f} | {rr:.2f}")
    np.savez(os.path.join(CACHE, f"{s}_confronto.npz"), rmsf=rmsf, half=half, D=D, cc=np.array(cc), evo=np.array(evo))

# ---------------------------------------------------------------- dove divergono le repliche
print("\n7. Dove divergono le repliche")
SEGMENTI = [(44, 47), (78, 85), (102, 108), (134, 138), (170, 179), (195, 199)]
STRAND2 = [k for k in STRAND if k >= 57]          # core beta senza lo strand 44-47
bsel2 = np.array([np.where(RES == k)[0][0] for k in STRAND2])

def segmenti(res):
    res = sorted(res)
    if not res:
        return "nessuno"
    out = []; a = p = res[0]
    for k in res[1:]:
        if k > p + 2:
            out.append((a, p)); a = k
        p = k
    out.append((a, p))
    return ", ".join(f"{x}-{y}" if x != y else f"{x}" for x, y in out)

xr = carica("xray", 1)
XREF = xr["ref"][[atomo(xr, "CA", k) for k in RES]].astype(np.float64)
print("   7a. Strand per strand: RMSD dalla partenza e frazione beta (DSSP) nei primi 100 ns e negli ultimi 150 ns")
for s in SISTEMI:
    u = mda.Universe(os.path.join(ROOT, s, "conf.gro"))
    for r in REPLICHE:
        d = carica(s, r)
        ca = np.array([atomo(d, "CA", k) for k in RES])
        P = d["X"][::10][:, ca].astype(np.float64); ref = d["ref"][ca].astype(np.float64)
        Pa = fit(P, ref, bsel)
        u.load_new(d["X"][::50].copy(), format=mda.coordinates.memory.MemoryReader)
        e = DSSP(u.select_atoms("protein")).run().results.dssp == "E"
        riga = []
        for a, b in SEGMENTI:
            i = np.where((RES >= a) & (RES <= b))[0]
            dev = np.sqrt(((Pa[:, i] - ref[i]) ** 2).sum(2).mean(1))
            riga.append(f"{a}-{b} {dev[:100].mean():.1f}->{dev[-150:].mean():.1f} A, beta {e[:20, a - 42:b - 41].mean():.2f}->{e[-30:, a - 42:b - 41].mean():.2f}")
        Px = fit(P, XREF, np.where(m)[0])
        rx = np.sqrt(((Px[:, m] - XREF[m]) ** 2).sum(2).mean(1))
        print(f"   {s:8s} r{r}: " + " | ".join(riga))
        print(f"   {'':12s}RMSD Ca 44-199 dalla struttura X-ray: primi 50 ns {rx[:50].mean():.2f} A, ultimi 100 ns {rx[-100:].mean():.2f} A")

print("\n   7b. Confronto fra repliche senza lo strand 44-47 e il loop 48-56 (Ca 57-199, fit sul core beta senza 44-47)")
m2 = (RES >= 57) & (RES <= 199)
for s in SISTEMI:
    Ys = []
    for r in REPLICHE:
        d = carica(s, r)
        ca = np.array([atomo(d, "CA", k) for k in RES])
        Ys.append(d["X"][:, ca].astype(np.float64)); M = d["ref"][ca].astype(np.float64)
    for _ in range(3):
        Ys = [fit(Y, M, bsel2) for Y in Ys]; M = np.concatenate(Ys).mean(0)
    rmsf = np.array([np.sqrt(((Y - Y.mean(0)) ** 2).sum(2).mean(0)) for Y in Ys])
    rr = [np.corrcoef(rmsf[a, m2], rmsf[b, m2])[0, 1] for a, b in ((0, 1), (0, 2), (1, 2))]
    V = []
    for Y in Ys:
        Zc = Y[:, m2].reshape(len(Y), -1); Zc = Zc - Zc.mean(0); w, v = np.linalg.eigh(Zc.T @ Zc); V.append(v[:, ::-1])
    rp = [rmsip(V[a], V[b]) for a, b in ((0, 1), (0, 2), (1, 2))]
    passo = int(round(1.0 / DT))
    F = np.concatenate([Y[::passo][:, m2] for Y in Ys]); idx = np.arange(m2.sum()); n3 = len(F) // 3
    D2 = np.array([np.sqrt(((fit(F, F[i], idx) - F[i]) ** 2).sum(2).mean(1)) for i in range(len(F))])
    lab = np.repeat([0, 1, 2], n3)
    ia = np.mean([D2[np.ix_(lab == k, lab == k)][np.triu_indices(n3, 1)].mean() for k in range(3)])
    ie = np.mean([D2[np.ix_(lab == a, lab == b)].mean() for a, b in ((0, 1), (0, 2), (1, 2))])
    print(f"   {s:8s}: RMSF r {min(rr):.2f}-{max(rr):.2f} | RMSIP {min(rp):.2f}-{max(rp):.2f} | RMSD fra/dentro {ie:.2f}/{ia:.2f} = {ie / ia:.2f}")
    # posizione media di ogni replica nella seconda meta': dove le repliche stanno in stati diversi
    h = len(Ys[0]) // 2
    mu = np.array([Y[h:].mean(0) for Y in Ys])
    fl = np.array([np.sqrt(((Y[h:] - Y[h:].mean(0)) ** 2).sum(2).mean(0)) for Y in Ys]).mean(0)
    dmax = np.max([np.linalg.norm(mu[a] - mu[b], axis=1) for a, b in ((0, 1), (0, 2), (1, 2))], axis=0)
    dentro = (RES >= CORE[0]) & (RES <= CORE[1])
    print(f"   {'':10s}residui con posizioni medie di replica distanti > 2 A nella seconda meta': "
          f"{segmenti(list(RES[(dmax > 2.0) & dentro]))} ({int(((dmax > 2.0) & dentro).sum())} residui)")
    top = np.argsort(np.where(dentro, dmax, 0))[::-1][:5]
    print(f"   {'':10s}massimi: " + ", ".join(f"{RES[i]} {dmax[i]:.1f} A (fluttuazione interna {fl[i]:.1f} A)" for i in top))
