"""Fase 5: RMSF dei Ca nelle 3 repliche X-ray (MD classiche, 500 ns) contro i B-factor del cristallo.

Uso (dalla cartella del progetto mrkA, con le cache di 01_estrai.py):
    python analisi_convergenza/05_rmsf_vs_bfactor.py . analisi_convergenza/cache mrka_bfactor.pdb --shift -35 > bfactor/risultati_rmsf_vs_bfactor.txt
Serve la cache delle repliche X-ray (01_estrai.py). Il pdb e' numerato 1-182 e conf.gro 42-202 (numerazione wt): lo script
non applica da solo lo sfasamento, lo propone (pdb = wt - 35) e si ferma finche' non lo si passa con --shift.
Opzioni: --catena X (default: la catena che coincide con la sequenza di conf.gro), --shift N. La tabella per residuo va
in <nome del pdb>_per_residuo.csv, oppure nel file dato da CSV_OUT=... . Confronto solo qualitativo: correlazione di rango (Spearman) per residuo sui Ca e
rango percentile medio per segmento. Nessuna conversione B <-> RMSF (B = 8 pi^2/3 RMSF^2 non e' un accordo
quantitativo: i B contengono anche disordine statico, impaccamento cristallino e TLS).
Richiede numpy e scipy, e provenienza.py (copia di analisi_rest2/provenienza.py) nella stessa cartella dello script.
"""
import sys, os, json, argparse, numpy as np
from scipy.stats import spearmanr, rankdata

ap = argparse.ArgumentParser()
ap.add_argument("root"); ap.add_argument("cache"); ap.add_argument("pdb")
ap.add_argument("--catena", default=None); ap.add_argument("--shift", type=int, default=0)
ap.add_argument("--identico-a", dest="identico", default=None, help="altro pdb (es. mrka_bfactor.pdb) di cui si verifica che coincida con la catena scelta")
a = ap.parse_args()

QUI = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, QUI)
from provenienza import prov, testo, sha256            # stesso modulo degli script REST2 (analisi_rest2)
ins = [a.pdb] + ([a.identico] if a.identico else []) + [os.path.join(a.root, "xray", "conf.gro")] + [os.path.join(a.cache, f"xray_run{r}.npz") for r in (1, 2, 3)] + \
      [os.path.join(a.root, "traj", f"xray_run{r}_pr.xtc") for r in (1, 2, 3)]
estr = next((p for p in (os.path.join(QUI, "01_estrai.py"), os.path.join(QUI, "..", "01_estrai.py")) if os.path.exists(p)), None)
if estr:
    ins.append(estr)
P = prov(__file__, ins, vars(a), {"ROOT": a.root, "CACHE": a.cache})
for riga in testo(P):
    print(riga)
print("# la cache non porta la provenienza (01_estrai.py di analisi_convergenza non la scrive): "
      + (f"estratta, per quanto risulta, con 01_estrai.py sha256 {sha256(estr)[:16]}" if estr else "01_estrai.py non trovato accanto allo script, versione non nota"))
print()

RES = np.arange(42, 203)
CORE = (44, 199)           # esclusi 42-43 e 200-202 (cariche artificiali ai terminali), come negli altri script
DT = 0.1                   # ns fra due frame
T0 = 50.0                  # ns scartati, come in 04_confronto_ristretto.py
STRAND = [44, 45, 46, 47, 78, 79, 80, 81, 82, 83, 84, 85, 88, 94, 102, 103, 104, 105, 107, 108, 122, 134,
          135, 136, 137, 138, 148, 162, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 195, 196, 197, 198, 199]
LOOP = [(59, 77), (87, 101), (112, 118), (124, 133), (138, 147), (159, 169), (180, 192), (48, 56)]
BASE = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS", "CYX": "CYS", "CYM": "CYS",
        "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "MSE": "MET"}
bs = lambda r: BASE.get(r, r)

def fit(P, Q, sel):
    Ps = P[:, sel]; Pm = Ps.mean(1, keepdims=True); Qm = Q[sel].mean(0)
    H = np.einsum('fni,nj->fij', Ps - Pm, Q[sel] - Qm)
    U, S, Vt = np.linalg.svd(H)
    D = np.tile(np.eye(3), (len(P), 1, 1)); D[:, 2, 2] = np.sign(np.linalg.det(U @ Vt))
    return np.einsum('fni,fij->fnj', P - Pm, U @ D @ Vt) + Qm

# ------------------------------------------------------------------ 1. lettura del pdb
righe = open(a.pdb).read().splitlines()
print(f"1. File {os.path.basename(a.pdb)}: {len(righe)} righe")
cr = [l.rstrip() for l in righe if l.startswith("CRYST1")]
print("   " + (cr[0] if cr else "CRYST1 assente"))
chiavi = ("PROGRAM", "RESOLUTION RANGE HIGH", "(WORKING SET)", "FREE R VALUE  ", "FROM WILSON PLOT", "Individual atomic B", "Overall:", "Protein:",
          "Water:", "Chain  ")
for l in [l.rstrip() for l in righe if l.startswith("REMARK   3") and any(k in l for k in chiavi)][:14]:
    print("   " + l)
print(f"   occorrenze di 'TLS' nel file: {sum('TLS' in l for l in righe)}; REMARK 200 (raccolta dati, temperatura): "
      f"{'presente' if any(l.startswith('REMARK 200') for l in righe) else 'assente'}; REMARK 290 / 350: "
      f"{'presenti' if any(l.startswith(('REMARK 290', 'REMARK 350')) for l in righe) else 'assenti'}")
print(f"   MODEL: {sum(l.startswith('MODEL') for l in righe)}")
at = []
for l in righe:
    if l.startswith(("ATOM", "HETATM")):
        at.append(dict(rec=l[:6].strip(), name=l[12:16].strip(), alt=l[16].strip(), res=l[17:20].strip(), ch=l[21].strip(),
                       n=int(l[22:26]), ic=l[26].strip(), xyz=(float(l[30:38]), float(l[38:46]), float(l[46:54])),
                       occ=float(l[54:60]) if l[54:60].strip() else np.nan,
                       b=float(l[60:66]) if l[60:66].strip() else np.nan, el=l[76:78].strip()))
prot = [x for x in at if x["rec"] == "ATOM"]
chains = sorted({x["ch"] for x in prot})
het = [x for x in at if x["rec"] == "HETATM"]
print(f"   record ATOM: catene {chains}, insertion code {sum(bool(x['ic']) for x in prot)}; "
      f"HETATM: {len(het)} (" + ", ".join(f"{k} {sum(x['res'] == k for x in het)}" for k in sorted({x['res'] for x in het})) + "), ignorati: nessun HETATM entra nel confronto")
for c in chains:
    xs = [x for x in prot if x["ch"] == c]; bb = np.array([x["b"] for x in xs]); rr = sorted({x["n"] for x in xs})
    gap = [(q + 1, r - 1) for q, r in zip(rr, rr[1:]) if r - q > 1]
    print(f"   catena {c}: {len(xs)} atomi, {len(rr)} residui {rr[0]}-{rr[-1]}, gap {gap if gap else 'nessuno'}, altloc {sum(bool(x['alt']) for x in xs)}; "
          f"B min {bb.min():.1f}, media {bb.mean():.1f}, max {bb.max():.1f}; B = 0: {int((bb == 0).sum())}")

# ------------------------------------------------------------------ 2. numerazione contro conf.gro
gro = open(os.path.join(a.root, "xray", "conf.gro")).read().splitlines()
seq = {}
for l in gro[2:2 + int(gro[1])]:
    seq.setdefault(int(l[0:5]), bs(l[5:10].strip()))
print(f"\n2. Controllo di numerazione e catena contro xray/conf.gro (residui {min(seq)}-{max(seq)}, {len(seq)} residui)")
by_ch = {c: {} for c in chains}
for x in prot:
    if x["name"] == "CA" and x["alt"] in ("", "A"):
        by_ch[x["ch"]][x["n"]] = x

def identita(c, shift):
    ok = [k for k in seq if k + shift in by_ch[c] and bs(by_ch[c][k + shift]["res"]) == seq[k]]
    return len(ok)

best = max(((identita(c, s), c, s) for c in chains for s in range(-150, 151)), key=lambda t: t[0])
for c in chains:
    print(f"   catena {c}: {len(by_ch[c])} Ca; residui uguali a conf.gro con shift 0: {identita(c, 0)}/{len(seq)}")
print(f"   migliore (catena, shift) fra tutte: catena {best[1]}, shift {best[2]:+d}, {best[0]}/{len(seq)} residui uguali")
ch = a.catena or best[1]
sh = a.shift if a.catena or a.shift else best[2]
if a.shift == 0 and best[2] != 0 and not a.catena:
    print(f"   ATTENZIONE: la numerazione del pdb e' sfalsata di {best[2]:+d} rispetto a conf.gro. Non applico lo shift da solo:")
    print(f"   rilancia con --shift {best[2]} se e' voluto (e controlla l'header).")
    sys.exit(1)
nid = identita(ch, sh)
print(f"   uso catena {ch}, shift {sh:+d}: {nid}/{len(seq)} residui con lo stesso amminoacido")
manc = [k for k in seq if k + sh not in by_ch[ch]]
dif = [(k, seq[k], bs(by_ch[ch][k + sh]["res"])) for k in seq if k + sh in by_ch[ch] and bs(by_ch[ch][k + sh]["res"]) != seq[k]]
print(f"   residui di conf.gro assenti nel pdb: {manc if manc else 'nessuno'}")
print(f"   residui con amminoacido diverso: {dif if dif else 'nessuno'}")
if nid < len(seq) - len(manc):
    print("   ATTENZIONE: sequenza diversa; non vado avanti."); sys.exit(1)

if a.identico:
    mia = {l[12:27]: l[30:66] for l in righe if l.startswith("ATOM") and l[21] == ch}
    alt = {l[12:27]: l[30:66] for l in open(a.identico).read().splitlines() if l.startswith("ATOM")}
    dd = [k for k in mia if k in alt and mia[k] != alt[k]]
    print(f"   catena {ch} contro {os.path.basename(a.identico)}: {len(mia)} atomi contro {len(alt)}, atomi solo da una parte {len(set(mia) ^ set(alt))}, "
          f"con coordinate, occupanza o B diversi {len(dd)}: " + ("identica" if not dd and set(mia) == set(alt) else "DIVERSA"))
import re
m3 = [re.findall(r"[-+]?\d+\.?\d*", l.split(":", 1)[1]) for l in righe if l.startswith("REMARK   3") and re.search(rf"Chain\s+{ch}:", l)]
if m3:
    xs = [x["b"] for x in prot if x["ch"] == ch]
    print(f"   REMARK 3, catena {ch}: {m3[0][3]} atomi, B min {m3[0][0]}, max {m3[0][1]}, media {m3[0][2]}; nei record ATOM della catena: "
          f"{len(xs)} atomi, B min {min(xs):.2f}, max {max(xs):.2f}, media {np.mean(xs):.2f}"
          + ("" if int(m3[0][3]) == len(xs) and abs(float(m3[0][1]) - max(xs)) < 0.01 else ": NON coincidono (secondo francesco: atomi di un altro monomero o del tag di cristallizzazione tolti dopo il raffinamento; si usano i record ATOM)"))

# Ca per residuo: se il Ca ha altloc si prende quello con occupanza maggiore (il primo se pari); nessuna media fra altloc
B = np.full(len(RES), np.nan); OCC = np.full(len(RES), np.nan); XYZ = np.full((len(RES), 3), np.nan)
for i, k in enumerate(RES):
    cs = [x for x in prot if x["ch"] == ch and x["n"] == k + sh and x["name"] == "CA"]
    if cs:
        c0 = max(cs, key=lambda x: (x["occ"], -cs.index(x)))
        B[i] = c0["b"]; OCC[i] = c0["occ"]; XYZ[i] = c0["xyz"]
        if len(cs) > 1:
            print(f"   altloc del Ca {k} (pdb {k + sh}): " + ", ".join(f"{x['alt'] or '-'} occ {x['occ']:.2f} B {x['b']:.2f}" for x in cs)
                  + f" -> uso l'altloc {c0['alt'] or '-'}")
bad = (np.isnan(B) | (B <= 0) | (OCC < 0.5))
print(f"   Ca con B mancante, <= 0 o occupanza < 0.5 (esclusi dal confronto): "
      f"{[int(k) for k in RES[bad]] if bad.any() else 'nessuno'}")
print(f"   B del Ca 44-199: media {np.nanmean(B[(RES >= 44) & (RES <= 199)]):.1f}, "
      f"min {np.nanmin(B[(RES >= 44) & (RES <= 199)]):.1f}, max {np.nanmax(B[(RES >= 44) & (RES <= 199)]):.1f} A^2")

# ------------------------------------------------------------------ 3. le coordinate sono quelle della partenza?
ref = []
for r in (1, 2, 3):
    d = np.load(os.path.join(a.cache, f"xray_run{r}.npz"), allow_pickle=True)
    nm = d["names"].astype(str); rs = d["resids"]
    ca = np.array([np.where((nm == "CA") & (rs == k))[0][0] for k in RES])
    ref = d["ref"][ca].astype(np.float64)
    break
ok = ~np.isnan(XYZ[:, 0])
sel = np.where(ok)[0]
Pa = fit(XYZ[None], ref, sel)[0]
dev = np.linalg.norm(Pa - ref, axis=1)
print(f"\n3. Ca del pdb contro Ca della partenza X-ray (conf.gro), dopo sovrapposizione su {len(sel)} Ca")
print(f"   RMSD {np.sqrt((dev[sel] ** 2).mean()):.2f} A; mediana {np.median(dev[sel]):.2f}; massimo {dev[sel].max():.2f} A "
      f"al residuo {RES[sel][dev[sel].argmax()]}")
lontani = [(int(k), round(float(x), 2)) for k, x in zip(RES[sel], dev[sel]) if x > 0.5]
print(f"   residui con Ca > 0.5 A dalla partenza: {lontani if lontani else 'nessuno'}")

# ------------------------------------------------------------------ 4. RMSF delle tre repliche
m = (RES >= CORE[0]) & (RES <= CORE[1])
bsel = np.array([np.where(RES == k)[0][0] for k in STRAND])
print(f"\n4. RMSF dei Ca nelle repliche X-ray (fit su {len(STRAND)} Ca di strand, media comune iterata)")
Ys = []; Mref = None
for r in (1, 2, 3):
    d = np.load(os.path.join(a.cache, f"xray_run{r}.npz"), allow_pickle=True)
    nm = d["names"].astype(str); rs = d["resids"]
    ca = np.array([np.where((nm == "CA") & (rs == k))[0][0] for k in RES])
    Ys.append(d["X"][:, ca].astype(np.float64)); Mref = d["ref"][ca].astype(np.float64)
    print(f"   replica {r}: {len(Ys[-1])} frame, {d['t'][0] / 1000:.1f}-{d['t'][-1] / 1000:.1f} ns")
FIN = {}
for nome, t0 in (("50-500 ns", T0), ("0-500 ns", 0.0)):
    k0 = int(round(t0 / DT))
    Z = [Y[k0:] for Y in Ys]; M = Mref
    for _ in range(3):
        Z = [fit(Y, M, bsel) for Y in Z]; M = np.concatenate(Z).mean(0)
    FIN[nome] = np.array([np.sqrt(((Y - Y.mean(0)) ** 2).sum(2).mean(0)) for Y in Z])
    print(f"   finestra {nome}: RMSF mediano {np.median(FIN[nome][:, m]):.2f} A, min {FIN[nome][:, m].min():.2f}, max {FIN[nome][:, m].max():.2f}")
    print("      r fra repliche (Ca 44-199): " + ", ".join(
        f"r{i + 1}-r{j + 1} {np.corrcoef(FIN[nome][i, m], FIN[nome][j, m])[0, 1]:.2f}" for i, j in ((0, 1), (0, 2), (1, 2))))
print("   (controllo: sulla finestra 0-500 ns i tre r devono coincidere con 0.82, 0.94, 0.89 di risultati_convergenza_500ns.txt)")

# ------------------------------------------------------------------ 5. correlazione di rango
def pct(x):                       # rango percentile 0-1 fra i residui validi
    return (rankdata(x) - 0.5) / len(x)

print("\n5. Correlazione di rango (Spearman) fra B-factor del Ca e RMSF dei Ca, residui 44-199 validi")
sott = [("tutti 44-199", (RES >= 44) & (RES <= 199)),
        ("senza 44-56 (vicino al taglio N-terminale)", (RES >= 57) & (RES <= 199)),
        ("senza i 5 residui per estremo", (RES >= 49) & (RES <= 194)),
        ("solo i loop (7 segmenti + 48-56)", np.isin(RES, [k for lo, hi in LOOP for k in range(lo, hi + 1)])),
        ("solo core beta (43 Ca di strand)", np.isin(RES, STRAND))]
print(f"   {'sottoinsieme':46s} {'n':>4s}  " + "  ".join(f"{n:>14s}" for n in ("media 3 rep.", "r1", "r2", "r3")))
for fin in FIN:
    print(f"   -- RMSF su {fin}")
    for nome, ms in sott:
        mm = ms & ~bad
        riga = []
        for v in (FIN[fin].mean(0), FIN[fin][0], FIN[fin][1], FIN[fin][2]):
            rho = spearmanr(B[mm], v[mm]).statistic
            riga.append(f"{rho:14.2f}")
        print(f"   {nome:46s} {int(mm.sum()):4d}  " + "  ".join(riga))
print("   (intervalli di confidenza nella sezione 6; il B non e' una sonda affidabile della mobilita', vedi sezione 7)")

# ------------------------------------------------------------------ 6. intervallo di confidenza (bootstrap a blocchi)
print("\n6. Intervallo di confidenza al 95% di Spearman: bootstrap a blocchi di residui consecutivi (lunghezza 5, 10 e 20; 5000 ricampionamenti),")
print("   RMSF = media delle 3 repliche su 50-500 ns. I residui vicini non sono indipendenti e con ~15 blocchi l'intervallo dipende dalla")
print("   lunghezza del blocco: si riporta il piu' largo")
rng = np.random.default_rng(1)
rm = FIN["50-500 ns"].mean(0)
def rho(x, y):
    return np.corrcoef(rankdata(x), rankdata(y))[0, 1]
IC = {}
for nome, ms in sott[:2] + sott[4:]:
    idx = np.where(ms & ~bad)[0]
    n = len(idx); r = {}
    for L in (5, 10, 20):
        nb = int(np.ceil(n / L)); out = []
        for _ in range(5000):
            st = rng.integers(0, n, nb)
            ii = idx[np.concatenate([(q + np.arange(L)) % n for q in st])[:n]]
            out.append(rho(B[ii], rm[ii]))
        r[L] = np.percentile(out, [2.5, 97.5])
    w = max(r, key=lambda L: r[L][1] - r[L][0])
    IC[nome] = dict(rho=float(rho(B[idx], rm[idx])), n=int(n), blocchi={str(L): [float(x) for x in r[L]] for L in r}, piu_largo=int(w))
    print(f"   {nome:46s} rho {rho(B[idx], rm[idx]):5.2f}   " + "   ".join(f"L={L:2d} [{r[L][0]:5.2f}, {r[L][1]:5.2f}]" for L in r)
          + f"   -> piu' largo (L={w}): [{r[w][0]:5.2f}, {r[w][1]:5.2f}]")

# ------------------------------------------------------------------ 7. plausibilita' del file di B-factor
print("\n7. Controllo di plausibilita': B-factor e RMSF contro due proxy strutturali (Ca della partenza, Ca 44-199 validi)")
print("   distanza del Ca dal centroide e numero di Ca entro 10 A (stesse 161 posizioni della MD; il pdb ha anche lo strand donatore)")
cen = ref[(RES >= 44) & (RES <= 199)].mean(0)
dist = np.linalg.norm(ref - cen, axis=1)
nco = np.array([((np.linalg.norm(ref - ref[i], axis=1) < 10.0).sum() - 1) for i in range(len(RES))])
mv2 = (RES >= 44) & (RES <= 199) & ~bad
print(f"   {'':34s} {'distanza dal centroide':>24s} {'Ca entro 10 A':>16s}")
for nome, v in (("B-factor del Ca (cristallo)", B), ("RMSF medio delle 3 repliche (MD)", rm)):
    print(f"   {nome:34s} {spearmanr(v[mv2], dist[mv2]).statistic:24.2f} {spearmanr(v[mv2], nco[mv2]).statistic:16.2f}")
allca = np.array([x["xyz"] for c in (by_ch[ch],) for _, x in sorted(c.items())])
cen2 = allca.mean(0)
dist2 = np.linalg.norm(XYZ - cen2, axis=1)
nco2 = np.array([((np.linalg.norm(allca - XYZ[i], axis=1) < 10.0).sum() - 1) if ok[i] else 0 for i in range(len(RES))])
print(f"   stessi proxy contando tutti i {len(allca)} Ca del pdb (anche lo strand al C-terminale e i residui 36-41):")
for nome, v in (("B-factor del Ca (cristallo)", B), ("RMSF medio delle 3 repliche (MD)", rm)):
    print(f"   {nome:34s} {spearmanr(v[mv2], dist2[mv2]).statistic:24.2f} {spearmanr(v[mv2], nco2[mv2]).statistic:16.2f}")
print("   Il B del cristallo non segue la burial (rho vicino a 0), l'RMSF delle MD si': in un cristallo ci si aspetta che il B cresca con la")
print("   distanza dal centro (atteso, non verificato qui su altri cristalli). Quindi il B di questo file non e' una sonda affidabile della")
print("   mobilita' locale e va letto cosi' prima di qualsiasi confronto per segmento. Cause possibili, non scelte: contatti di reticolo,")
print("   partner del complesso (le altre catene del cristallo, assenti nelle MD), strand donatore. Il file non riporta TLS: i B sono individuali isotropi.")

# ------------------------------------------------------------------ 8. per segmento
print("\n8. Rango percentile medio per segmento (0 = piu' rigido, 1 = piu' mobile fra i Ca 44-199 validi)")
mv = m & ~bad
pB = np.full(len(RES), np.nan); pR = np.full(len(RES), np.nan)
pB[mv] = pct(B[mv]); pR[mv] = pct(FIN["50-500 ns"].mean(0)[mv])
print(f"   {'segmento':10s} {'n':>3s} {'B medio':>8s} {'rango B':>8s} {'RMSF medio':>11s} {'rango RMSF':>11s} {'differenza':>11s}   RMSF r1 / r2 / r3 (A)")
for lo, hi in LOOP + [(78, 85), (102, 108), (134, 137), (170, 179), (195, 199)]:
    ms = (RES >= lo) & (RES <= hi) & mv
    if not ms.any():
        continue
    print(f"   {lo:3d}-{hi:<6d} {int(ms.sum()):3d} {np.mean(B[ms]):8.1f} {np.mean(pB[ms]):8.2f} {FIN['50-500 ns'].mean(0)[ms].mean():11.2f} "
          f"{np.mean(pR[ms]):11.2f} {np.mean(pR[ms]) - np.mean(pB[ms]):+11.2f}   "
          + " / ".join(f"{FIN['50-500 ns'][r][ms].mean():.2f}" for r in range(3)))
print("   (descrittivo: dato che il B non segue la burial, sezione 7, non e' una misura di accordo)")
print("   differenza = rango RMSF - rango B: positiva = la MD e' piu' mobile del cristallo rispetto al resto della proteina")

# ------------------------------------------------------------------ 9. dettaglio 134-150
rm_seg = FIN["50-500 ns"].mean(0)
print("\n9. Dettaglio per residuo 134-150 (RMSF su 50-500 ns)")
print(f"   {'res':>4s} {'aa':>4s} {'B Ca':>6s} {'rango B':>8s} {'RMSF':>6s} {'rango RMSF':>11s}   RMSF r1 / r2 / r3")
for i, k in enumerate(RES):
    if 134 <= k <= 150:
        print(f"   {k:4d} {seq[k]:>4s} {B[i]:6.1f} {pB[i]:8.2f} {FIN['50-500 ns'].mean(0)[i]:6.2f} {pR[i]:11.2f}   "
              + " / ".join(f"{FIN['50-500 ns'][r][i]:.2f}" for r in range(3)))
s1 = (RES >= 138) & (RES <= 147) & mv
print(f"   Spearman B vs RMSF dentro 138-147: {spearmanr(B[s1], rm_seg[s1]).statistic:.2f} (n = {int(s1.sum())}: solo descrittivo, non e' un test)")
s2 = (RES >= 143) & (RES <= 146)
print(f"   RMSF in 143-146 per replica: r1 {FIN['50-500 ns'][0][s2].min():.2f}-{FIN['50-500 ns'][0][s2].max():.2f}, "
      f"r2 {FIN['50-500 ns'][1][s2].min():.2f}-{FIN['50-500 ns'][1][s2].max():.2f}, "
      f"r3 {FIN['50-500 ns'][2][s2].min():.2f}-{FIN['50-500 ns'][2][s2].max():.2f} A: le repliche non concordano in questo tratto")

# ------------------------------------------------------------------ 10. tabella per residuo
csv = os.path.splitext(os.path.abspath(sys.argv[3] if len(sys.argv) > 3 else "x"))[0] + "_per_residuo.csv"
csv = os.environ.get("CSV_OUT", csv)
with open(csv, "w") as f:
    f.write("resid_wt,resid_pdb,aa,B_Ca,rango_B,RMSF_r1,RMSF_r2,RMSF_r3,RMSF_media,rango_RMSF\n")
    for i, k in enumerate(RES):
        if CORE[0] <= k <= CORE[1]:
            f.write(f"{k},{k + sh},{seq[k]},{B[i]:.2f},{pB[i]:.3f},{FIN['50-500 ns'][0][i]:.3f},{FIN['50-500 ns'][1][i]:.3f},"
                    f"{FIN['50-500 ns'][2][i]:.3f},{rm[i]:.3f},{pR[i]:.3f}\n" if mv[i] else
                    f"{k},{k + sh},{seq[k]},,,{FIN['50-500 ns'][0][i]:.3f},{FIN['50-500 ns'][1][i]:.3f},{FIN['50-500 ns'][2][i]:.3f},{rm[i]:.3f},\n")
print(f"\n10. Tabella per residuo (44-199) scritta in {os.path.basename(csv)} (B e rango B vuoti dove il Ca non e' nel pdb)")

# ------------------------------------------------------------------ 11. limiti e provenienza dei file scritti
print("\n11. Limiti da riportare insieme ai risultati")
print("   - il B non segue la burial: non e' una sonda affidabile della mobilita' locale (sezione 7); cause possibili non scelte:")
print("     contatti di reticolo, partner del complesso (le altre catene del cristallo), strand donatore")
print("   - i contatti con le altre catene e con le copie di simmetria sono in 06_contatti_cristallo.py; gli atomi tolti dopo il raffinamento (altro monomero o tag) non ci sono, quindi i contatti possono essere sottostimati")
print("   - temperatura di raccolta del cristallo non riportata (nessun REMARK 200); le MD sono a 310 K")
print("   - le statistiche di REMARK 3 includono atomi di un altro monomero o del tag di cristallizzazione tolti dopo il raffinamento (secondo francesco, vedi sezione 2): i B usati sono quelli dei record ATOM")
print("   - le MD sono senza partner e senza strand donatore (23-42 rimosso)")
print("   - 3 repliche per sistema; nessuna conversione B <-> RMSF (B = 8 pi^2/3 RMSF^2 non e' un accordo quantitativo)")
P["output"] = {os.path.basename(csv): sha256(csv)}
P["intervalli_confidenza"] = IC
js = os.path.splitext(csv)[0] + "_provenienza.json"
with open(js, "w") as f:
    json.dump(P, f, indent=1, ensure_ascii=False)
print(f"\nProvenienza completa (sha256 interi) in {os.path.basename(js)}")
