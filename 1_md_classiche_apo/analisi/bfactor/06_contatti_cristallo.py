"""Fase 6: contatti della catena A (MrkA) nel cristallo, per i tratti 87-101 e 138-147 e per il resto della proteina.

Uso:  python 06_contatti_cristallo.py originale/mrka-Fab.pdb [--shift -35] [--catena A] [--csv-out file.csv]
Criteri fissati dal revisore prima di vedere i dati (29 ott 2026):
  - file: originale con header; catena A = MrkA (numerazione wt = pdb + 35); catene B e C = partner nella stessa unita' asimmetrica;
    acqua e HETATM esclusi; per gli altloc si prende quello a occupanza maggiore (il primo se pari);
  - copie di simmetria da CRYST1 (C 1 2 1) generate qui con numpy, senza gemmi (non installato): operatori di C2 e traslazioni di cella,
    con tre controlli (matrice SCALE, assenza di sovrapposizioni, ricalcolo diretto dei contatti nell'unita' asimmetrica);
  - contatto = due atomi pesanti entro 4.0 A; sensibilita' a 3.5 e 5.0 A;
  - per ogni residuo 44-199 della catena A: contatto si'/no separato per partner: Fab nella stessa unita' asimmetrica (B, C),
    copie di simmetria di MrkA, copie di simmetria del Fab; riferimento = frazione di residui in contatto nel resto della proteina;
  - il B dei Ca in contatto e non in contatto e' solo descrittivo, senza test.
Aggiunta non prevista dai criteri e dichiarata: contatti con i residui della stessa catena A fuori dal costrutto 42-202 (pdb 1-6 e 168-182,
l'N-terminale e lo strand al C-terminale), che le MD non hanno.
Richiede numpy e scipy, e provenienza.py nella stessa cartella.
"""
import sys, os, json, argparse, itertools, numpy as np
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("pdb"); ap.add_argument("--catena", default="A"); ap.add_argument("--shift", type=int, default=-35)
ap.add_argument("--csv-out", default=None)
a = ap.parse_args()

QUI = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, QUI)
from provenienza import prov, testo, sha256
P = prov(__file__, [a.pdb], vars(a), {"PDB": os.path.dirname(os.path.abspath(a.pdb))})
for riga in testo(P):
    print(riga)
print()

CUTS = (3.5, 4.0, 5.0); CUT = 4.0
SEG = {"87-101": (87, 101), "138-147": (138, 147)}
CORE = (44, 199)
CAT = ("fab_au", "mrka_sim", "fab_sim", "mrka_fuori")       # mrka_fuori: aggiunta dichiarata, stessa catena A fuori da 42-202
NOMI = {"fab_au": "Fab, stessa unita' asimmetrica (B, C)", "mrka_sim": "copie di simmetria di MrkA",
        "fab_sim": "copie di simmetria del Fab", "mrka_fuori": "catena A fuori da 42-202 (pdb 1-6, 168-182; aggiunta)"}

# ------------------------------------------------------------------ 1. lettura: solo ATOM, un altloc per atomo
righe = open(a.pdb).read().splitlines()
cr = [l for l in righe if l.startswith("CRYST1")]
assert cr, "CRYST1 assente"
c = cr[0]
cella = [float(c[6:15]), float(c[15:24]), float(c[24:33]), float(c[33:40]), float(c[40:47]), float(c[47:54])]
sg = c[55:66].strip()
print(f"1. {os.path.basename(a.pdb)}: CRYST1 a {cella[0]} b {cella[1]} c {cella[2]}, alfa {cella[3]} beta {cella[4]} gamma {cella[5]}, gruppo {sg}")
OPS = {"C 1 2 1": [(np.diag([1, 1, 1]), (0, 0, 0)), (np.diag([-1, 1, -1]), (0, 0, 0)),
                   (np.diag([1, 1, 1]), (.5, .5, 0)), (np.diag([-1, 1, -1]), (.5, .5, 0))]}
assert sg in OPS, f"gruppo {sg} non previsto"
ops = OPS[sg]
al, be, ga = np.radians(cella[3:])
V = cella[0] * cella[1] * cella[2] * np.sqrt(1 - np.cos(al) ** 2 - np.cos(be) ** 2 - np.cos(ga) ** 2 + 2 * np.cos(al) * np.cos(be) * np.cos(ga))
O = np.array([[cella[0], cella[1] * np.cos(ga), cella[2] * np.cos(be)],
              [0, cella[1] * np.sin(ga), cella[2] * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)],
              [0, 0, V / (cella[0] * cella[1] * np.sin(ga))]])
Oi = np.linalg.inv(O)
sc = np.array([[float(x) for x in l.split()[1:4]] for l in righe if l.startswith(("SCALE1", "SCALE2", "SCALE3"))])
print(f"   controllo 1: inversa della matrice di ortogonalizzazione contro SCALE1-3 del file: differenza massima {np.abs(Oi - sc).max():.1e}"
      + (" (ok)" if np.abs(Oi - sc).max() < 2e-6 else " (NON coincidono)"))
print(f"   operatori del gruppo {sg} usati: " + "; ".join(f"R={R.diagonal().astype(int).tolist()} t={t}" for R, t in ops))

brut = {}
for l in righe:
    if l.startswith("ATOM"):
        el = l[76:78].strip() or l[12:16].strip()[0]
        if el in ("H", "D"):
            continue
        k = (l[21], int(l[22:26]), l[26], l[12:16].strip())
        occ = float(l[54:60]) if l[54:60].strip() else 1.0
        if k not in brut or occ > brut[k][0]:                 # occupanza maggiore; a parita' resta il primo
            brut[k] = (occ, l)
AT = [v[1] for v in brut.values()]
ch = np.array([l[21] for l in AT]); rn = np.array([int(l[22:26]) for l in AT]); rname = np.array([l[17:20] for l in AT])
an = np.array([l[12:16].strip() for l in AT]); X = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in AT])
Bf = np.array([float(l[60:66]) for l in AT])
print(f"   atomi pesanti ATOM dopo la scelta degli altloc: " + ", ".join(f"catena {q} {int((ch == q).sum())}" for q in sorted(set(ch)))
      + "; HETATM (acqua) esclusi")
F = X @ Oi.T                                               # frazionarie

# ------------------------------------------------------------------ 2. immagini per simmetria e traslazione di cella
def immagini(raggio):
    """(indice operatore, traslazione, coordinate cartesiane) per tutte le immagini che entrano nel riquadro della catena A +/- raggio."""
    qa = ch == a.catena
    lo, hi = X[qa].min(0) - raggio, X[qa].max(0) + raggio
    out = []
    for io, (R, t) in enumerate(ops):
        Fo = F @ R.T + np.array(t)
        for tr in itertools.product(range(-4, 5), repeat=3):
            Xi = (Fo + np.array(tr)) @ O.T
            vicini = np.all((Xi >= lo) & (Xi <= hi), axis=1)
            if vicini.any():
                out.append((io, tr, Xi))
    return out

IM = immagini(max(CUTS))
print(f"   immagini (operatore, traslazione di cella) con atomi nel riquadro della catena {a.catena} +/- {max(CUTS)} A: {len(IM)}")

# controllo 2: nessuna sovrapposizione fra unita' asimmetrica e le sue copie
tauto = cKDTree(X)
dmin = np.inf; n22 = 0
for io, tr, Xi in IM:
    if io == 0 and tr == (0, 0, 0):
        continue
    d, _ = tauto.query(Xi, k=1)
    dmin = min(dmin, d.min()); n22 += int((d < 2.2).sum())
print(f"   controllo 2: distanza minima fra un atomo dell'unita' asimmetrica e una sua copia di simmetria/traslazione: {dmin:.2f} A; "
      f"atomi a meno di 2.2 A: {n22}" + (" (ok)" if dmin > 2.2 and n22 == 0 else " (SOVRAPPOSIZIONI: operatori sbagliati)"))
print("   (le sovrapposizioni attese con operatori giusti sono nulle: legami idrogeno e contatti vdW restano sopra 2.4 A)")

# ------------------------------------------------------------------ 3. contatti
qa = np.where(ch == a.catena)[0]
wt = rn + (-a.shift)                                      # numerazione wt
res_q = sorted({int(w) for w in wt[qa] if CORE[0] <= w <= CORE[1]})
fuori_costrutto = (wt < 42) | (wt > 202)                  # atomi della catena A che le MD non hanno
tq = cKDTree(X[qa])
MIN = {cat: {w: np.inf for w in res_q} for cat in CAT}
PART = {cat: {w: {} for w in res_q} for cat in CAT}       # partner entro CUT: {(catena, residuo): distanza minima}
for io, tr, Xi in IM:
    identita = io == 0 and tr == (0, 0, 0)
    for q in sorted(set(ch)):
        sel = np.where(ch == q)[0]
        if q == a.catena:
            if identita:
                sel = sel[fuori_costrutto[sel]]
                cat = "mrka_fuori"
            else:
                cat = "mrka_sim"
        else:
            cat = "fab_au" if identita else "fab_sim"
        if len(sel) == 0:
            continue
        tp = cKDTree(Xi[sel])
        liste = tq.query_ball_tree(tp, max(CUTS))
        for iq, lst in enumerate(liste):
            if not lst:
                continue
            w = int(wt[qa[iq]])
            if w not in MIN[cat]:
                continue
            d = np.linalg.norm(Xi[sel[lst]] - X[qa[iq]], axis=1)
            MIN[cat][w] = min(MIN[cat][w], d.min())
            for j, dj in zip(sel[lst], d):
                if dj <= CUT:
                    kk = (q, int(rn[j]), str(rname[j]), io, tr)
                    PART[cat][w][kk] = min(PART[cat][w].get(kk, 99), float(dj))

# controllo 3: ricalcolo diretto (tutte le coppie, senza alberi) dei contatti catena A - Fab nell'unita' asimmetrica
dd = np.linalg.norm(X[qa][:, None, :] - X[ch != a.catena][None, :, :], axis=2)
dir_au = {w: float(dd[wt[qa] == w].min()) for w in res_q}
diff = 0.0; stonati = []
for w in res_q:
    alb = MIN["fab_au"][w]
    if dir_au[w] <= max(CUTS) or alb <= max(CUTS):
        if abs(dir_au[w] - alb) > 1e-6:
            stonati.append(w)
        else:
            diff = max(diff, abs(dir_au[w] - alb))
print(f"\n3. Controllo 3: distanza minima catena A - Fab nell'unita' asimmetrica, alberi contro tutte le coppie, residui entro {max(CUTS)} A: "
      f"{sum(1 for w in res_q if dir_au[w] <= max(CUTS))}; residui con risultati diversi: {stonati if stonati else 'nessuno'}")

# ------------------------------------------------------------------ 4. B dei Ca
bca = {int(wt[i]): float(Bf[i]) for i in qa if an[i] == "CA"}
res_valid = [w for w in res_q if w in bca]
assert len(res_valid) == len(res_q), f"Ca mancanti: {[w for w in res_q if w not in bca]}"
mancanti = [w for w in range(CORE[0], CORE[1] + 1) if w not in res_q]
print(f"   residui 44-199 senza atomi nella catena {a.catena} (gap, esclusi senza interpolare): {mancanti}")

def flag(cat, w, cut):
    return MIN[cat][w] <= cut

def qualunque(w, cut, cats=("fab_au", "mrka_sim", "fab_sim")):
    return any(flag(cat, w, cut) for cat in cats)

segs = {n: [w for w in res_q if lo <= w <= hi] for n, (lo, hi) in SEG.items()}
resto = [w for w in res_q if not any(lo <= w <= hi for lo, hi in SEG.values())]
gruppi = dict(segs); gruppi["resto 44-199"] = resto; gruppi["tutti 44-199"] = res_q

print("\n4. Frazione di residui con almeno un contatto (n in contatto / n residui), per partner")
print("   'qualunque' = Fab nella stessa unita' asimmetrica, copie di simmetria di MrkA o del Fab (la categoria aggiunta non conta)")
for cut in CUTS:
    print(f"\n   soglia {cut} A" + ("  <- riferimento" if cut == CUT else "  (sensibilita')"))
    print(f"   {'gruppo':14s} {'n':>3s} " + " ".join(f"{c:>12s}" for c in ("qualunque",) + CAT))
    for nome, g in gruppi.items():
        riga = [f"{sum(qualunque(w, cut) for w in g):3d} ({100 * sum(qualunque(w, cut) for w in g) / len(g):3.0f}%)"]
        for cat in CAT:
            n = sum(flag(cat, w, cut) for w in g)
            riga.append(f"{n:3d} ({100 * n / len(g):3.0f}%)")
        print(f"   {nome:14s} {len(g):3d} " + " ".join(f"{r:>12s}" for r in riga))

print(f"\n5. Dettaglio per residuo nei due tratti (soglia {CUT} A; distanza minima in A per partner, '-' oltre {max(CUTS)} A)")
for nome in SEG:
    print(f"\n   {nome}")
    print(f"   {'wt':>4s} {'aa':>4s} {'B Ca':>6s}  " + " ".join(f"{c:>10s}" for c in CAT) + "   partner entro 4.0 A (catena residuo, operatore, traslazione)")
    for w in segs[nome]:
        i = [k for k in qa if wt[k] == w][0]
        dist = [f"{MIN[cat][w]:10.2f}" if MIN[cat][w] < 99 else f"{'-':>10s}" for cat in CAT]
        pt = []
        for cat in CAT[:3]:
            for (q, r, rnm, io, tr), dj in sorted(PART[cat][w].items(), key=lambda t: t[1])[:3]:
                pt.append(f"{q}{r} {rnm} op{io}{list(tr)} {dj:.2f}")
        print(f"   {w:4d} {rname[i]:>4s} {bca[w]:6.1f}  " + " ".join(dist) + "   " + "; ".join(pt[:4]))

print(f"\n6. Partner per tratto (soglia {CUT} A): catene e residui di contatto")
for nome in list(SEG) + ["resto 44-199"]:
    for cat in CAT[:3]:
        ps = {}
        for w in gruppi[nome]:
            for (q, r, rnm, io, tr) in PART[cat][w]:
                ps.setdefault(q if cat != "mrka_sim" else "A", set()).add((r, io, tr))
        riga = ", ".join(f"{q}: {len({r for r, _, _ in v})} residui in {len({(io, tr) for _, io, tr in v})} copie" for q, v in sorted(ps.items()))
        print(f"   {nome:14s} {NOMI[cat]:42s} {riga if riga else 'nessuno'}")

print("\n7. B-factor dei Ca (A^2), solo descrittivo, senza test; contatto = 'qualunque' a 4.0 A")
for nome, g in gruppi.items():
    si = [bca[w] for w in g if qualunque(w, CUT)]; no = [bca[w] for w in g if not qualunque(w, CUT)]
    f = lambda v: f"n {len(v):3d}, media {np.mean(v):5.1f}, mediana {np.median(v):5.1f}" if v else "n   0"
    print(f"   {nome:14s} in contatto: {f(si):38s} non in contatto: {f(no)}")

# ------------------------------------------------------------------ 8. limiti
print("\n8. Limiti da riportare insieme ai risultati")
print("   - contatto = prossimita' di atomi pesanti nel cristallo (4.0 A, sensibilita' 3.5 e 5.0): non misura energia ne' superficie sepolta,")
print("     e non dice se il contatto stabilizza il tratto")
print("   - gli atomi tolti dopo il raffinamento (secondo francesco appartenevano a un altro monomero o al tag di cristallizzazione) erano nel")
print("     cristallo e i loro contatti non sono contati: i contatti possono essere sottostimati")
print("   - copie di simmetria generate con numpy dagli operatori di C 1 2 1 (gemmi non installato), con i tre controlli sopra; acqua esclusa")
print("   - i B dei Ca sono quelli dei record ATOM; nessuna conclusione su epitopo o sull'obiettivo 2; mai 'causato da'")

# ------------------------------------------------------------------ 9. csv
csv = a.csv_out or os.path.splitext(os.path.abspath(a.pdb))[0] + "_contatti_per_residuo.csv"
with open(csv, "w") as f:
    f.write("resid_wt,resid_pdb,aa,B_Ca,tratto," + ",".join(f"d_min_{c}" for c in CAT) + "," +
            ",".join(f"{c}_{str(u).replace('.', '')}" for u in CUTS for c in ("qualunque",) + CAT) + "\n")
    for w in res_q:
        i = [k for k in qa if wt[k] == w][0]
        tr = next((n for n, (lo, hi) in SEG.items() if lo <= w <= hi), "resto")
        f.write(f"{w},{w + a.shift},{rname[i]},{bca[w]:.2f},{tr}," + ",".join(f"{MIN[c][w]:.2f}" if MIN[c][w] < 99 else "" for c in CAT) + "," +
                ",".join(str(int(qualunque(w, u))) if c == "qualunque" else str(int(flag(c, w, u))) for u in CUTS for c in ("qualunque",) + CAT) + "\n")
print(f"\n9. Tabella per residuo scritta in {os.path.basename(csv)}")
P["output"] = {os.path.basename(csv): sha256(csv)}
js = os.path.splitext(csv)[0] + "_provenienza.json"
with open(js, "w") as f:
    json.dump(P, f, indent=1, ensure_ascii=False)
print(f"   provenienza completa in {os.path.basename(js)}")
