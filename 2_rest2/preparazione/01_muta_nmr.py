"""Porta il medoide NMR ossidato (9HW9, modello 6) alla sequenza del costrutto X-ray.

Uso (dalla cartella del progetto mrkA):
    python rest2/01_muta_nmr.py mrka_nmr_medoide_ox_nocap.pdb xray/conf.gro rest2/mrka_nmr_seqxray_ox_nocap.pdb

Le 8 mutazioni NMR -> X-ray si fanno senza modellare rotameri:
- N86D: isosterica, ND2 diventa OD2;
- T106A, S128A, T143A, K160A: la catena laterale si tronca al CB;
- L145V, A167Q, A182T: la catena laterale X-ray dello stesso residuo viene innestata
  sovrapponendo N, CA, C X-ray a quelli NMR, cosi' il rotamero e' quello cristallografico.
Il backbone NMR non viene toccato. Solo atomi pesanti, CYX 43/87 e terminali liberi come
nel file di partenza.
"""
import sys, numpy as np

NMR, XRAY, OUT = sys.argv[1:4]
TRONCA = {106: "THR", 128: "SER", 143: "THR", 160: "LYS"}
INNESTA = {145: ("LEU", "VAL"), 167: ("ALA", "GLN"), 182: ("ALA", "THR")}
ORDINE = {"ALA": ["N", "CA", "C", "O", "CB"], "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
          "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"], "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
          "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"]}

def leggi_pdb(path):
    atomi = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            atomi.append({"name": l[12:16].strip(), "resn": l[17:20].strip(), "chain": l[21],
                          "resi": int(l[22:26]), "xyz": np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]),
                          "elem": l[76:78].strip()})
    return atomi

def leggi_gro(path):
    righe = open(path).read().splitlines()
    atomi = []
    for l in righe[2:2 + int(righe[1])]:
        atomi.append({"name": l[10:15].strip(), "resn": l[5:10].strip(), "resi": int(l[0:5]),
                      "xyz": 10.0 * np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])])})
    return atomi

def kabsch(P, Q):
    """Rotazione R e traslazione t che portano P su Q (righe = punti)."""
    pm, qm = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - pm).T @ (Q - qm))
    D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))])
    R = U @ D @ Vt
    return R, qm - pm @ R

def diedro(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))

def ruota(p, origine, asse, gradi):
    """Rotazione di Rodrigues del punto p attorno all'asse che passa per origine."""
    k = asse / np.linalg.norm(asse); v = p - origine; t = np.radians(gradi)
    return origine + v * np.cos(t) + np.cross(k, v) * np.sin(t) + k * (k @ v) * (1 - np.cos(t))

def distanza_minima(sc, ambiente):
    """Distanza minima corretta per il tipo di coppia: fra atomi polari N/O si accetta un legame H
    da 2.6 A, per tutte le altre coppie si chiede 3.2 A. Ritorna il margine peggiore e la coppia."""
    peggio = (np.inf, None)
    for n, p in sc:
        for m, q, tag in ambiente:
            d = np.linalg.norm(p - q)
            soglia = 2.6 if n[0] in "NO" and m[0] in "NO" else 3.2
            if d - soglia < peggio[0]:
                peggio = (d - soglia, (n, tag, round(d, 2)))
    return peggio

# rotameri canonici di chi1 in ordine di frequenza (librerie backbone-indipendenti)
ROTAMERI = {"VAL": [175, -60, 60], "THR": [60, -60, 180], "GLN": [-65, 180, 60]}
GAMMA = {"VAL": "CG1", "THR": "OG1", "GLN": "CG"}

nmr, xr = leggi_pdb(NMR), leggi_gro(XRAY)
residui = sorted({a["resi"] for a in nmr})
out = []
for k in residui:
    res = [a for a in nmr if a["resi"] == k]
    nome = res[0]["resn"]
    if k == 86:
        assert nome == "ASN"
        nuovi = []
        for a in res:
            a = dict(a, resn="ASP")
            if a["name"] == "ND2":
                a["name"], a["elem"] = "OD2", "O"
            nuovi.append(a)
        res = sorted(nuovi, key=lambda a: ORDINE["ASP"].index(a["name"]))
    elif k in TRONCA:
        assert nome == TRONCA[k], (k, nome)
        res = [dict(a, resn="ALA") for a in res if a["name"] in ORDINE["ALA"]]
        res = sorted(res, key=lambda a: ORDINE["ALA"].index(a["name"]))
    elif k in INNESTA:
        da, a_ = INNESTA[k]
        assert nome == da, (k, nome)
        rx = {a["name"]: a["xyz"] for a in xr if a["resi"] == k and not a["name"].startswith("H")}
        assert set(rx) == set(ORDINE[a_]), (k, sorted(rx))
        bb = {a["name"]: a for a in res if a["name"] in ("N", "CA", "C", "O")}
        R, t = kabsch(np.array([rx[n] for n in ("N", "CA", "C")]), np.array([bb[n]["xyz"] for n in ("N", "CA", "C")]))
        sc = [(n, rx[n] @ R + t) for n in ORDINE[a_][4:]]
        # ambiente: atomi pesanti degli altri residui, esclusi C(k-1) e N(k+1) legati al residuo
        amb = [(b["name"], b["xyz"], f"{b['resi']}{b['name']}") for b in nmr if b["resi"] != k
               and not (b["resi"] == k - 1 and b["name"] == "C") and not (b["resi"] == k + 1 and b["name"] == "N")]
        ca, cb = bb["CA"]["xyz"], dict(sc)["CB"]
        chi_x = diedro(bb["N"]["xyz"], ca, cb, dict(sc)[GAMMA[a_]])
        def con_chi1(chi):
            return [(n, p if n == "CB" else ruota(p, ca, cb - ca, chi - chi_x)) for n, p in sc]
        scelta, prove = None, []
        for nome_r, chi in [("X-ray", chi_x)] + [(f"canonico {c}", c) for c in ROTAMERI[a_]]:
            margine, coppia = distanza_minima(con_chi1(chi), amb)
            prove.append((margine, nome_r, chi))
            print(f"   {k}{a_} chi1 {chi:7.1f} ({nome_r}): margine {margine:+.2f} A, coppia peggiore {coppia}")
            if scelta is None and margine >= 0:
                scelta = (nome_r, chi)
        if scelta is None:
            # nessun rotamero pulito: si tiene il rotamero (X-ray o canonico) col contatto piu' lieve,
            # purche' entro 0.5 A dalla soglia, che la minimizzazione risolve; un chi1 eclissato no
            margine, nome_r, chi = max(prove)
            if margine < -0.5:
                raise SystemExit(f"{k}{a_}: contatto troppo stretto anche col rotamero migliore")
            scelta = (f"{nome_r}, contatto lieve {margine:+.2f} A", chi)
        print(f"   {k}{a_}: scelto chi1 {scelta[1]:.1f} ({scelta[0]})")
        res = [dict(bb[n], resn=a_) for n in ("N", "CA", "C", "O")]
        for n, p in con_chi1(scelta[1]):
            res.append({"name": n, "resn": a_, "chain": bb["N"]["chain"], "resi": k, "xyz": p, "elem": n[0]})
    out.extend(res)

with open(OUT, "w") as f:
    f.write("REMARK   medoide NMR 9HW9 (modello 6, ossidato) con la sequenza del costrutto X-ray\n")
    f.write("REMARK   N86D T106A S128A T143A L145V K160A A167Q A182T; 145/167/182 innestati dall'X-ray\n")
    for i, a in enumerate(out, 1):
        nm = a["name"] if len(a["name"]) == 4 else " " + a["name"]
        x, y, z = a["xyz"]
        f.write(f"ATOM  {i:5d} {nm:<4s} {a['resn']:>3s} {a['chain']}{a['resi']:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a['elem']:>2s}\n")
    f.write("TER\nEND\n")
print(f"scritti {len(out)} atomi in {OUT}")
