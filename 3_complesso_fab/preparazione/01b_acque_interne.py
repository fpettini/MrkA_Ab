#!/usr/bin/env python3
"""Seleziona le acque cristallografiche sepolte nella sola MrkA, da mettere identiche in complesso e libero.

Uso: python3 01b_acque_interne.py <originale/mrka-Fab.pdb> <mrka-fab_completo.pdb> <acque_interne.txt>

Criterio (fissato il 31 ott 2026, su richiesta del revisore: acque interne di MrkA in tutti e due i sistemi o in nessuno):
 - acqua entro 4.0 A da un atomo pesante della catena A;
 - sepolta nella sola MrkA: frazione accessibile dell'ossigeno < 5% (Shrake-Rupley, 400 punti, r acqua 1.4 + sonda 1.4,
   atomi pesanti 1.8 A), calcolata con la sola catena A del modello completo (loop 116-117 compreso);
 - nessun atomo polare (N, O, S) del modello completo (A, B, C) a meno di 2.5 A e nessun carbonio a meno di 3.0 A;
 - almeno un contatto polare fra 2.5 e 3.3 A;
 - escluse le 7 acque a ponte MrkA-Fab (gia' nel complesso); distanza acqua-acqua >= 2.4 A fra tutte quelle tenute (7 comprese).
Soglie 2.5/3.0/3.3/2.4 dal revisore (31 ott 2026). Sensibilita': quante sarebbero con accessibilita' < 10%.
Altloc: si tiene la A. Uscita: una riga per acqua con id originale e coordinate in nm.
"""
import sys
import numpy as np
import gemmi

SETTE = {("E", 11), ("E", 43), ("E", 129), ("E", 146), ("S", 501), ("S", 2103), ("S", 4277)}


def main(orig, completo, out):
    o = gemmi.read_structure(orig); o.remove_alternative_conformations()
    c = gemmi.read_structure(completo)[0]
    pa = np.array([a.pos.tolist() for r in c["A"] for a in r if a.element.name != "H"])
    pall = np.array([a.pos.tolist() for ch in "ABC" for r in c[ch] for a in r if a.element.name != "H"])
    polare = np.array([a.element.name in ("N", "O", "S") for ch in "ABC" for r in c[ch] for a in r if a.element.name != "H"])
    ponte = [np.array(r[0].pos.tolist()) for ch in c for r in ch if r.name == "HOH"]
    i = np.arange(400) + 0.5; phi = np.arccos(1 - 2 * i / 400); th = np.pi * (1 + 5 ** 0.5) * i
    sfera = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]

    def accessibile(x):
        punti = x + 2.8 * sfera
        vic = pa[np.linalg.norm(pa - x, axis=1) < 6.1]
        return (np.linalg.norm(punti[:, None] - vic[None], axis=2) > 3.2).all(1).mean()

    sepolte, sens10, scelte, escluse = 0, 0, [], []
    for ch in o[0]:
        for r in ch:
            if r.name != "HOH" or (ch.name, r.seqid.num) in SETTE:
                continue
            x = np.array(r[0].pos.tolist())
            if np.linalg.norm(pa - x, axis=1).min() > 4.0:
                continue
            acc = accessibile(x)
            sens10 += acc < 0.10
            if acc >= 0.05:
                continue
            sepolte += 1
            d = np.linalg.norm(pall - x, axis=1)
            dp, dc = d[polare].min(), d[~polare].min()
            ok = dp >= 2.5 and dc >= 3.0 and ((d[polare] >= 2.5) & (d[polare] <= 3.3)).any()
            (scelte if ok else escluse).append((ch.name, r.seqid.num, x, r[0].b_iso, dp, dc))
    # acqua-acqua >= 2.4 A (con le 7 a ponte e fra le tenute)
    for w in scelte:
        vic = [v for v in ponte + [z[2] for z in scelte if z is not w] if np.linalg.norm(v - w[2]) < 2.4]
        assert not vic, f"{w[0]}{w[1]}: acqua a meno di 2.4 A"
    scelte = [(ch, n, x / 10, b, dp, dc) for ch, n, x, b, dp, dc in scelte]
    with open(out, "w") as f:
        for ch, n, x, b, dp, dc in scelte:
            f.write(f"{x[0]:.4f} {x[1]:.4f} {x[2]:.4f}  # {ch}{n} B {b:.1f} polare {dp:.2f} C {dc:.2f}\n")
    print(f"sepolte nella sola MrkA (<5%): {sepolte} (sensibilita' <10%: {sens10}); tenute: {len(scelte)}")
    print("tenute:", " ".join(f"{ch}{n}" for ch, n, *_ in scelte))
    print("escluse:", " ".join(f"{ch}{n}(pol {dp:.2f}, C {dc:.2f})" for ch, n, _, _, dp, dc in escluse))


if __name__ == "__main__":
    main(*sys.argv[1:4])
