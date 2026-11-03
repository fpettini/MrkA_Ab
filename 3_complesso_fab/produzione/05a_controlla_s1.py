#!/usr/bin/env python3
"""Controllo dopo equil_s1 chiesto dal revisore (1 nov 2026). A fine min_c3 la forza massima e' sul CG di Asp40 di MrkA
(2.6e3 kJ/mol/nm nel complesso, 2.4e3 nel vecchio complesso, 2.0e3 nel libero). Per gli atomi pesanti di Asp40 stampa i tre
contatti pesanti piu' vicini (esclusi i residui 39-41) prima (min_c3) e dopo s1, e la distanza minima. Gli avvisi LINCS li conta 05.

Uso (dalla cartella del sistema): python3 05a_controlla_s1.py min_c3.gro equil_s1.gro      (solo libreria standard)
"""
import sys


def leggi(p):
    r = open(p).read().splitlines()
    n = int(r[1])
    at = [(int(l[0:5]), l[5:10].strip(), l[10:15].strip(), (float(l[20:28]), float(l[28:36]), float(l[36:44])))
          for l in r[2:2 + n]]
    b = [float(x) for x in r[2 + n].split()] + [0.0] * 6
    # formato gro: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
    return at, ((b[0], b[3], b[4]), (b[5], b[1], b[6]), (b[7], b[8], b[2]))


def distanza(p, q, scatola):
    """immagine minima per scatola triclina (come GROMACS: prima c, poi b, poi a), in A"""
    d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]]
    for v, k in ((scatola[2], 2), (scatola[1], 1), (scatola[0], 0)):
        if v[k] > 0:
            s = round(d[k] / v[k])
            if s:
                d = [d[0] - s * v[0], d[1] - s * v[1], d[2] - s * v[2]]
    return 10 * (d[0] ** 2 + d[1] ** 2 + d[2] ** 2) ** 0.5


def pesante(nome):
    return not nome.startswith("H") and nome != "MW"


def contatti(gro):
    at, scatola = leggi(gro)
    # Asp40 di MrkA = primo residuo 40 del file (MrkA e' la prima catena); esclusi i residui contigui 39-41
    i40 = next(i for i, a in enumerate(at) if a[0] == 40 and a[1].endswith("ASP"))
    ini = i40
    while ini > 0 and at[ini - 1][0] in (39, 40):
        ini -= 1
    fin = i40
    while fin + 1 < len(at) and at[fin + 1][0] in (40, 41):
        fin += 1
    asp = [i for i in range(ini, fin + 1) if at[i][0] == 40 and pesante(at[i][2])]
    altri = [i for i in range(len(at)) if (i < ini or i > fin) and pesante(at[i][2])]
    out = {}
    for i in asp:
        d = sorted((distanza(at[i][3], at[j][3], scatola), j) for j in altri)[:3]
        out[at[i][2]] = [(f"{at[j][1]}{at[j][0]}:{at[j][2]}", x) for x, j in d]
    return out


def main(prima, dopo):
    a, b = contatti(prima), contatti(dopo)
    print(f"Asp40 di MrkA, atomi pesanti: 3 contatti pesanti piu' vicini (residui 39-41 esclusi), in A")
    print(f"{'atomo':6s} {prima:38s} {dopo}")
    for nome in a:
        fa = ", ".join(f"{n} {x:.2f}" for n, x in a[nome]); fb = ", ".join(f"{n} {x:.2f}" for n, x in b[nome])
        print(f"{nome:6s} {fa:38s} {fb}")
    ma = min(x for v in a.values() for _, x in v); mb = min(x for v in b.values() for _, x in v)
    print(f"distanza minima pesante-pesante: {ma:.2f} A ({prima}) -> {mb:.2f} A ({dopo})"
          + ("   SOTTO 2.5 A: da guardare prima di s2" if mb < 2.5 else ""))


if __name__ == "__main__":
    main(*sys.argv[1:3])
