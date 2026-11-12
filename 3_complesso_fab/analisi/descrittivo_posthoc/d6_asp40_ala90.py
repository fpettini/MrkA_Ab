#!/usr/bin/env python3
"""Descrittivo POST HOC (12 nov 2026): il legame idrogeno del cristallo fra un ossigeno carbossilico di Asp40 e l'N di
Ala90 nelle simulazioni. Non e' un criterio e non c'e' nessun test. Legge le serie per fotogramma d1_<sistema>_run<r>.npz
(distanze D40OD1_A90N e D40OD2_A90N, in A) prodotte da d1_serie.py dalle traiettorie estratte.
Uso: python3 d6_asp40_ala90.py <cartella descrittivo> <pdb del modello completato>"""
import sys, numpy as np
d, pdb = sys.argv[1], sys.argv[2]
SOGLIA = 3.5
at = {}
for riga in open(pdb):
    if riga.startswith('ATOM') and riga[21] == 'A':
        at[(int(riga[22:26]), riga[12:16].strip())] = np.array([float(riga[30:38]), float(riga[38:46]), float(riga[46:54])])
# numerazione del pdb: catena A numerata come nel cristallo; Asp40 e Ala90
def dist(a, b): return float(np.linalg.norm(at[a] - at[b]))
print('POST HOC 12 nov 2026, descrittivo, nessun test. Soglia N-O %.1f A.' % SOGLIA)
print('cristallo (modello completato): Asp40 OD1 - Ala90 N %.2f A; Asp40 OD2 - Ala90 N %.2f A' % (dist((40, 'OD1'), (90, 'N')), dist((40, 'OD2'), (90, 'N'))))
print()
print('replica        finestra ns   fraz. (OD1 o OD2) <= 3.5   fraz. OD1 <= 3.5   distanza minima media   mediana   primo fotogramma')
def riga(sis, r, t0, t1):
    z = np.load(f'{d}/d1_{sis}_run{r}.npz')
    t = z['t_ns']; m = (t >= t0) & (t <= t1)
    d1, d2 = z['D40OD1_A90N'], z['D40OD2_A90N']; dm = np.minimum(d1, d2)
    print(f'{sis+" run"+str(r):<14} {t0:>4.0f}-{t1:<6.0f}   {(dm[m] <= SOGLIA).mean():>10.3f}               {(d1[m] <= SOGLIA).mean():>8.3f}           {dm[m].mean():>8.2f}            {np.median(dm[m]):>6.2f}    {dm[0]:>6.2f}')
    return (dm[m] <= SOGLIA).mean()
fc = [riga('complesso', r, 50, 500) for r in range(1, 6)]
fl = [riga('libero', r, 50, 1000) for r in range(1, 11)]
fl5 = [riga('libero', r, 50, 500) for r in range(1, 11)]
print()
print('complesso 50-500 ns: frazione per replica min %.3f max %.3f' % (min(fc), max(fc)))
print('libero 50-1000 ns:   frazione per replica min %.3f max %.3f' % (min(fl), max(fl)))
print('libero 50-500 ns:    frazione per replica min %.3f max %.3f' % (min(fl5), max(fl5)))
