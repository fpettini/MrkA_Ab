# 3 – Complesso MrkA–Fab e MrkA libera

Struttura X-ray del complesso con il Fab completato (catena intera, S–S intercatena L214–H220), MrkA libera
dalla stessa struttura. amber19sb.ff + OPC (GROMACS 2026), 298 K, KCl 0.2 M, 7 acque a ponte + 12 interne.
Produzione a 2 fs (HMR 1): complesso 5 × 500 ns, libero 10 × 1 µs, a tratti alternati (06, `TURNI=alternati`).

```
preparazione/  01_ricostruisci_fab ─► 01b_acque_interne ─► 02_prepara_sistemi (+02a_aiuto)
               ─► controlli_ff/04, 04b, 04c (GROMACS = OpenMM)  ─► densita/ (OPC+KCl vs OpenMM)
produzione/    03_test_velocita ─► 05_equil_s1_s4 (+05a) ─► 07_test_deriva (2/3/4 fs) ─► 06_repliche
analisi/       08_verifica_numerazione ─► 09_estrai (gmx) ─► 10_serie ─► 11_controlli
               ─► 12 scopo 1 (contatti), 13 scopo 2 (epitopo), 14 scopo 3 (RMSF) ─► 15_figure; 16_verifica_immagini (gmx mindist -pi)
               (criteri fissati il 2 nov 2026, vedi analisi/README.md)
```
