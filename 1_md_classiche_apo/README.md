# 1 – MD classiche di MrkA libera (X-ray vs NMR)

Disegno: X-ray ossidato, NMR (9HW9, modello 6 = medoide) ossidato e ridotto; 3 repliche indipendenti per sistema
(prima 200 ns, poi estese a 500 ns). amber14sb + TIP3P, NaCl 0.15 M, dodecaedro `-d 1.3`, 310 K, pH 7.4,
terminali liberi e carichi (NSER 42, CGLN 202), S–S Cys43–Cys87 esplicito nei sistemi ossidati.

```
preparazione/   protonazione (PDB2PQR/PROPKA, pH 7.4) ─► pdb2gmx ─► solvatazione/ioni ─► min_c1-c3
produzione/     equil_s1-s4 ─► equil_s5 (20 ns, velocità nuove per replica) ─► produzione
analisi/        01_estrai ─► 02_analisi ─► 03_figure ─► 04_confronto_ristretto
                bfactor/05_rmsf_vs_bfactor , bfactor/06_contatti_cristallo
```

- [`preparazione/`](preparazione/) – vedi README (contiene solo l'input APBS: gli script di setup sono sul PC).
- `produzione/` – **non ancora caricata**: lo script `eq.sh` e gli mdp delle MD classiche sono solo sul PC.
- [`analisi/`](analisi/) – convergenza, confronto ristretto, RMSF vs B-factor, contatti cristallini.
