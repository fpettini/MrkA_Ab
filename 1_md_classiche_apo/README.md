# 1 – MD classiche di MrkA libera (X-ray vs NMR)

Disegno: X-ray ossidato, NMR (9HW9, modello 6 = medoide) ossidato e ridotto; 3 repliche indipendenti per sistema
(prima 200 ns, poi estese a 500 ns). amber14sb + TIP3P, NaCl 0.15 M, dodecaedro `-d 1.3`, 310 K, pH 7.4,
terminali liberi e carichi (NSER 42, CGLN 202), S–S Cys43–Cys87 esplicito nei sistemi ossidati.

```
preparazione/   PDB2PQR/PROPKA (pH 7.4) ─► eq_<sistema>.sh: pdb2gmx ─► box, solvatazione, ioni ─► min_c1-c3 (doppia precisione)
                ─► equil_s1-s4 ─► per replica equil_s5 (20 ns, velocità nuove) ─► produzione 200 ns
produzione/     eq_estensione_500ns.sh: estensione delle repliche a 500 ns
analisi/        01_estrai ─► 02_analisi ─► 03_figure ─► 04_confronto_ristretto ; disulfide_trajectory
                bfactor/05_rmsf_vs_bfactor , bfactor/06_contatti_cristallo
```

- [`preparazione/`](preparazione/) – script di setup ed equilibrazione per sistema, mdp di minimizzazione e ioni.
- [`produzione/`](produzione/) – mdp di equilibrazione e produzione, estensione a 500 ns.
- [`analisi/`](analisi/) – convergenza, confronto ristretto, ponte S–S, RMSF vs B-factor, contatti cristallini.

Nota: alcuni commenti negli mdp e negli `eq.sh` ("zuccheri", "NLN/OLT", "844k atomi") vengono dal template
di un altro progetto con glicani e non descrivono questo sistema; i valori sono quelli giusti.
