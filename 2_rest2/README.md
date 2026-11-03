# 2 – REST2 su MrkA (X-ray e NMR con la sequenza X-ray)

16 repliche, temperature efficaci 310.15–480 K in progressione geometrica, atomi "caldi" = loop
59-77, 87-101, 112-118, 124-133, 138-147, 159-169, 180-192. GROMACS 2024.6 + PLUMED 2.10.2 (`-hrex`),
scambio ogni 1 ps, `-notunepme`. Produzione 200 ns, poi estesa a 400 ns. REST3 preparato e collaudato ma non lanciato.

```
preparazione/  00_compila ─► 01_muta_nmr ─► 02_controlla ─► (pdb2gmx, equil come le MD classiche; grompp -pp)
               ─► 03_marca_caldi ─► 04_prepara_rest2  [REST3: 07_kappa_rest3 + 08_prepara_rest3]
produzione/    05_lancia_rest2 ─► 06_controlla_scambi (in corsa) ─► 09_allunga_rest2 ─► 05 (ripresa)
analisi/       01_estrai ─► 02_controlli ─► 03_popolazioni
```

Il file `produzione/COMANDI.md` raccoglie la sequenza completa dei comandi usati sulla workstation.
