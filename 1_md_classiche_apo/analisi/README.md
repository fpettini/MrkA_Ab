# Analisi – MD classiche apo

Workflow (dalla cartella del progetto `mrkA`; traiettorie attese in `traj/<sistema>_run<r>_pr.xtc`):

```
01_estrai.py ─► cache/*.npz (sola proteina)
      │
      ├─► 02_analisi.py ─► cache/*_serie.npz + numeri di RIEPILOGO_convergenza ─► 03_figure.py (fig. 1-5)
      ├─► 04_confronto_ristretto.py ─► confronto_ristretto.txt/.json
      └─► bfactor/05_rmsf_vs_bfactor.py (solo X-ray)
bfactor/06_contatti_cristallo.py (solo dal pdb del cristallo, indipendente dalle MD)
disulfide_trajectory.py (distanza SG–SG lungo qualunque xtc)
```

| Script | Scopo |
|---|---|
| `01_estrai.py` | Estrae dalle xtc i soli atomi della proteina di xray, nmr_ox, nmr_red (3 repliche) e li salva in `.npz` (in parallelo). |
| `02_analisi.py` | Verifica di convergenza: RMSD Cα 44–199 e del core β, Rg, RMSF per replica, RMSD incrociato fra repliche, DSSP, tempi di equilibrazione (pymbar). Stampa i numeri del riepilogo. |
| `03_figure.py` | Figure di convergenza (serie temporali, RMSF per replica, RMSD incrociato, struttura secondaria, convergenza vs durata). |
| `04_confronto_ristretto.py` | Confronto fra sistemi sulle sole grandezze convergenti (core β, distanza S–S, Rg), media 50–500 ns per replica, unità statistica = replica. |
| `disulfide_trajectory.py` | Monitora i ponti S–S dichiarati nella topologia lungo una o più traiettorie (immagine minima, legge solo gli SG); uscita a schermo o CSV. |
| `bfactor/05_rmsf_vs_bfactor.py` | RMSF dei Cα X-ray (500 ns) contro i B-factor del cristallo: Spearman per residuo e rango percentile per segmento; richiede `--shift -35` (pdb = wt − 35). Solo qualitativo. |
| `bfactor/06_contatti_cristallo.py` | Contatti della catena A nel cristallo (Fab nella stessa unità asimmetrica + copie di simmetria da CRYST1, C 1 2 1) per i tratti 87–101 e 138–147 contro il resto; soglie 4.0 Å (sensibilità 3.5/5.0). |
| `bfactor/provenienza.py` | Modulo di provenienza (copia di quello di `2_rest2/analisi/`), usato da 05 e 06. |

Convenzioni: residui 42–43 e 200–202 esclusi dalle metriche per residuo (cariche ai siti di troncamento);
116–117 ricostruiti; allineamento sul dominio pilin.
