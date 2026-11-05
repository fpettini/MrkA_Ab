# Analisi – complesso MrkA–Fab e libero

Gli script applicano il documento [MrkA–Fab: criteri di analisi](https://claude.ai/artifact/DPDeZEeE17qTTZxD3ukXdX)
(chiuso il 2 nov 2026, prima di aprire le traiettorie). Le poche scelte che il documento non fissa sono marcate
`NON NEL DOC` nel codice e nei txt, e vanno confermate dal revisore.

```
09_estrai.sh ─► 10_serie.py ─► 11_controlli.py ─► 12_scopo1_contatti.py
  (gmx: solo     (serie per      (validità delle     13_scopo2_epitopo.py
   proteina)      fotogramma)     repliche)           14_scopo3_rmsf.py ─► 15_figure.py
```

| Script | Cosa fa | Uscita |
|---|---|---|
| `08_verifica_numerazione.py` | Numerazione e sequenza di ogni catena nei `.gro` = `mrka-fab_completo.pdb` (prerequisito). | a schermo |
| `09_estrai.sh` | Con `gmx trjconv` estrae il gruppo `Protein_Other` da ogni `run<r>.xtc`, dall'inizio della produzione (tinit del tpr) fino a `FINO` ns; salta le repliche che non ci sono arrivate. Non tocca i file di produzione e funziona anche mentre mdrun allunga l'xtc. | `traiettorie/` |
| `10_serie.py` | Ricompone le catene spezzate dalle PBC e riassembla A, L e H; per fotogramma calcola RMSD dei 43 Cα, RMSD del backbone dell'epitopo, Cα della catena A dopo il fit, distanza dalle immagini (come `gmx mindist -pi`), contatti A–Fab entro 5 Å, centri, RMSD del Fab e dell'Fv, omega delle Pro cis. | `cache/serie_<sistema>_run<r>.npz` |
| `11_controlli.py` | Validità delle repliche con scarto 50 e 100 ns; durate confermative T3 (complesso) e T2 (libero); stato INTERMEDIA/FINALE per scopo. | `risultati/controlli.txt`, `.json` |
| `12_scopo1_contatti.py` | Scopo 1: contatti persistenti, persi e nuovi a 4.0 Å (3.5 e 5.0 di sensibilità). | `scopo1_contatti.txt`, `.json`, `.csv` |
| `13_scopo2_epitopo.py` | Scopo 2: soglia P95 dal complesso, frazione entro soglia, uscite e rientri (1 ns) nel libero. | `scopo2_epitopo.txt`, `.json` |
| `14_scopo3_rmsf.py` | Scopo 3: RMSF dei Cα e medie di segmento 87-101 e 138-147; Welch + Holm solo nell'analisi finale. | `scopo3_rmsf.txt`, `.json`, `scopo3_rmsf_per_residuo.csv` |
| `15_figure.py` | Figure descrittive (serve matplotlib). | `fig_*.png` |
| `comune_fab.py` | Definizioni dal documento (residui, soglie, scarti), lettura del gro, PBC, fit, provenienza. | – |

## Comandi (sulla workstation, dalla cartella `complesso_fab`)

```
FINO=250 bash analisi/09_estrai.sh
python3 analisi/10_serie.py
python3 analisi/11_controlli.py
python3 analisi/12_scopo1_contatti.py
python3 analisi/13_scopo2_epitopo.py
python3 analisi/14_scopo3_rmsf.py
python3 analisi/15_figure.py
```

- Pacchetti: Python 3 con numpy, scipy (≥ 1.7), MDAnalysis (≥ 2.0); matplotlib solo per le figure.
- Per allungare l'analisi (per esempio a 500 ns) si rilancia tutto con `FINO=500`: 10 usa per ogni replica la
  traiettoria estratta più lunga.
- Analisi INTERMEDIA finché non ci sono tutte le repliche del disegno alla durata massima (5 × 500 ns, 10 × 1 µs):
  i numeri sono descrittivi e i test dello scopo 3 non si fanno. `--finale` in 11 dichiara finale un disegno
  ridotto, solo se lo decide francesco prima di guardare le traiettorie.
- Ogni txt ha in testa la provenienza (sha256 di script, `comune_fab.py` e input, versioni, argomenti); 09 la
  scrive in `traiettorie/provenienza_09.txt`.
- Mai il tpr come topologia in MDAnalysis (rinumera i residui): 10 usa `proteina_<sistema>.gro` e controlla
  all'avvio A87 CYX, A94 LYS, L214 CCYX, H220 CYX e la sequenza di ogni catena contro il pdb.

## Collaudo (5 nov 2026, prima dei dati veri)

- 09 con GROMACS 2023 su due MD brevi vere della catena A (tinit 0 e 20000 ps): estrazione, salto delle repliche
  corte, controllo dell'ultimo fotogramma.
- Distanza dalle immagini di 10 = `gmx mindist -pi` sulla stessa traiettoria entro 0.005 Å.
- 10-15 su 15 traiettorie sintetiche da 250 ns (gro con i nomi della produzione, molecole spezzate dalle PBC in un
  dodecaedro) con esito noto: replica non valida per contatti persi e centri, replica "da decidere", replica non
  valida per le immagini, uscite e rientri, RMSF diverso solo in 87-101.
- Verifica indipendente con MDAnalysis: contatti identici, RMSD entro 1e-5 Å, RMSF entro 5e-4 Å.
