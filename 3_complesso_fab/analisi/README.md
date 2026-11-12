# Analisi – complesso MrkA–Fab e libero

Gli script applicano il documento [MrkA–Fab: criteri di analisi](https://claude.ai/artifact/DPDeZEeE17qTTZxD3ukXdX)
(chiuso il 2 nov 2026, prima di aprire le traiettorie). Le scelte operative che il documento non fissava sono
marcate `NON NEL DOC` nel codice: le ha decise il revisore il 5 nov 2026, prima dei dati, e sono nel registro del
documento.

```
09_estrai.sh ─► 10_serie.py ─► 11_controlli.py ─► 12_scopo1_contatti.py
  (gmx: solo     (serie per  │   (validità delle     13_scopo2_epitopo.py
   proteina)      fotogramma)│    repliche)           14_scopo3_rmsf.py ─► 15_figure.py
                             ├─► 16_verifica_immagini.py (gmx mindist -pi su alcuni fotogrammi)
                             └─► 17_cerniera.py (diagnostica descrittiva Fv-core, solo complesso)
12, 13, 14, 17 (json finali) ─► 18_figura_articolo.py (figure per l'articolo, solo disegno)
```

| Script | Cosa fa | Uscita |
|---|---|---|
| `08_verifica_numerazione.py` | Numerazione e sequenza di ogni catena nei `.gro` = `mrka-fab_completo.pdb` (prerequisito). | a schermo |
| `09_estrai.sh` | Con `gmx trjconv` estrae il gruppo `Protein_Other` da ogni `run<r>.xtc`, dall'inizio della produzione (tinit del tpr) fino a `FINO` ns; salta le repliche che non ci sono arrivate. Non tocca i file di produzione e funziona anche mentre mdrun allunga l'xtc. | `traiettorie/` |
| `10_serie.py` | Ricompone le catene spezzate dalle PBC e riassembla A, L e H; per fotogramma calcola RMSD dei 43 Cα, RMSD del backbone dell'epitopo, Cα della catena A dopo il fit, distanza dalle immagini (come `gmx mindist -pi`), contatti A–Fab entro 5 Å e distanza minima A–Fab, distanza fra i centri dell'Fv e dei 43 Cα del core, RMSD del Fab e dell'Fv, omega delle Pro cis; solo nel complesso anche la misura secondaria dello scopo 2 (RMSD dell'epitopo dopo fit sull'epitopo stesso). | `cache/serie_<sistema>_run<r>.npz` |
| `11_controlli.py` | Validità delle repliche con scarto 50 e 100 ns (contatti persi = occupanza < 0.5; centri Fv–core > 5 Å in > 50% dei fotogrammi); controllo del riassemblaggio (A–Fab < 4 Å in ogni fotogramma); durate confermative T3 (complesso) e T2 (libero); stato INTERMEDIA/FINALE per scopo. | `risultati/controlli.txt`, `.json` |
| `12_scopo1_contatti.py` | Scopo 1: contatti persistenti, persi e nuovi a 4.0 Å (3.5 e 5.0 di sensibilità). | `scopo1_contatti.txt`, `.json`, `.csv` |
| `13_scopo2_epitopo.py` | Scopo 2: soglia P95 dal complesso, frazione entro soglia, uscite e rientri (1 ns) nel libero. Misura primaria (fit sul core, doc) e secondaria (fit sull'epitopo stesso, aggiunta il 6 nov 2026 dopo l'intermedia, approvata da francesco e dal revisore) affiancate, con "esiti concordi / discordi" nell'analisi finale. La secondaria del libero si calcola dagli xtc solo con lo stato FINALE per lunghezze (5 × 500 ns, 10 × 1 µs; con il solo `--finale` di 11, 13 si ferma); controllo: una replica del complesso ricalcolata dall'xtc = cache di 10 entro 1e-6 Å. | `scopo2_epitopo.txt`, `.json` |
| `14_scopo3_rmsf.py` | Scopo 3: RMSF dei Cα e medie di segmento 87-101 e 138-147; Welch + Holm solo nell'analisi finale. | `scopo3_rmsf.txt`, `.json`, `scopo3_rmsf_per_residuo.csv` |
| `15_figure.py` | Figure descrittive (serve matplotlib). | `fig_*.png` |
| `16_verifica_immagini.py` | Controllo del revisore: su 8 fotogrammi per replica `gmx trjconv -pbc cluster` (Protein_Other) dalla traiettoria di produzione, poi `gmx mindist -pi`, contro la distanza di 10 ricalcolata senza tetto; tolleranza 0.005 Å. Serve GROMACS. | `verifica_immagini.txt`, `.json` |
| `17_cerniera.py` | Diagnostica descrittiva, NON criterio (6 nov 2026, dopo l'intermedia a 250 ns; nota nel registro del doc): solo complesso, ogni 1 ns, RMSD del backbone dell'epitopo dopo fit sul core (= cache, controllo), sull'Fv e sull'epitopo stesso; Cα 90-93 dopo fit sul core e sull'Fv; rotazione e spostamento dell'Fv rispetto al core; scostamento per residuo dei Cα dell'epitopo. Nel libero il fit sull'epitopo non si calcola prima dei dati finali (revisore, 6 nov). Dall'11 nov 2026 scrive lo stato (INTERMEDIA o FINALE) dalle lunghezze delle repliche del complesso, non più un'etichetta fissa. | `diagnostica_cerniera.txt`, `.json`, `fig_diagnostica_cerniera.png` |
| `18_figura_articolo.py` | Figure per l'articolo, in inglese, SOLO DISEGNO (condizioni del revisore, 11 nov 2026): non ricalcola e non riaggrega niente; ogni numero disegnato è una foglia dei json finali di 12, 13, 14 e 17, le cui sha256 stanno nel txt. Si ferma se uno scopo non è FINALE, se manca la secondaria del libero o se un esito dei contatti cambia fra scarto 50 e 100 ns. | `fig_articolo_md.png/.pdf` (principale), `_S1` (scopo 2 nel tempo, blocchi di 50 ns), `_S2` (rotazione dell'Fv), `fig_articolo_md.txt` |
| `rilancia_finale.sh`, `SHA256SUMS` | Analisi finale con una sola versione: controlla le sha256 degli script contro `SHA256SUMS` (si ferma se una è diversa), mette da parte `analisi/risultati` e rilancia 10-18 nell'ordine. `SHA256SUMS` va rigenerato a ogni modifica degli script (`sha256sum 0*.py 0*.sh 1*.py comune_fab.py > SHA256SUMS`). | – |
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
python3 analisi/16_verifica_immagini.py
python3 analisi/17_cerniera.py
python3 analisi/18_figura_articolo.py      # solo con i dati finali
```

Per l'analisi finale, dopo aver copiato in `analisi/` tutti i file di questa cartella, basta un comando, che
controlla le versioni e rilancia 10-18:

```
bash analisi/rilancia_finale.sh
```

- Pacchetti: Python 3 con numpy, scipy (≥ 1.7), MDAnalysis (≥ 2.0); matplotlib solo per le figure.
- Per allungare l'analisi (per esempio a 500 ns) si rilancia tutto con `FINO=500`: 10 usa per ogni replica la
  traiettoria estratta più lunga.
- Analisi INTERMEDIA finché non ci sono tutte le repliche del disegno alla durata massima (5 × 500 ns, 10 × 1 µs):
  solo numeri per replica, senza esiti (persistente, popolata, ...) né regole di conteggio (ceil 0.8, 4/5, 8/10,
  "non concludibile"), senza test; nessuna frase del paper da questi numeri; nessuna replica si rilancia o si
  sostituisce per quello che si vede (condizioni del revisore, 5 nov 2026). `--finale` in 11 dichiara finale un
  disegno ridotto, solo se lo decide francesco prima di guardare le traiettorie.
- 10 rifà una serie se è cambiato uno qualsiasi di xtc, `proteina_<sistema>.gro`, pdb di riferimento, `10_serie.py`
  o `comune_fab.py` (sha256 confrontati con la provenienza della serie, chiave per chiave).
- Una sola versione per tutte le serie (revisore, 7 nov 2026): 11-17 si fermano se una cache non viene dal
  `10_serie.py` e dal `comune_fab.py` attuali; basta rilanciare 10, che rifà solo quelle.
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
- 5 nov, dopo il verdetto del revisore: stesse traiettorie sintetiche (replica "da decidere" ora valida, replica
  non valida per contatti persi e centri Fv–core); 16 con GROMACS 2023 su fotogrammi spezzati in una cella stretta
  (distanze 10-52 Å): differenza massima 0.0047 Å da `gmx mindist -pi`. Con l'uscita di `-pbc cluster` in xtc la
  differenza arrivava a 0.013 Å (riarrotondamento a 0.001 nm degli atomi spostati): per questo 16 scrive un trr.
- 6 nov, misura secondaria dello scopo 2: 10 rifatto sui 4 complessi veri (campi vecchi uguali entro 2e-14 Å, primaria
  identica: P95 3.904 Å); secondaria del complesso = MDAnalysis `rms.rmsd(superposition=True)` entro 1e-14 Å; controllo
  di 13 (ricalcolo dall'xtc = cache) 0 Å. Sintetici: intermedia senza secondaria del libero; `--finale` senza lunghezze
  -> 13 si ferma; cache senza secondaria -> 13 si ferma; FINALE per lunghezze (disegno ridotto solo nella copia di
  collaudo) -> secondaria del libero calcolata, esiti affiancati con "concordi/discordi", nessun file della secondaria
  del libero su disco.
- 7 nov, una sola versione delle serie: sui sintetici con le cache del giorno prima 11 si ferma (15 su 15 di versione
  diversa); dopo 10 l'intera catena 11-15 e 17 gira; con una sola cache vecchia rimessa (libero run3) 11 e 13 si
  fermano e 10 rifà solo quella.
- 7 nov, ingressi di 10: con xtc, gro, pdb e script invariati 10 salta la serie; cambiando il gro del complesso
  rifà solo il complesso; cambiando il pdb rifà tutte e due; rilanciato subito dopo, salta.
