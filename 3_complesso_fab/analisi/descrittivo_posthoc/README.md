# Analisi descrittive dell'11-12 nov 2026 (dopo gli esiti)

Questi script NON decidono nessun esito. Sono stati scritti e lanciati l'11 e il 12 nov 2026, DOPO aver visto i risultati
finali (complesso 5 x 500 ns, libero 10 x 1 us). Sono di due tipi, sempre etichettati nei testi:

- POST HOC (d1-d8, d10 e parte di d11): descrivono che cosa c'e' dietro gli esiti degli scopi 1 e 2. Non sono criteri e
  non portano test.
- PREVISTE DAI CRITERI ma eseguite dopo gli esiti (d9 e il punto 1 di d11): interfaccia (legami H, ponti salini,
  superficie sepolta) e RMSF del libero oltre 500 ns. Le definizioni geometriche non erano nei criteri: sono fissate nel
  docstring di d9, prima di vedere i numeri.

Regola del revisore (12 nov): nessun confronto a fotogrammi uniti di repliche diverse. Per questo d3 e' superato da d8,
e il punto 6 di d11 (legame Asp40-Ala90 e involucro) non e' riportato in nessun testo. Le uscite (.txt) stanno nella cartella del progetto `risultati_finali_1us/descrittivo_posthoc/`.

| Script | Che cosa descrive | Ingressi |
|---|---|---|
| `d1_serie.py` | serie per fotogramma: RMSD dell'epitopo con i due fit, scostamento per residuo, diedri phi/psi di 86-102, distanze scelte (fra cui Asp40 OD1/OD2 - Ala90 N) | xtc e gro di 09, pdb |
| `d2_riassunto.py`, `d3_dentro_fuori.py` | libero dentro e fuori dall'involucro secondario: per replica e per residuo; correlazione di rango con l'RMSF di 87-101 (descrittiva) | serie di d1, json di 13 e 14, cache di 10 |
| `d4_isteresi.py` | uscite dall'involucro per replica e conteggi con permanenza minima di 1 ns (registrata), 2 e 5 ns | serie di d1, json di 13 |
| `d5_run2_prima_dopo.py` | contatti con catena pesante e leggera e punta 90-93 nelle finestre 50-300 e 300-500 ns | cache di 10, json di 17 |
| `d6_asp40_ala90.py` | occupanza del legame H del cristallo Asp40 - Ala90 N | serie di d1, pdb |
| `d7_punta_core.py` | scostamento dei Ca dell'epitopo dopo fit sul core (complesso e libero); RMSF dei residui 36-40 | cache di 10, json di 14, pdb |
| `d8_diedri_per_replica.py` | diedri phi/psi di 86-102 nel libero, dentro e fuori dall'involucro secondario NELLA STESSA REPLICA, e confronto fra repliche; corregge d3 | serie di d1, json di 13, pdb |
| `d9_interfaccia.py` | PREVISTA DAI CRITERI: legami H, ponti salini e superficie sepolta all'interfaccia per replica del complesso (Shrake-Rupley, sonda 1.4 A, raggi di Bondi) | xtc e gro di 09, pdb |
| `d10_persistenza_scan.py` | repliche del libero con un ciclo uscita-rientro in funzione della permanenza richiesta; tratti consecutivi fuori dall'involucro; ultimo rientro; cicli nella sola finestra 500-1000 ns | serie di d1, json di 13 |
| `d11_verifiche.py` | RMSF del libero oltre 500 ns (PREVISTA DAI CRITERI); differenza rilevabile in 138-147; RMSF dei residui 36-41; frazione dei fotogrammi con la seconda disposizione di phi95/phi97; correlazione di rango con scarto 100 ns; Welch bilaterale per 87-101 (controllo di robustezza) | cache di 10, json di 13 e 14, serie di d1 |
| `genera_tabelle.py` | tabelle del supplementare e del report, generate dai json registrati e dai json di d9, d10 e d11 (nessun ricalcolo) | json di 11-14 e 17, json di d9-d11, pdb |

Come sono stati lanciati: da una cartella di lavoro con `trj_tutti/` (xtc e gro di 09), `cache/` (di 10), `risultati/`
(di 11-17), `descrittivo/` (uscite di d1) e `mrka-fab_completo.pdb`. d1-d3 si lanciano da quella cartella senza
argomenti; d4, d5 e d10 prendono la cartella come argomento; d6-d9, d11 e genera_tabelle prendono i percorsi come argomenti
(vedi il docstring di ognuno). Le uscite di d8-d11 stanno in `descrittivo/posthoc12/`. Sulla workstation la cartella equivalente e' `complesso_fab/analisi`, con
`traiettorie/` al posto di `trj_tutti/`.

Versione degli script di analisi usata: analisi_fab 2026-11-06 (comune_fab.py 3d90ec74).
