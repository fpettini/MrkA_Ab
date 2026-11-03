# Specifica degli script di analisi REST2 di MrkA (01, 02, 03)

Viene dal piano di analisi approvato il 24 ottobre 2026. Le soglie e le regole qui sotto sono fissate: non
vanno cambiate per far tornare un risultato. Se una regola risulta inapplicabile o ambigua, scegli
l'interpretazione più semplice, scrivila nel commento dello script e segnalala nel resoconto.

## Contesto

- Proteina: dominio pilinico di MrkA, residui wt 42-202 (161 residui), sequenza X-ray in entrambi i sistemi.
- Due run REST2 con la stessa Hamiltoniana, che partono da strutture diverse: `xray` (struttura X-ray) e
  `nmr` (modello 6 di 9HW9, mutato alla sequenza X-ray). Ogni run ha 16 repliche; si analizzano
  `r0` (310.15 K, lambda = 1) e `r15` (la più calda, lambda = 0.646).
- La domanda: per ogni loop, che frazione di frame somiglia all'X-ray e che frazione all'NMR, e se i due run
  danno le stesse popolazioni (convergenza).

## Dati (cartella del progetto `ROOT`, stessa struttura del PC)

- `ROOT/rest2/<sistema>/rest2_16x480/r0/rest2.xtc` e `.../r15/rest2.xtc`, con `<sistema>` in {xray, nmr}.
  Solo proteina (2286 atomi), un frame ogni 10 ps, 0-200 ns (20001 frame). Oggi esiste solo `xray`:
  il run `nmr` è ancora in produzione, quindi ogni script deve saltare con un messaggio un sistema
  senza traiettoria.
- **Le molecole nell'xtc sono spezzate dalle condizioni periodiche** (box dodecaedrico triclino; in
  80 frame su 101 controllati c'è un "legame" di backbone lungo fino a 123 Å). L'estrazione deve
  ricomporre la proteina.
- Topologia: `ROOT/rest2/<sistema>/conf.gro` (uscita di pdb2gmx, solo proteina, stesso ordine di atomi
  dell'xtc, numerazione 42-202). **Mai usare un .tpr come topologia in MDAnalysis**: rinumera i residui da 1.
- Riferimenti:
  - X-ray: `ROOT/rest2/xray/conf.gro`.
  - NMR: i 20 modelli di `ROOT/9hw9_20modelli.pdb` (record MODEL 1-20, catena A, residui 1-201 in
    numerazione del costrutto: **wt = costrutto + 27**, cioè Ser15 = Ser42, Cys16 = Cys43, Gln175 = Gln202).
    Hanno la sequenza NMR in 8 posizioni, ma si usa solo il backbone, che ha gli stessi nomi.
  - Controllo: `ROOT/rest2/nmr/conf.gro` (modello 6 con la sequenza X-ray, partenza del run nmr).

## Definizioni comuni

- CORE: i 43 Cα usati nelle MD classiche (lista `STRAND` in `analisi_convergenza/04_confronto_ristretto.py`):
  44-47, 78-85, 88, 94, 102-105, 107, 108, 122, 134-138, 148, 162, 170-179, 195-199.
- LOOP (backbone N, CA, C, O di ogni residuo):
  59-77, 87-101, 112-118 **senza 116 e 117**, 124-133, 138-147, 159-169, 180-192, 48-56.
  I primi sette sono riscaldati in REST2, 48-56 no.
- d(frame, R) per un loop: si sovrappone il frame alla struttura R sui Cα del CORE (Kabsch con
  traslazione), poi RMSD degli atomi del loop **senza** nuovo fit. Unità: Å.

## 01_estrai.py

Uso: `python analisi_rest2/01_estrai.py ROOT CACHE`

- Per ogni sistema presente e per r0 e r15: legge xtc + conf.gro, controlla che il numero e l'ordine degli
  atomi coincidano con conf.gro.
- Ricompone la proteina in ogni frame percorrendo la catena: ogni atomo va all'immagine minima rispetto
  all'atomo a cui è legato (N(i) → CA(i) → C(i) → N(i+1); O e H su C e N; SG su CA; scatola triclina, per
  esempio con `MDAnalysis.lib.distances.minimize_vectors`).
- Salva in `CACHE/<sistema>_<replica>.npz` solo gli atomi che servono: N, CA, C, O, H di backbone e SG di
  Cys43 e Cys87. Campi: `X` (frame × atomi × 3, float32, Å), `t` (ps), `names`, `resids`, `resnames` e
  `ref` (le stesse coordinate prese da conf.gro).
- Controllo obbligatorio dopo la ricomposizione: in ogni frame la distanza C(i)-N(i+1) massima è sotto
  2.0 Å. Se no, errore.
- Stampa per ogni file: frame, intervallo di tempo, passo, atomi.

## 03_popolazioni.py (analisi principale)

Uso: `python analisi_rest2/03_popolazioni.py ROOT CACHE OUT`

1. **Tabella dei riferimenti.** Stampa in `OUT/riferimenti.txt`, per ogni loop:
   - d(X-ray, modello NMR più vicino), con il numero del modello;
   - mediana di d(X-ray, modelli NMR);
   - dispersione NMR, cioè la mediana di d fra le 190 coppie di modelli;
   - d(X-ray, rest2/nmr/conf.gro).

   Stampa anche l'RMSD dei 43 Cα del CORE fra X-ray e rest2/nmr/conf.gro. Questi numeri devono
   riprodurre il calcolo di riferimento entro 0.02 Å:

   | loop | X-NMR più vicino [modello] | mediana X-NMR | dispersione NMR | X-modello 6 mutato |
   |---|---|---|---|---|
   | 59-77 | 2.13 [17] | 2.64 | 1.21 | 2.52 |
   | 87-101 | 2.12 [10] | 2.83 | 1.71 | 2.71 |
   | 112-118 | 1.10 [2] | 1.76 | 1.05 | 1.56 |
   | 124-133 | 1.45 [3] | 1.99 | 0.88 | 1.92 |
   | 138-147 | 2.88 [18] | 3.08 | 1.07 | 2.97 |
   | 159-169 | 1.49 [10] | 1.70 | 0.74 | 1.66 |
   | 180-192 | 1.24 [3] | 1.65 | 0.91 | 1.87 |
   | 48-56 | 1.24 [14] | 1.44 | 0.55 | 1.43 |

   Valore atteso per il CORE: 1.04 Å. Se non tornano, non ritoccare il metodo per farli tornare: segnala la
   differenza e la causa probabile.
2. **Distanze per frame**, su r0 di ogni sistema presente:
   - dX = d(frame, X-ray);
   - dN = minimo su k di d(frame, modello k), con l'indice del modello;
   - dN6 = d(frame, rest2/nmr/conf.gro).

   Salvale in `OUT/distanze_<sistema>.npz`, per le mappe 2D che verranno dopo.
3. **Assegnazione** (corretta il 24 ottobre: il riferimento principale è il modello 6, scelta di francesco).
   Per ogni loop sia D6 = d(X-ray, modello 6 mutato) del punto 1.
   - Se min(dX, dN6) > D6, il frame è "altro".
   - Altrimenti è "X" se dX < dN6, "N" se no.

   Assegnazione di controllo: stessa regola con dN al posto di dN6 e D = d(X-ray, modello NMR più vicino).
4. **Scarto iniziale: 20 ns fissi** (decisione di francesco del 24 ottobre; la regola precedente, il più lungo fra
   20 ns e il t0 di Chodera, sull'X-ray dava 181 ns per le fluttuazioni lente). I t0 di
   `pymbar.timeseries.detect_equilibration(A, nskip=10)` sulle serie 0/1 "è X" e "è N" a 100 ps si stampano solo
   come diagnostica. `--t0 NS` cambia lo scarto.
5. **Popolazioni ed errori** nella finestra [T0, 200 ns], per sistema e loop: frazioni pX, pN, pAltro.
   - Errore standard dai blocchi da 20 ns contati dalla fine: si usano solo blocchi interi e l'avanzo si
     toglie dall'inizio. SE = deviazione standard (ddof=1) delle n medie di blocco / √n.
   - Stampa anche il SE con blocchi da 5, 10, 20 e 30 ns. Se SE(30) > SE(20), scrivi
     "SE probabilmente sottostimato".
6. **Metà.** Dividi la finestra in due metà uguali e calcola pX e pN con il loro SE in ciascuna (blocchi da
   20 ns dentro ogni metà, stessa regola).
   - Le metà sono "coerenti" se |Δ| < 2·√(SE1² + SE2²) sia per pX sia per pN.
7. **Confronto fra i run**, solo se ci sono entrambi i sistemi: Δ = p(xray) − p(nmr) per pX e per pN, con
   SEc = √(SE_xray² + SE_nmr²). Esito per loop:
   - "indeterminato" se uno dei SEc è ≥ 0.10;
   - altrimenti "convergente" se |Δ| < 2·SEc per pX e per pN e le metà sono coerenti in entrambi i run;
   - altrimenti "non convergente in 200 ns".
8. Ripeti i punti 5-7 con l'assegnazione di controllo (20 modelli), in una sezione separata.
9. Scrivi `OUT/popolazioni.txt`, leggibile, con tutti i numeri usati, e `OUT/popolazioni.json`, con gli stessi
   numeri. Nessun numero calcolato altrove.

## 02_controlli.py (controlli preliminari)

Uso: `python analisi_rest2/02_controlli.py ROOT CACHE OUT [--t0 NS]` (default 20 ns), finestra [T0, 200 ns].

Per ogni sistema presente, su r0 e r15:
- RMSD dei 43 Cα del CORE dalla propria conf.gro (fit sul CORE): media e massimo. Flag se la media in r15
  supera 2.0 Å.
- Frazione β con `MDAnalysis.analysis.dssp.DSSP` (un frame ogni 1 ns) sui residui 44-199: media in r0 e
  in r15.
- Raggio di girazione dei Cα di ciascuno dei sette loop caldi: media e SE a blocchi (regola del punto 5) in
  r0 e r15. Flag se Rg(r15) < Rg(r0) − 2·√(SE0² + SE15²).
- Distanza SG43-SG87: media e massimo in r0 (controllo di integrità del ponte).

Scrivi `OUT/controlli.txt` e `OUT/controlli.json`. Il controllo degli scambi c'è già (06_controlla_scambi.py
sulla workstation) e qui non si rifà.

## Stile

- Come gli script esistenti in `analisi_convergenza/`: docstring in italiano con la riga "Uso:",
  argomenti posizionali ROOT/CACHE/OUT, numpy vettoriale, commenti brevi in italiano, nessuna classe.
- Si usano solo i pacchetti già installati (MDAnalysis 2.10, numpy, scipy, pymbar 4, matplotlib).
  Non installare niente.
