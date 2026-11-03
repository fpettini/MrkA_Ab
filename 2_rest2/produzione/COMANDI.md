# REST2 su MrkA: comandi per la workstation xe7745-1

xe7745-1 ha 8 RTX PRO 6000 Blackwell e 2 EPYC 9655 (192 core fisici). Tutti i comandi partono da `/mnt/f.pettini/mrkA`, e gli script stanno in `rest2/`.
`<sistema>` è la cartella di un sistema, cioè quella che contiene `equil_s4.gro` e `processed.top` (X-ray: `rest2/xray`).

## 0. Copiare la cartella rest2 (da PowerShell sul PC)

Sul PC, `rest2` contiene già gli script aggiornati e `xray/` con `equil_s4.gro`, `processed.top` e `amber14sb.ff`.

```powershell
ssh f.pettini@<indirizzo_xe7745> mkdir -p /mnt/f.pettini/mrkA
cd C:\Users\francesco.pettini\Documents\project\mrkA
scp -r rest2 f.pettini@<indirizzo_xe7745>:/mnt/f.pettini/mrkA/
```

## 1. Compilare GROMACS 2024.6 + PLUMED 2.10.2 (una volta sola)

```bash
cd /mnt/f.pettini/mrkA
df -h /mnt/f.pettini
nohup bash rest2/00_compila_gromacs_plumed.sh /mnt/f.pettini/sw > rest2/compila.log 2>&1 &
tail -f rest2/compila.log
```

- Su xe7745-1 `/usr/local/cuda-12.8` ha solo le librerie, senza nvcc, e il CUDA 13.1 non va bene per GROMACS 2024. Lo script lo vede e mette il CUDA 12.8 di conda nell'ambiente `gmx2024-plumed-cuda`, insieme a PLUMED con MPI, Open MPI, cmake e gcc 14. Se trova un nvcc in `/usr/local/cuda-12.8/bin`, usa quello di sistema e l'ambiente `gmx2024-plumed`.
- Lo script non usa torch, non tocca l'ambiente base e riusa l'ambiente se esiste già.
- Durante la compilazione nvcc stampa molte volte "incompatible redefinition for option 'compiler-bindir'". Non è un problema: è lo stesso compilatore indicato due volte.
- Installa in `/mnt/f.pettini/sw/gromacs-2024.6-plumed`.
- Alla fine deve stampare `OK: gmx_mpi mdrun ha -hrex`.
- Crea `rest2/attiva_gmx2024.sh`. Da quel momento, in ogni nuova shell, prima di tutto:

```bash
source rest2/attiva_gmx2024.sh
```

**Perché 2024.6 e non 2024.3:**
- La 2024.3 non compila per le Blackwell RTX. L'ho riprodotto: ptxas si ferma sul kernel di pruning con "Value of threads per SM ... is out of range", come segnalato sul forum GROMACS. La correzione è nella 2024.6.
- La 2024.4 corregge anche le interazioni mancanti vicino al cutoff con il pruning su GPU (issue 5138).
- La patch PLUMED con `-hrex` esiste solo per la 2024.3. Sulla 2024.6 si applica senza conflitti, con spostamenti di una riga, e le modifiche di PLUMED risultano identiche.

## 2. Topologia espansa di ogni sistema

Per l'X-ray non serve: `rest2/xray/processed.top` è già fatta e vale anche qui.
Per un nuovo sistema (NMR con la sequenza X-ray):

```bash
source rest2/attiva_gmx2024.sh
cd <sistema>
cp -r ../xray/amber14sb.ff .
sed -e "s/@NSTEPS@/0/" -e "s/@SEED@/1/" ../rest2.mdp.in > pp.mdp
gmx_mpi grompp -f pp.mdp -c equil_s4.gro -p topol.top -pp processed.top -o pp.tpr
cd -
```

## 3. Test con tutte le repliche a lambda = 1 (pochi minuti)

Le Hamiltoniane sono identiche, quindi ogni scambio deve essere accettato.

```bash
source rest2/attiva_gmx2024.sh
nvidia-smi topo -m > rest2/topo.txt
bash rest2/04_prepara_rest2.sh rest2/xray 8 310.15 0.1 test_lam1_b
bash rest2/05_lancia_rest2.sh rest2/xray/test_lam1_b
python3 rest2/06_controlla_scambi.py rest2/xray/test_lam1_b
grep -H Performance rest2/xray/test_lam1_b/r*/rest2.log
```

Senza altri argomenti, 05 usa una GPU per replica (`-gpu_id 01234567`) e 24 thread per replica. Il rango i va sui core 24·i…24·i+23, quindi le repliche 0-3 stanno sul primo processore e le 4-7 sul secondo.

Controlla tre cose:
- accettazione **1.00** su tutte le coppie;
- in `r0/rest2.log`, la riga sotto "Mapping of GPU IDs" mette una replica per GPU;
- i ns/giorno per replica.

Mandami `topo.txt`, l'output di 06 e le righe Performance. Oppure copia `test_lam1_b` e `topo.txt` sul PC e li analizzo io. `topo.txt` dice se le GPU 0-3 sono davvero collegate al primo processore.

## 4. Taratura della scala (1 ns)

```bash
bash rest2/04_prepara_rest2.sh rest2/xray 16 480 1 taratura_16x480
bash rest2/05_lancia_rest2.sh rest2/xray/taratura_16x480
python3 rest2/06_controlla_scambi.py rest2/xray/taratura_16x480
```

Risultati X-ray (2026-10-20): 8 repliche fino a 480 K danno il 6-11% e nessun viaggio completo; 16 repliche (2 per GPU, 12 thread ciascuna) danno il 35-48%. Scala scelta: **16x480**.

L'obiettivo è il 20-30% di accettazione su ogni coppia.
- Se una coppia è sotto il 15%, aumenta le repliche (12 o 16, con più repliche per GPU) o abbassa Tmax.
- Se sono tutte sopra il 40%, si può salire con Tmax.

Mandami l'output e decidiamo la scala.

## 5. Produzione

```bash
bash rest2/04_prepara_rest2.sh rest2/xray 16 480 200 rest2_16x480
nohup bash rest2/05_lancia_rest2.sh rest2/xray/rest2_16x480 > /dev/null 2>&1 &
python3 rest2/06_controlla_scambi.py rest2/xray/rest2_16x480    # in qualsiasi momento
```

- Dopo il lancio puoi chiudere la connessione. 05 avvia mpirun con `setsid`, in una sessione sua: mpirun di Open MPI ignora `nohup` e, quando il terminale si chiude, ferma tutte le repliche con TERM. Era questo a fermare la produzione del 23 ottobre.
- Per fermare il run: `pkill -TERM -u $USER gmx_mpi`. Ctrl+C non lo ferma più.
- Se il run si ferma, rilancia lo stesso `05_lancia_rest2.sh`: riparte da `rest2.cpt`.
- `r0/rest2.xtc` è l'ensemble a 310.15 K (lambda = 1). Le cartelle tengono l'Hamiltoniana, le coordinate si scambiano.
- I primi 20 ns vanno scartati: si parte da `equil_s4`, con il backbone ancora vincolato.
- Un sistema alla volta usa tutte e 8 le GPU. Il secondo sistema si lancia dopo.
- **Rumore della GPU:** a lambda = 1 su GPU l'accettazione è 0.88-0.96 invece di 1.00; su CPU è 1.00. La differenza viene dagli arrotondamenti della GPU (qualche decimo di kT sull'energia della stessa configurazione). È piccola rispetto agli scambi veri, ma va citata nei metodi.

## Note

- **Argomenti di 05:** `bash rest2/05_lancia_rest2.sh <run> [gpu_id] [thread_per_replica]`. `gpu_id` può essere anche `cpu`, oppure `gputasks=<2 cifre per replica>` per decidere replica per replica (PP e PME, in ordine di replica).
- **Disconnessione:** `setsid` basta se logind non uccide i processi al logout. Si controlla con
  `busctl get-property org.freedesktop.login1 /org/freedesktop/login1 org.freedesktop.login1.Manager KillUserProcesses`:
  se risponde `b true`, serve `loginctl enable-linger` e il lancio con `systemd-run --user`.
- **Update su GPU:** con PLUMED non è possibile (il log dice "PLUMED is not supported"), quindi update e vincoli girano sulla CPU. Per questo conta dare a ogni replica core suoi.
- **Atomi caldi:** 1140 atomi in 85 residui (59-77, 87-101, 112-118, 124-133, 138-147, 159-169, 180-192), carica totale −4. Per cambiarli: `--caldi` in `03_marca_caldi.py`, oppure la variabile `CALDI` in cima al file.
- **Carica netta:** con lambda < 1 il box non è più neutro, perché le cariche calde sono moltiplicate per sqrt(lambda). grompp lo segnala con un avviso, e lo script ammette solo quell'avviso. La replica a lambda = 1 resta neutra.
- **Cosa ho provato io:**
  - su CPU, con la tua topologia X-ray e GROMACS 2024.3 + PLUMED: lambda = 1 dà energie identiche alla topologia non scalata; lambda = 0.5 su tutti gli atomi dimezza esattamente diedri, LJ, Coulomb, termini 1-4 e correzione di dispersione, e lascia uguali legami e angoli; REST2 a 2 repliche con `-hrex` funziona;
  - per la 2024.6: patch PLUMED applicata, compilazione completa per sm_120 con il CUDA 12.8 di conda, `-hrex` presente, REST2 a 2 repliche su CPU (lambda = 1: 5 scambi su 5 accettati).
- **Cosa non ho potuto provare:** l'esecuzione sulle GPU. Lo verifica il punto 3.

## 6. Sistema NMR con la sequenza X-ray

Stesso protocollo dell'X-ray: 16 repliche, 310.15-480 K, stessi residui caldi, 200 ns. Non servono il test a lambda = 1 né la taratura.
Controlli fatti sui file del PC (2026-10-24):
- la proteina in `nmr/topol.top` è identica a quella X-ray (atomi, cariche, legami, diedri, ponte S-S 21-656, posre.itp);
- 1140 atomi caldi con carica -4, come l'X-ray;
- mdp di equilibrazione identici, equil_s4 a 310.1 K e 1 bar;
- il box è più grande: 49957 atomi contro 47310.

Da PowerShell sul PC:

```powershell
cd C:\Users\francesco.pettini\Documents\project\mrkA\rest2
scp -r nmr f.pettini@<indirizzo_xe7745>:/mnt/f.pettini/mrkA/rest2/
```

Su xe7745-1:

```bash
cd /mnt/f.pettini/mrkA
source rest2/attiva_gmx2024.sh
pgrep -a gmx_mpi        # deve essere vuoto

cd rest2/nmr
sed -e "s/@NSTEPS@/0/" -e "s/@SEED@/1/" ../rest2.mdp.in > pp.mdp
gmx_mpi grompp -f pp.mdp -c equil_s4.gro -p topol.top -pp processed.top -o pp.tpr
cd ../..

bash rest2/04_prepara_rest2.sh rest2/nmr 16 480 200 rest2_16x480
nohup bash rest2/05_lancia_rest2.sh rest2/nmr/rest2_16x480 > /dev/null 2>&1 &
python3 rest2/06_controlla_scambi.py rest2/nmr/rest2_16x480    # dopo circa un'ora
```

- La prima riga di 04 deve dire `1140 atomi caldi in 85 residui ..., carica -4.000`.
- Dopo un'ora tutte le coppie dovrebbero stare vicino a 0.4, come l'X-ray. Se una scende sotto 0.20, fermati e mandami l'output.

## 7. REST3 (preparato il 24 ottobre; si lancia solo se lo decidiamo, a GPU libere)

Rispetto a REST2 cambia solo il vdW fra atomi caldi e acqua/ioni, moltiplicato per kappa in ogni replica (Zhang, Liu, Chen, JCTC 2023, 19, 1602). Atomi caldi, scala di lambda, mdp e lancio restano quelli di REST2.
- Nell'articolo tutta la proteina è calda. Qui sono caldi solo i loop, quindi kappa va solo fra atomi caldi e solvente: il vdW fra loop e core resta x sqrt(lambda) come in REST2. Lo script del loro tutorial scalerebbe anche loop-core.
- r0-r2 hanno kappa = 1 e la topologia identica a REST2, quindi r0 resta l'Hamiltoniana non scalata. Kappa cambia l'efficienza del campionamento, non l'ensemble di r0.
- I file del run si chiamano `rest2.*` anche con REST3, perché 05 usa `-deffnm rest2`.

Scala di kappa con kappa_max = 1.05 (curva kappa(lambda) del tutorial REST3, che coincide con la formula dell'articolo per 16 repliche; al nostro lambda minimo 0.646 la curva dà 1.049):

| replica | 0-2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| kappa | 1 | 1.0015 | 1.0030 | 1.0065 | 1.0100 | 1.0150 | 1.0200 | 1.0237 | 1.0272 | 1.0312 | 1.0351 | 1.0401 | 1.0451 | 1.0500 |

Da PowerShell sul PC:

```powershell
cd C:\Users\francesco.pettini\Documents\project\mrkA\rest2
scp 07_kappa_rest3.py 08_prepara_rest3.sh COMANDI.md f.pettini@<indirizzo_xe7745>:/mnt/f.pettini/mrkA/rest2/
```

### 7.1 Test degli scambi (1 ns, pochi minuti)

```bash
cd /mnt/f.pettini/mrkA
source rest2/attiva_gmx2024.sh
pgrep -a gmx_mpi        # deve essere vuoto
bash rest2/08_prepara_rest3.sh rest2/xray 16 480 1 1.05 taratura_rest3_k105
cat rest2/xray/taratura_rest3_k105/r15/kappa.txt
bash rest2/05_lancia_rest2.sh rest2/xray/taratura_rest3_k105
python3 rest2/06_controlla_scambi.py rest2/xray/taratura_rest3_k105
```

- 08 deve dire `1140 atomi caldi in 85 residui ..., carica -4.000`, e la riga della replica 15 deve finire con `kappa 1.0500`.
- `kappa.txt` deve dire `20 tipi caldi x 4 tipi del solvente (Cl-_tip3p, HW_tip3p, Na+_tip3p, OW_tip3p; molecole SOL, NA, CL) = 80 coppie`.
- La stessa taratura con REST2 dava 0.35-0.48 su tutte le coppie. Se una coppia scende sotto 0.20, fermati e mandami l'output.

### 7.2 Produzione X-ray, con un controllo a 100 ns

```bash
bash rest2/08_prepara_rest3.sh rest2/xray 16 480 200 1.05 rest3_16x480_k105
nohup bash rest2/05_lancia_rest2.sh rest2/xray/rest3_16x480_k105 > /dev/null 2>&1 &
```

Quando 06 dice che il run ha passato i 100 ns (circa 10 ore), taglia r0 e r15 a 100 ns senza fermare il run:

```bash
cd /mnt/f.pettini/mrkA/rest2/xray
for r in r0 r15; do
  mkdir -p rest3_k105_100ns/$r
  gmx_mpi trjcat -f rest3_16x480_k105/$r/rest2.xtc -e 100000 -o rest3_k105_100ns/$r/rest2.xtc
done
python3 ../06_controlla_scambi.py rest3_16x480_k105 > rest3_k105_100ns/scambi.txt
```

Copia la cartella `rest3_k105_100ns` sul PC in `rest2\xray\` e dimmelo. Il controllo è 02 sulla finestra 20-100 ns, con una cache separata:

```
python analisi_rest2/01_estrai.py . analisi_rest2/cache_rest3 rest3_k105_100ns
python analisi_rest2/02_controlli.py . analisi_rest2/cache_rest3 analisi_rest2/risultati_rest3
```

Criteri, fissati prima di vedere i dati:
- REST2 X-ray nella stessa finestra (20-100 ns) segnalava r15 più compatto di r0 in 87-101 (-3.5%), 112-118 (-3.9%) e 124-133 (-1.3%).
- Si prosegue fino a 200 ns se: tutte le coppie restano sopra 0.20, 87-101 non ha più il flag di compattazione e il core di r15 resta sotto 2 Å. Poi si lancia l'NMR con gli stessi comandi su `rest2/nmr`.
- Se 87-101 resta compattato, kappa 1.05 non risolve il problema: si ferma il run (`pkill -TERM -u $USER gmx_mpi`) e si torna ad allungare REST2. Kappa non si alza alla cieca.
- I viaggi completi si confrontano con REST2 (240 in 200 ns), ma non decidono da soli.

### Perché kappa non è tarato come nell'articolo

- L'articolo ha scelto kappa (1.06 a 500 K) confrontando una proteina disordinata, tutta calda, con una MD normale a 500 K, su raggio di girazione, distanza fra le estremità ed elica.
- Per MrkA quel confronto non si può fare: a 480 K una MD normale srotolerebbe anche il core, che in REST2 resta a 310 K.
- Una taratura nostra sul raggio di girazione costerebbe almeno 100 ns per ogni valore di kappa, perché in REST2 la compattazione di 87-101 in r15 si forma in 50-100 ns. Per questo si parte dalla curva pubblicata e si controlla a 100 ns.

### Cosa ho provato io (CPU, GROMACS 2024.3 + PLUMED 2.10.2, topologia X-ray della produzione)

- Kappa = 1 a lambda 0.646: energia potenziale su `equil_s4.gro` identica a REST2 (-639428.56 kJ/mol); matrice LJ uguale entro 1.3e-6, cioè l'arrotondamento in singola precisione.
- Kappa 1.05 e 1.10: cambia solo il vdW caldi-solvente, esattamente x1.05 e x1.10 (da -655.21 a -687.97 e -720.74 kJ/mol). Caldo-caldo, caldo-freddo, freddo-solvente e solvente-solvente restano identici.
- Cambia anche la correzione di dispersione (-10.4 kJ/mol con kappa 1.05), perché dipende dalla C6 media del sistema: fa parte dell'Hamiltoniana di ogni replica.
- 08 con 16 repliche: nessun avviso nuovo di grompp; `processed_caldi.top` e r0 identiche byte per byte alla produzione REST2; lambda identiche a `rest2_16x480/scala.txt`.
- 2 repliche (r14 e r15) con `-hrex`: girano e tentano gli scambi.
- **Non provato:** le GPU e l'accettazione vera. Lo dice il punto 7.1.

## 8. Allungare REST2: X-ray e NMR fino a 400 ns

Perché 400 ns: nell'X-ray i loop più lenti (59-77 e 180-192) cambiano popolazione fra le due metà di 0.20. Per avere in ogni run un errore di circa 0.07, e quindi un errore combinato sotto la soglia di 0.10 del piano, servono circa 380 ns dopo lo scarto. È una stima da una sola coppia di metà, quindi incerta: se a 400 ns non basta, si allunga ancora con lo stesso comando.

Un sistema alla volta, come la produzione. Ordine scelto da francesco il 24 ottobre: l'NMR si ferma a 200 ns e si allunga prima l'X-ray. Così l'X-ray a 400 ns dice se la stima regge, mentre si analizzano i 200 ns dell'NMR. Nei due ordini i run arrivano a 400 ns alla stessa data; cambia solo quale dei due è pronto prima.

Quando l'NMR ha finito i 200 ns (si ferma da solo), si lancia a mano l'X-ray:

```bash
cd /mnt/f.pettini/mrkA
source rest2/attiva_gmx2024.sh
pgrep -a gmx_mpi        # deve essere vuoto: l'NMR ha finito i 200 ns
for s in nmr xray; do python3 rest2/06_controlla_scambi.py rest2/$s/rest2_16x480 > rest2/$s/rest2_16x480/scambi_200ns.txt; done
head -1 rest2/nmr/rest2_16x480/scambi_200ns.txt     # ultimo scambio al passo 99999500 (199.999 ns), come l'X-ray
R=rest2/xray/rest2_16x480
bash rest2/09_allunga_rest2.sh $R 400
nohup bash rest2/05_lancia_rest2.sh $R > /dev/null 2>&1 &
```

- 09 deve stampare `r0: da 200 a 400 ns` e così via fino a r15. Si rifiuta se c'è un gmx_mpi attivo o se manca un `rest2.cpt`. Il tpr vecchio resta come `rest2_fino_200ns.tpr`.
- 05 riparte dal checkpoint e aggiunge in coda a `rest2.xtc`, `rest2.log` e `rest2.edr`. L'ho provato su CPU con 2 repliche e `-hrex`: riparte dal passo del checkpoint e gli scambi continuano.
- I 200 ns dell'X-ray sono già sul PC (copia del 24 ottobre). Sulla workstation, dopo l'allungamento, si ritagliano con `gmx_mpi trjcat -f r0/rest2.xtc -e 200000 -o r0/rest2_200ns.xtc`.

Copia sul PC i file NMR a 200 ns (da PowerShell). L'NMR è fermo, quindi `rest2.xtc` è già quello dei 200 ns:

```powershell
cd C:\Users\francesco.pettini\Documents\project\mrkA\rest2\nmr
mkdir -Force rest2_16x480\r0, rest2_16x480\r15 | Out-Null
$W = "f.pettini@<indirizzo_xe7745>:/mnt/f.pettini/mrkA/rest2/nmr/rest2_16x480"
scp "$W/scambi_200ns.txt" rest2_16x480\scambi_200ns.txt
scp "$W/r0/rest2.xtc" rest2_16x480\r0\rest2.xtc
scp "$W/r15/rest2.xtc" rest2_16x480\r15\rest2.xtc
```

Quando l'X-ray arriva a 400 ns e si ferma da solo, si copiano sul PC `r0/rest2.xtc`, `r15/rest2.xtc` e l'output di 06, e si analizza. L'NMR si allunga poi con gli stessi comandi e `R=rest2/nmr/rest2_16x480`. Il revisore ha accettato l'allungamento a condizione che l'analisi confermativa si faccia a 400 ns su entrambi i sistemi: l'X-ray a 400 ns può portare l'NMR oltre i 400 ns, non fermarlo prima.

- Tempo: circa 20 ore per sistema. L'X-ray fa 236 ns/giorno per replica (misurati). L'NMR fa circa 225 ns/giorno: l'ho dedotto dagli orari dei messaggi del 24 ottobre, non dal log.
- Alla fine servono `r0/rest2.xtc` e `r15/rest2.xtc` di entrambi i sistemi, con l'output di 06. Si rilanciano 01-03 come in `analisi_rest2/LEGGIMI.md`: lo scarto resta 20 ns e la finestra arriva alla fine della traiettoria. 03 avvisa se i due run non hanno la stessa lunghezza.
- Se dopo l'analisi dei 200 ns scegliamo REST3, si ferma l'allungamento con `pkill -TERM -u $USER gmx_mpi`: si perdono solo le ore già fatte.

**Com'è andata (27 ottobre).** L'X-ray è stato allungato a 400 ns dal 25 ottobre (`scambi_400ns.txt` del 26), poi l'NMR (`scambi_400ns.txt` del 27, ultimo scambio a 399.999 ns, 403 viaggi). Il 27 ottobre francesco ha deciso che l'esito a 400 ns è definitivo: nessun altro allungamento. REST3 (sezione 7) non è stato lanciato.
