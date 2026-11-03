# Analisi REST2: script, analisi e risultati

Tutti i comandi vanno lanciati dalla cartella del progetto `mrkA`. I risultati finiscono in `analisi_rest2/risultati/`.

| Passo | Script | Analisi | Input | Output | Esito X-ray (24 ott) |
|---|---|---|---|---|---|
| 0 | `rest2/06_controlla_scambi.py` | Scambi: accettazione per coppia e viaggi completi per walker (controllo 1) | `rest2/<sistema>/rest2_16x480/r*/rest2.log`, `scala.txt` | a schermo | Accettazione 0.39-0.44, 240 viaggi, almeno 3 per walker |
| 1 | `analisi_rest2/01_estrai.py` | Ricompone la proteina spezzata dalle condizioni periodiche; salva il backbone e gli SG 43/87 di r0 e r15 | `rest2/<sistema>/rest2_16x480/r0/rest2.xtc`, `r15/rest2.xtc`, `rest2/<sistema>/conf.gro` | `analisi_rest2/cache/<sistema>_r0.npz`, `<sistema>_r15.npz` | Legame C-N massimo 1.47 Å |
| 2 | `analisi_rest2/02_controlli.py` | r15 contro r0: RMSD del core, frazione β, Rg dei loop caldi, S-S (controllo 2) | cache | `controlli.txt`, `controlli.json` | Core in r15 1.58 Å; flag Rg su 87-101 (-6%) e 124-133 (-3%) |
| 3 | `analisi_rest2/03_popolazioni.py` | Distanze dai riferimenti, popolazioni X/N/altro, metà, confronto fra i run (controllo 3 e analisi principale) | cache, `9hw9_20modelli.pdb`, `rest2/xray/conf.gro`, `rest2/nmr/conf.gro` | `riferimenti.txt`, `popolazioni.txt`, `popolazioni.json`, `distanze_<sistema>.npz` | Distanze del piano riprodotte; metà non coerenti in 59-77, 112-118, 159-169, 180-192 |

## Comandi per l'NMR

```
python3 rest2/06_controlla_scambi.py rest2/nmr/rest2_16x480
python analisi_rest2/01_estrai.py . analisi_rest2/cache
python analisi_rest2/02_controlli.py . analisi_rest2/cache analisi_rest2/risultati
python analisi_rest2/03_popolazioni.py . analisi_rest2/cache analisi_rest2/risultati
```

- 01 elabora da solo tutti i sistemi che hanno la traiettoria, quindi rifà anche l'X-ray. Serve, perché il confronto in 03 vuole le cache di entrambi i run.
- Lo scarto iniziale è di 20 ns fissi in 02 e 03 (decisione del 24 ottobre). `--t0 NS` lo cambia, solo per un controllo.
- Dal 25 ottobre 03 fa anche: (1) la sensibilità con scarto di 120 ns (`--t0-sens`, 0 la toglie); un loop è convergente solo se lo è con entrambi gli scarti (sezione ESITO FINALE); (2) il test di equivalenza |Δ| + t·SEc < 0.15 per pX e pN, fra i run e fra le metà, con t di Student a n−1 gradi di libertà; blocchi da 20 ns sotto 300 ns di finestra e da 40 ns sopra (`--blocco` li fissa), SE preso come il massimo fra il blocco scelto e i blocchi più lunghi con almeno 4 blocchi; (3) la sensibilità alla soglia "altro" con 0.75·D e 1.25·D, per l'assegnazione principale e quella di controllo (sezione RIEPILOGO ESITI). pX e pN sono frazioni di frame entro la soglia; per 87-101 e 112-118 conta solo il confronto fra i run.
- A 200 ns lo scarto di 120 ns lascia metà da 40 ns, cioè 2 blocchi e 1 grado di libertà (t = 12.7): con quello scarto l'equivalenza in pratica non si può dimostrare prima dei 400 ns.
- Provenienza: `provenienza.py` (stessa cartella degli script) mette in testa a ogni txt e sotto la chiave "provenienza" di ogni json lo sha256 dello script e degli input, data, versioni di Python e pacchetti e gli argomenti. 01 la scrive anche dentro le cache npz. Le cache fatte prima del 25 ottobre non l'hanno: vanno rifatte con 01. I percorsi sono scritti relativi a `<ROOT>`, `<CACHE>` e `<OUT>`; i txt di 02 e 03 riportano anche lo sha256 delle xtc da cui vengono le cache.
- Esiti di 03: "convergente" se tutte le prove sono equivalenti; "differenza rilevata" se almeno una prova è significativa e non equivalente; altrimenti "equivalenza non dimostrata". 03 avvisa quando i gradi di libertà minimi sono sotto 3.
- Risultati a 200 ns (X-ray e NMR, 25 ottobre) in `analisi_rest2/risultati_200ns/`.
- La finestra va dallo scarto alla fine della traiettoria, quindi dopo un allungamento 02 e 03 usano tutto. Il confronto fra i run vale solo a pari lunghezza: 03 lo segnala se le lunghezze differiscono.
- 01 accetta come terzo argomento il nome del run (default `rest2_16x480`). Per un run REST3 si usano una cache e una cartella di risultati separate, vedi `rest2/COMANDI.md`, sezione 7.
- Pacchetti: Python con MDAnalysis 2.10 (basta la 2.7 o successiva), numpy, scipy e pymbar 4. Non so se sulla workstation ci sono.
- Tempi nel mio ambiente: 01 circa 36 s per sistema, 02 circa 22 s, 03 circa 40 s.
