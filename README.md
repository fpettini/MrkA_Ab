# MrkA – script delle simulazioni e delle analisi

Script e file di input (mdp, template) usati per le dinamiche molecolari su MrkA
(subunità fimbriale di tipo 3 di *K. pneumoniae*, UniProt P12267) e sul complesso
MrkA–Fab. Niente traiettorie, strutture o risultati: solo ciò che serve a rifare i calcoli.

| Cartella | Filone | Stato |
|---|---|---|
| [`1_md_classiche_apo/`](1_md_classiche_apo/) | MD classiche di MrkA libera, X-ray vs NMR (3 × 500 ns), convergenza, RMSF vs B-factor, contatti nel cristallo | concluso |
| [`2_rest2/`](2_rest2/) | REST2 (16 repliche, 310–480 K) su X-ray e NMR con la sequenza X-ray, analisi delle popolazioni dei loop | concluso (400 ns) |
| [`3_complesso_fab/`](3_complesso_fab/) | Complesso MrkA–Fab intero e MrkA libera (amber19sb + OPC, 298 K, KCl 0.2 M) | produzione in corso |

Ogni filone è diviso in `preparazione/`, `produzione/` e `analisi/`, ciascuna con un README
che descrive il workflow e lo scopo di ogni script.

## Workflow generale

```
struttura di partenza ──► preparazione (protonazione, topologia, solvatazione, minimizzazione, controlli)
                       ──► produzione (equilibrazione, repliche / REST2, controlli in corsa)
                       ──► analisi (estrazione, convergenza, grandezze per residuo/loop)
```

## Note d'uso

- Sulla workstation gli script di un filone stanno tutti nella **stessa cartella** (p.es. `complesso_fab/`
  con la sottocartella `mdp/`). La divisione in preparazione/produzione/analisi è solo per documentazione:
  per rilanciarli, copiare insieme i file del filone e rimettere gli mdp in `mdp/`.
- Gli script si lanciano dalla cartella del progetto `mrkA` (vedi il blocco "Uso" in testa a ogni script).
- Software: GROMACS 2026 (complesso Fab), GROMACS 2024.6 + PLUMED 2.10.2 (REST2), Python con
  MDAnalysis ≥ 2.7, numpy, scipy, pymbar 4, OpenMM, gemmi, mdtraj (dove indicato).
- Gli script di analisi scrivono la **provenienza** (sha256 di script e input, versioni dei pacchetti)
  in testa ai risultati tramite `provenienza.py`.
- Numerazione dei residui: wt 42–202 (costrutto) in tutti gli script; vedi i README dei singoli filoni per le conversioni.
