# Preparazione – MD classiche apo

| File | Scopo |
|---|---|
| `xray_apbs.in` | Input APBS per il potenziale elettrostatico della struttura X-ray protonata a pH 7.4 (pqr PARSE da PDB2PQR 3.7.1 + PROPKA). Usato per confrontare la carica/elettrostatica X-ray (−5) e NMR (−3). |

Workflow (comandi, non script):

1. `pdb2pqr` (PROPKA, pH 7.4) su X-ray e medoide NMR → stati di protonazione.
2. Rimozione di ACE/NME; `gmx pdb2gmx -ff amber14sb -water tip3p -ignh -ss -ter` (prima il ponte S–S, poi i due terminali).
3. `editconf` dodecaedro `-d 1.3`, `solvate` con `spc216.gro`, `genion -neutral` NaCl 0.15 M.
4. Minimizzazioni min_c1–c3.

Gli script/mdp originali di questi passi sono sul PC (cartella `mrkA/`) e non sono ancora in questo repository.
