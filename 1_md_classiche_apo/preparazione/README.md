# Preparazione – MD classiche apo

Gli `eq_<sistema>.sh` sono gli `eq.sh` delle cartelle `xray/`, `nmr_ox/`, `nmr_red/` sul PC; ognuno si lancia
dalla cartella del sistema e fa tutto il tronco comune e le repliche. Le parti già eseguite sono commentate (`#`).

| File | Scopo |
|---|---|
| `eq_xray.sh` | X-ray ossidato: minimizzazioni min_c1–c3 (gmx_d), equil_s1–s4, poi 3 repliche (equil_s5 20 ns con `gen-vel`, semi r·11111; produzione). |
| `eq_nmr_ox.sh` | NMR ossidato: in più `pdb2gmx -ff amber14sb -water tip3p -ignh -ss -ter`, `editconf -bt dodecahedron -d 1.3`, `solvate` (spc216), `genion -neutral -conc 0.15`, `make_ndx`; poi come sopra. |
| `eq_nmr_red.sh` | NMR ridotto (senza S–S): stesso protocollo. |
| `mdp/min_c1-3.mdp` | Minimizzazioni steepest descent con restraint sul backbone. |
| `mdp/ions.mdp` | mdp per genion. |
| `xray_apbs.in` | Input APBS per il potenziale elettrostatico dell'X-ray protonata a pH 7.4 (pqr PARSE da PDB2PQR 3.7.1 + PROPKA), usato per confrontare la carica X-ray (−5) e NMR (−3). |

La protonazione (PDB2PQR/PROPKA a pH 7.4) è stata fatta da riga di comando; non c'è uno script.
