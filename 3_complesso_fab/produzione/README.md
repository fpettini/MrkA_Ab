# Produzione – complesso MrkA–Fab e libero

```
03_test_velocita.sh   (ns/giorno, 2 vs 4 fs, Fab vs Fv)
05_equil_s1_s4.sh     (tronco comune, +05a dopo s1)
07_test_deriva.sh     (scelta del passo: esito 2 fs)
06_repliche.sh        (equil_s5 + produzione di tutte le repliche, controlli automatici)
```

| Script / file | Scopo |
|---|---|
| `03_test_velocita.sh` | Test di velocità (NVT con restraint) di complesso, libero e fv a 2 e 4 fs; serve solo a decidere Fab/Fv e il passo. |
| `05_equil_s1_s4.sh` | Equilibrazione comune s1 (NVT, 1 fs) – s4 (NPT), restraint 1000→1000→400→100; conta gli avvisi LINCS e si ferma se ce ne sono. |
| `05a_controlla_s1.py` | Controllo dopo s1 dei contatti di Asp40 (atomo con forza massima a fine minimizzazione). |
| `07_test_deriva.sh` | Deriva dell'energia conservata e temperature per 2 fs, 3 fs e varianti a 4 fs con HMR; criteri del revisore fissati prima. Esito: 2 fs. |
| `06_repliche.sh` | Repliche indipendenti da s4: equil_s5 20 ns (velocità nuove, seme r·11111) e produzione a tratti (`BLOCCO`), code per GPU, `TURNI=alternati`, regola LINCS, controllo T (AVVISO/STOP) a fine tratto, `controlla` per lo stato. Versione in uso: sha256 8077b6ec…. |
| `mdp/equil_s1-s4.mdp` | Equilibrazione con restraint. |
| `mdp/equil_s5.mdp.in`, `mdp/md.mdp.in` | Template di s5 e produzione (06 sostituisce passo, HMR, semi, durata, uscite). |
| `mdp/test.mdp.in` | Template del test di velocità (03). |

Lancio usato (2 nov 2026, 4 GPU):
`DT=2 NT=60 GPU_COMPLESSO="<4 id>" GPU_LIBERO="" TURNI=alternati setsid nohup bash 06_repliche.sh > 06.out 2>&1 < /dev/null &`
