# Produzione – MD classiche apo

| File | Scopo |
|---|---|
| `mdp/equil_s1.mdp` | NVT 1 ns, dt 1 fs, restraint backbone 1000. |
| `mdp/equil_s2-s4.mdp` | NPT, restraint backbone 1000 → 400 → 100. |
| `mdp/equil_s5.mdp.in` | Template di equil_s5 (20 ns NPT senza restraint, velocità rigenerate): `eq.sh` sostituisce `${CONT} ${GVEL} ${GSEED} ${LSEED}`. |
| `mdp/md.mdp.in` | Template della produzione (continuazione da equil_s5, `ld-seed` per replica). |
| `eq_estensione_500ns.sh` | Versione più recente di `eq.sh` (cartella principale `mrkA/`, 19 ott): estende le 3 repliche da 200 a 500 ns con `convert-tpr -until 500000` e `mdrun -cpi -append`. Il resto è commentato perché già eseguito. |

Gli mdp sono identici per xray, nmr_ox, nmr_red (e per i sistemi di partenza di REST2).
