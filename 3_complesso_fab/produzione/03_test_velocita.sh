#!/usr/bin/env bash
# Test di velocita' (ns/giorno) di complesso, libero e fv, con passo 2 fs (come le apo) e 4 fs con hydrogen mass
# repartitioning (mass-repartition-factor = 3). NVT a 310.15 K con restraint 1000 sugli atomi pesanti (solo per stabilita'
# partendo dalla minimizzazione: il costo per passo e' lo stesso della produzione). Un run per GPU, tutti in parallelo.
# Serve solo a decidere Fab/Fv e 2/4 fs: nessun risultato scientifico da questi run.
#
# Uso: bash 03_test_velocita.sh "0 1 2 3 4 5"      (GPU da usare, una per run, nell'ordine dei 6 test)
# Variabili: GMX (default gmx), NT (thread per run, default 40), PS (ps simulati per test, default 100)
set -euo pipefail
QUI=$(cd "$(dirname "$0")" && pwd)
read -r -a GPU <<< "${1:-0 1 2 3 4 5}"
GMX=${GMX:-gmx}; NT=${NT:-40}; PS=${PS:-100}
mkdir -p "$QUI/test_velocita"; cd "$QUI/test_velocita"
i=0
for S in complesso libero fv; do
  [ -f "$QUI/$S/min_c3.gro" ] || { echo "manca $S/min_c3.gro: lancia prima 02"; exit 1; }
  for DT in 2 4; do
    T=${S}_${DT}fs
    NS=$(( PS * 1000 / DT ))
    HMR=1; [ $DT = 4 ] && HMR=3
    sed -e "s/@DT@/0.00$DT/" -e "s/@NSTEPS@/$NS/" -e "s/@HMR@/$HMR/" "$QUI/mdp/test.mdp.in" > $T.mdp
    $GMX grompp -f $T.mdp -c "$QUI/$S/min_c3.gro" -r "$QUI/$S/min_c3.gro" -p "$QUI/$S/topol.top" -n "$QUI/$S/index.ndx" \
         -o $T.tpr -maxwarn 0 > grompp_$T.log 2>&1 || { echo "grompp fallito: test_velocita/grompp_$T.log"; exit 1; }
    $GMX mdrun -deffnm $T -nt $NT -gpu_id ${GPU[$i]} -resethway -noconfout > mdrun_$T.out 2>&1 &
    i=$((i + 1))
  done
done
wait
{
  printf "%-14s %10s %12s\n" sistema atomi "ns/giorno"
  for S in complesso libero fv; do for DT in 2 4; do
    T=${S}_${DT}fs
    printf "%-14s %10s %12s\n" $T "$(sed -n 2p "$QUI/$S/min_c3.gro" | tr -d ' ')" "$(grep Performance: $T.log | awk '{print $2}')"
  done; done
} | tee risultati.txt
grep -l "LINCS WARNING\|Fatal error" *.log mdrun_*.out 2>/dev/null && echo "ATTENZIONE: avvisi nei log sopra" || echo "nessun avviso LINCS o errore"
