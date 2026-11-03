#!/usr/bin/env bash
# Equilibrazione comune equil_s1-s4 di complesso e libero (s1 a 1 fs, s2-s4 a 2 fs, senza HMR), restraint sugli atomi
# pesanti 1000 -> 1000 -> 400 -> 100 kJ/mol/nm2 con riferimento min_c3.gro. Un sistema per GPU, in parallelo.
# Dopo s4 le repliche si separano (06_repliche.sh): fino a s4 il tronco e' comune e identico per le repliche.
#
# Uso (dalla cartella complesso_fab, GROMACS caricato):
#   setsid nohup bash 05_equil_s1_s4.sh "0 1" > 05.out 2>&1 < /dev/null &      # GPU 0 al complesso, GPU 1 al libero
# Variabili: GMX (default gmx), NT (thread OpenMP per run, default 40: come nel test di velocita')
# Riprende da solo: salta gli stadi con il .gro gia' scritto. Conta gli avvisi LINCS di ogni stadio in avvisi_lincs.txt e si ferma
# (quel sistema) se uno stadio ne ha; dopo s1 scrive controllo_s1.txt (05a, contatti di Asp40). IGNORA_LINCS=1 per proseguire dopo il controllo.
set -euo pipefail
QUI=$(cd "$(dirname "$0")" && pwd)
read -r -a G <<< "${1:-0 1}"
GMX=${GMX:-gmx}; NT=${NT:-40}

lancia() {
  local S=$1 GPU=$2 prev=min_c3
  cd "$QUI/$S"
  : > avvisi_lincs.txt
  for s in equil_s1 equil_s2 equil_s3 equil_s4; do
    if [ ! -f $s.gro ]; then
      if [ $s = equil_s1 ]; then T=""; else T="-t $prev.cpt"; fi
      $GMX grompp -f "$QUI/mdp/$s.mdp" -c $prev.gro $T -r min_c3.gro -p topol.top -n index.ndx -o $s.tpr -maxwarn 0 > grompp_$s.log 2>&1 \
        || { echo "$S: grompp $s fallito, vedi $S/grompp_$s.log"; return 1; }
      $GMX mdrun -deffnm $s -ntmpi 1 -ntomp $NT -gpu_id $GPU > mdrun_$s.out 2>&1 \
        || { echo "$S: mdrun $s fallito, vedi $S/mdrun_$s.out"; return 1; }
    fi
    n=$(grep -c 'LINCS WARNING' mdrun_$s.out || true)
    echo "$s avvisi_LINCS $n $(grep Performance: $s.log | awk '{print $2" ns/giorno"}')" | tee -a avvisi_lincs.txt
    if [ $s = equil_s1 ]; then
      # revisore 1 nov: Fmax di min_c3 sul CG di Asp40 di MrkA -> contatti di Asp40 prima e dopo s1 (05a)
      python3 "$QUI/05a_controlla_s1.py" min_c3.gro equil_s1.gro > controllo_s1.txt 2>&1 || echo "$S: 05a non riuscito, vedi $S/controllo_s1.txt"
      tail -1 controllo_s1.txt
    fi
    # revisore 1 nov: con avvisi LINCS ci si ferma e si guarda prima dello stadio successivo (dopo il controllo: IGNORA_LINCS=1)
    if [ "$n" != 0 ] && [ "${IGNORA_LINCS:-0}" != 1 ]; then
      echo "$S: FERMO dopo $s: $n avvisi LINCS in $S/mdrun_$s.out, da guardare prima di proseguire"; return 1
    fi
    prev=$s
  done
  echo "$S: equilibrazione s1-s4 finita (equil_s4.gro, equil_s4.cpt)"
}

lancia complesso "${G[0]}" > "$QUI/05_complesso.out" 2>&1 & P1=$!
lancia libero "${G[1]}" > "$QUI/05_libero.out" 2>&1 & P2=$!
RC=0; wait $P1 || RC=1; wait $P2 || RC=1
cat "$QUI/05_complesso.out" "$QUI/05_libero.out"
exit $RC
