#!/usr/bin/env bash
# Test della deriva dell'energia conservata (revisore, 2 nov 2026), solo log/energie, niente traiettoria.
# Stesso punto di partenza (equil_s4.gro), nessun restraint, velocita' generate con lo stesso seme, mdp di produzione
# (mdp/md.mdp.in), NS ns per variante:
#   rif_2fs        2 fs, HMR 1 (riferimento)
#   a_4fs          4 fs, HMR 3, come in produzione (solo controllo: il criterio 2 e' gia' fallito in s5)
#   b_4fs_lincs8   come a) con lincs-order 8 e lincs-iter 2
#   c_4fs_vbt      come a) con verlet-buffer-tolerance 0.0005
#   d_3fs          3 fs, HMR 3, il resto come in produzione
# Criteri fissati dal revisore prima del test, entrambi per ogni variante:
#   1. |deriva| al massimo 3 volte quella di rif_2fs;  2. T media di Protein_Other e di Water_and_ions entro 298 +- 1 K.
#   3. (revisore, 11:13) per ciascun gruppo: T della seconda meta' entro 298 +- 1 K e (seconda - prima meta') < 0.5 K.
# Decisione (revisore): fra le varianti che passano in ENTRAMBI i sistemi vale la piu' veloce; si applica uguale ai due
# sistemi e si riparte da s5. Se non passa nessuna, entrambi a 2 fs.
# Uso (dalla cartella complesso_fab, con 06 fermo):
#   bash 07_test_deriva.sh <libero|complesso> "<gpu> [gpu ...]" [ns=2] [ntomp=60]
#   con piu' GPU le varianti si dividono fra le GPU e girano in parallelo; rilanciato, salta le varianti gia' finite.
set -euo pipefail
QUI=$(cd "$(dirname "$0")" && pwd)
S=$1; read -r -a GPUS <<< "$2"; NS=${3:-2}; NTOMP=${4:-60}; GMX=${GMX:-gmx}
case $S in libero|complesso) ;; *) echo "sistema: libero o complesso"; exit 1 ;; esac
[ ${#GPUS[@]} -gt 0 ] || { echo "manca la GPU"; exit 1; }
[ -f "$QUI/$S/equil_s4.gro" ] || { echo "manca $S/equil_s4.gro"; exit 1; }
D="$QUI/$S/test_deriva"; mkdir -p "$D"; cd "$D"
VARIANTI="rif_2fs a_4fs b_4fs_lincs8 c_4fs_vbt d_3fs"
for V in $VARIANTI; do
  case $V in
    rif_2fs) DTPS=0.002; HMR=1; NSTE=10000; PASSI=$(( NS * 500000 )) ;;
    d_3fs)   DTPS=0.003; HMR=3; NSTE=5000;  PASSI=$(( NS * 1000000 / 3 )) ;;
    *)       DTPS=0.004; HMR=3; NSTE=5000;  PASSI=$(( NS * 250000 )) ;;
  esac
  M=deriva_$V.mdp
  sed -e 's|${CONT}|no|g' -e 's|${GVEL}|yes|g' -e 's|${GSEED}|12345|g' -e 's|${LSEED}|12345|g' -e "s|\${DT}|$DTPS|g" \
      -e "s|\${NSTEPS}|$PASSI|g" -e "s|\${HMR}|$HMR|g" -e 's|${NSTX}|0|g' -e "s|\${NSTE}|$NSTE|g" "$QUI/mdp/md.mdp.in" > $M
  case $V in
    b_4fs_lincs8) sed -i -E 's/^(lincs-order +)= *6 *$/\1= 8/' $M; echo "lincs-iter               = 2" >> $M
                  grep -qE '^lincs-order += 8$' $M || { echo "lincs-order non sostituito in $M"; exit 1; } ;;
    c_4fs_vbt)    grep -q '^verlet-buffer-tolerance' $M && { echo "verlet-buffer-tolerance gia' nel template"; exit 1; }
                  echo "verlet-buffer-tolerance  = 0.0005" >> $M ;;
  esac
  ! grep -q '\${' $M || { echo "segnaposto non sostituiti in $M"; exit 1; }
done
finita() { [ -f deriva_$1.log ] && grep -q 'Conserved energy drift' deriva_$1.log; }
corri() {   # variante gpu
  local V=$1 g=$2
  finita $V && return 0
  $GMX grompp -f deriva_$V.mdp -c ../equil_s4.gro -p ../topol.top -n ../index.ndx -o deriva_$V.tpr -maxwarn 0 > grompp_$V.log 2>&1 \
    || { echo "grompp $V fallito: $D/grompp_$V.log"; return 1; }
  echo "$(date '+%F %T') $S: $NS ns, variante $V, GPU $g"
  rm -f deriva_$V.log deriva_$V.edr deriva_$V.cpt   # un test interrotto si rifa' da capo
  $GMX mdrun -deffnm deriva_$V -ntmpi 1 -ntomp $NTOMP -gpu_id $g > mdrun_$V.out 2>&1 || { echo "mdrun $V fallito: $D/mdrun_$V.out"; return 1; }
}
declare -A LISTA; i=0
for V in $VARIANTI; do g=${GPUS[$(( i % ${#GPUS[@]} ))]}; LISTA[$g]="${LISTA[$g]:-}$V "; i=$(( i + 1 )); done
PIDS=()
for g in "${!LISTA[@]}"; do ( for V in ${LISTA[$g]}; do corri $V $g || exit 1; done ) & PIDS+=($!); done
err=0; for p in "${PIDS[@]}"; do wait $p || err=1; done
[ $err = 0 ] || { echo "almeno una variante e' fallita: vedi sopra"; exit 1; }

deriva() { awk '/Conserved energy drift/{d=$4} END{print d}' "$1"; }
tgruppi() {   # medie finali delle T dei gruppi dal log: GROMACS scrive le intestazioni attaccate ("T-Protein_OtherT-Water_and_ions"),
             # i valori sono nell'ordine di tc-grps (Protein_Other Water_and_ions)
  awk '/A V E R A G E S/{f=1} f && /T-Protein_Other/{getline; printf "%s %s\n", $1, $2; exit}' "$1"
}
meta() {      # variante: T di proteina e acqua+ioni nella prima e nella seconda meta' (gmx energy sull'edr): "p1 p2 w1 w2"
  printf "T-Protein_Other\nT-Water_and_ions\n\n" | $GMX energy -f deriva_$1.edr -o T_$1.xvg > energy_$1.log 2>&1 || { echo "- - - -"; return; }
  awk 'BEGIN{n=0} /^@ s[0-9]+ legend/{k=substr($2,2)+2; if ($0 ~ /T-Protein_Other/) cp=k; if ($0 ~ /T-Water_and_ions/) cw=k}
       /^[#@]/{next} $1 > 0 {t[n]=$1; p[n]=$cp; w[n]=$cw; n++}   # t=0 escluso: T non fisiche al passo 0 (614-631 K a 4 fs)
       END{if (!cp || !cw || n < 4) {print "- - - -"; exit}; h=t[n-1]/2
           for(i=0;i<n;i++) if (t[i] < h) {p1+=p[i]; w1+=w[i]; na++} else {p2+=p[i]; w2+=w[i]; nb++}
           printf "%.2f %.2f %.2f %.2f\n", p1/na, p2/nb, w1/na, w2/nb}' T_$1.xvg
}
R=$(deriva deriva_rif_2fs.log)
{
  echo "$S, $NS ns per variante, partenza equil_s4.gro, seme 12345."
  echo "Criteri: 1 |deriva| <= 3x rif_2fs; 2 T media di entrambi i gruppi in 298 +- 1 K; 3 per ogni gruppo 2a meta' in 298 +- 1 K e 2a - 1a < 0.5 K"
  printf "%-14s %11s %6s %8s %8s %15s %15s %6s %8s  %s\n" variante deriva "x rif" "T prot" "T acq" "prot 1a/2a" "acq 1a/2a" LINCS "ns/g" esito
  for V in $VARIANTI; do
    d=$(deriva deriva_$V.log); read -r tp tw <<< "$(tgruppi deriva_$V.log)"; read -r p1 p2 w1 w2 <<< "$(meta $V)"
    w=$(grep -c 'LINCS WARNING' mdrun_$V.out || true); p=$(awk '/^Performance:/{print $2}' deriva_$V.log)
    awk -v V=$V -v d="$d" -v r="$R" -v tp="$tp" -v tw="$tw" -v p1="$p1" -v p2="$p2" -v w1="$w1" -v w2="$w2" -v w="$w" -v p="$p" 'function a(x){return x<0?-x:x}
      BEGIN{x = a(d)/a(r); c1 = (V=="rif_2fs") || x <= 3; c2 = a(tp-298) <= 1 && a(tw-298) <= 1
            c3 = (p2 != "-") && a(p2-298) <= 1 && (p2-p1) < 0.5 && a(w2-298) <= 1 && (w2-w1) < 0.5
            m = (c1 ? "" : " deriva") (c2 ? "" : " T") (c3 ? "" : " T-nel-tempo")
            e = (V=="rif_2fs") ? ((c2 && c3) ? "riferimento, T ok" : "riferimento, FUORI:" m) : ((c1 && c2 && c3) ? "PASSA" : "NON PASSA:" m)
            printf "%-14s %11.3e %6.1f %8.2f %8.2f %7s/%-7s %7s/%-7s %6d %8s  %s\n", V, d, x, tp, tw, p1, p2, w1, w2, w, p, e}'
  done
} | tee risultato.txt
