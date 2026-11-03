#!/usr/bin/env bash
# Allunga un run REST2 (o REST3) gia' finito: porta la durata di tutte le repliche a <ns_totali>.
#
# Uso (ambiente attivo: source rest2/attiva_gmx2024.sh):
#   bash rest2/09_allunga_rest2.sh <cartella_run> <ns_totali>
#   es. bash rest2/09_allunga_rest2.sh rest2/xray/rest2_16x480 400
# Poi si rilancia con 05, che riparte da rest2.cpt e aggiunge in coda a rest2.xtc, .log ed .edr:
#   nohup bash rest2/05_lancia_rest2.sh rest2/xray/rest2_16x480 > /dev/null 2>&1 &
# Il tpr vecchio resta come rest2_fino_<ns>ns.tpr. Si puo' ripetere con un <ns_totali> piu' grande.
# Collaudato su CPU (2 repliche, -hrex): riparte dal passo del checkpoint e il log continua.
set -euo pipefail
[ $# -eq 2 ] || { sed -n 4,6p "$0"; exit 1; }
command -v gmx_mpi > /dev/null || { echo "gmx_mpi non trovato: prima  source rest2/attiva_gmx2024.sh"; exit 1; }
if pgrep -u "$(id -u)" -x gmx_mpi > /dev/null; then
  echo "C'e' un gmx_mpi in esecuzione: il run va allungato da fermo (pgrep -a gmx_mpi)"; exit 1
fi
RUN=$(cd "$1" && pwd); NS=$2
N=$(wc -l < "$RUN/scala.txt")
for i in $(seq 0 $(( N - 1 ))); do
  d=$RUN/r$i
  [ -f "$d/rest2.cpt" ] || { echo "manca $d/rest2.cpt: senza checkpoint non si allunga"; exit 1; }
done
PS=$(python3 -c "print(round($NS * 1000))")
for i in $(seq 0 $(( N - 1 ))); do
  d=$RUN/r$i
  # awk legge tutto il dump: uscire prima chiuderebbe la pipe e con pipefail lo script si fermerebbe
  VECCHI=$(gmx_mpi dump -s "$d/rest2.tpr" 2> /dev/null | awk '$1 == "nsteps" && !f {print $3 / 500000; f = 1}')
  gmx_mpi convert-tpr -s "$d/rest2.tpr" -until "$PS" -o "$d/rest2_nuovo.tpr" > "$d/convert.log" 2>&1 \
    || { echo "convert-tpr fallito in $d: vedi $d/convert.log"; exit 1; }
  mv "$d/rest2.tpr" "$d/rest2_fino_${VECCHI}ns.tpr"
  mv "$d/rest2_nuovo.tpr" "$d/rest2.tpr"
  echo "r$i: da $VECCHI a $NS ns"
done
echo "pronto: rilancia con  nohup bash rest2/05_lancia_rest2.sh ${1%/} > /dev/null 2>&1 &"
