#!/usr/bin/env bash
# Lancia (o riprende) un run REST2 preparato da 04_prepara_rest2.sh, tutte le repliche in un solo mpirun.
#
# Uso:  bash rest2/05_lancia_rest2.sh <cartella_run> [gpu_id] [thread_per_replica]
#   es. nohup bash rest2/05_lancia_rest2.sh xray/rest2_8x480 > /dev/null 2>&1 &
# Default: tutte le GPU fisiche (nvidia-smi -L) e thread = core fisici / n_repliche, cosi' ogni replica
# ha core suoi senza hyperthreading; gpu_id "cpu" per girare senza GPU.
# Per decidere replica per replica: gputasks=<2 cifre per replica, PP e PME, in ordine di replica>
#   es. gputasks=3333334444556677 = r0-r2 sulla GPU 3, r3-r4 sulla 4, r5 sulla 5, r6 sulla 6, r7 sulla 7.
# Se trova rest2.cpt riprende dal checkpoint. Durante il run mdrun.out non cambia (mdrun scrive solo all'avvio):
# l'avanzamento si vede con  python3 rest2/06_controlla_scambi.py <cartella_run>  (ultimo passo e ns).
# -replex 500: un tentativo di scambio ogni 1 ps; -hrex: scambio fra Hamiltoniane diverse (patch PLUMED).
# -notunepme: su GPU mdrun ritocca cutoff Coulomb e griglia PME replica per replica, anche con il replica
# exchange, e le Hamiltoniane smetterebbero di essere quelle della scala.
# setsid: mpirun di Open MPI ignora nohup (rimette un suo gestore di SIGHUP) e quando il terminale si chiude
# (logout, ssh caduto) ferma tutte le repliche con TERM. In una sessione sua non riceve il SIGHUP del terminale.
# Per lo stesso motivo Ctrl+C non ferma il run: per fermarlo  pkill -TERM -u $USER gmx_mpi
set -euo pipefail
# controlli: ambiente giusto e nessun altro run REST2 gia' in corso (due run nelle stesse cartelle si rovinano a vicenda)
command -v gmx_mpi > /dev/null || { echo "gmx_mpi non trovato: prima  source rest2/attiva_gmx2024.sh"; exit 1; }
mpirun --version 2>&1 | grep -q "Open MPI" || { echo "mpirun non e' quello di Open MPI dell'ambiente: prima  source rest2/attiva_gmx2024.sh"; exit 1; }
if pgrep -u "$(id -u)" -x gmx_mpi > /dev/null; then
  echo "C'e' gia' un gmx_mpi in esecuzione (pgrep -a gmx_mpi). Per fermarlo: pkill -TERM -u \$USER gmx_mpi"; exit 1
fi
RUN=$(cd "$1" && pwd)
N=$(wc -l < "$RUN/scala.txt")
GPU=${2:-$(seq -s '' 0 $(( $(nvidia-smi -L | grep -c '^GPU') - 1 )))}
NCORE=$(lscpu -p=core,socket | grep -v '^#' | sort -u | wc -l)
NT=${3:-$(( NCORE / N ))}
cd "$RUN"
DIRS=$(seq -f "r%g" 0 $(( N - 1 )) | tr '\n' ' ')
REPLEX=${REPLEX:-500}
case "$GPU" in
  cpu)        DEV="-nb cpu -pme cpu" ;;
  gputasks=*) M=${GPU#gputasks=}
              [ ${#M} -eq $(( 2 * N )) ] || { echo "gputasks: servono $(( 2 * N )) cifre, ne ho $(( ${#M} ))"; exit 1; }
              DEV="-nb gpu -pme gpu -gputasks $M" ;;
  *)          DEV="-nb gpu -pme gpu -gpu_id $GPU" ;;
esac
CPI=""
[ -f r0/rest2.cpt ] && CPI="-cpi rest2.cpt"
setsid -w mpirun -np "$N" --bind-to none gmx_mpi mdrun -multidir $DIRS -deffnm rest2 $CPI \
  -plumed plumed.dat -replex "$REPLEX" -hrex -dlb no -notunepme -ntomp "$NT" -pin on $DEV < /dev/null >> mdrun.out 2>&1
