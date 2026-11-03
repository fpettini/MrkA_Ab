#!/usr/bin/env bash
# Prepara un run REST2: una cartella per replica con topologia scalata, mdp e tpr.
#
# Uso (ambiente conda attivo e GMXRC della 2024.3 caricato, vedi COMANDI.md):
#   bash rest2/04_prepara_rest2.sh <cartella_sistema> <n_repliche> <Tmax_K> <ns> <nome_run>
#   es. bash rest2/04_prepara_rest2.sh xray 8 480 200 rest2_8x480
#
# <cartella_sistema> deve contenere equil_s4.gro e processed.top (topologia espansa, grompp -pp).
# Le temperature efficaci vanno da 310.15 K a Tmax in progressione geometrica; lambda = 310.15/T.
# Con Tmax = 310.15 tutte le repliche hanno lambda = 1: e' il test in cui gli scambi vanno accettati tutti.
set -euo pipefail
QUI=$(cd "$(dirname "$0")" && pwd)
[ $# -eq 5 ] || { sed -n 4,6p "$0"; exit 1; }
SIS=$(cd "$1" && pwd); N=$2; TMAX=$3; NS=$4; NOME=$5
T0=310.15
for f in equil_s4.gro processed.top; do
  [ -f "$SIS/$f" ] || { echo "manca $SIS/$f"; exit 1; }
done
for c in gmx_mpi plumed; do
  command -v $c > /dev/null || { echo "$c non trovato: attiva l'ambiente conda e carica GMXRC"; exit 1; }
done
RUN=$SIS/$NOME
[ -e "$RUN" ] && { echo "$RUN esiste gia': scegli un altro nome"; exit 1; }
mkdir -p "$RUN"

python3 "$QUI/03_marca_caldi.py" "$SIS/processed.top" "$RUN/processed_caldi.top" | tee "$RUN/caldi.txt"
NSTEPS=$(python3 -c "print(round($NS * 500000))")
python3 - "$N" "$T0" "$TMAX" > "$RUN/scala.txt" <<'PY'
import sys
n, t0, tmax = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
for i in range(n):
    t = t0 * (tmax / t0) ** (i / (n - 1))
    print(f"{i} {t:.2f} {t0 / t:.6f}")
PY

while read -r i T LAM; do
  d=$RUN/r$i
  mkdir -p "$d"
  plumed partial_tempering "$LAM" < "$RUN/processed_caldi.top" > "$d/topol.top"
  sed -e "s/@NSTEPS@/$NSTEPS/" -e "s/@SEED@/$(( 1000 + 7919 * (i + 1) ))/" "$QUI/rest2.mdp.in" > "$d/rest2.mdp"
  # l'unico avviso ammesso e' la carica netta del box: e' inevitabile, perche' le cariche calde
  # (in tutto -4 e) sono moltiplicate per sqrt(lambda); qualsiasi altro avviso ferma lo script
  if ! gmx_mpi grompp -f "$d/rest2.mdp" -c "$SIS/equil_s4.gro" -p "$d/topol.top" -o "$d/rest2.tpr" \
       -po "$d/mdout.mdp" -maxwarn 1 > "$d/grompp.log" 2>&1; then
    echo "grompp fallito in $d: vedi $d/grompp.log"; exit 1
  fi
  if grep -q "^WARNING" "$d/grompp.log" && ! grep -A3 "^WARNING" "$d/grompp.log" | grep -q "Ewald electrostatics in a system with net charge"; then
    echo "avviso inatteso di grompp in $d: vedi $d/grompp.log"; exit 1
  fi
  : > "$d/plumed.dat"
  printf "replica %2d  T_eff %7.2f K  lambda %s\n" "$i" "$T" "$LAM"
done < "$RUN/scala.txt"
echo "pronto: $RUN ($N repliche, $NS ns, $NSTEPS passi)"
