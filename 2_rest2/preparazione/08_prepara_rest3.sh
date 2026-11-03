#!/usr/bin/env bash
# Prepara un run REST3: come 04_prepara_rest2.sh (stessi atomi caldi, stessa scala di lambda, stesso mdp),
# in piu' il vdW fra atomi caldi e acqua/ioni di ogni replica e' moltiplicato per kappa (07_kappa_rest3.py).
#
# Uso (ambiente attivo: source rest2/attiva_gmx2024.sh):
#   bash rest2/08_prepara_rest3.sh <cartella_sistema> <n_repliche> <Tmax_K> <ns> <kappa_max> <nome_run>
#   es. bash rest2/08_prepara_rest3.sh rest2/xray 16 480 200 1.05 rest3_16x480_k105
#
# Scala di kappa: la curva kappa(lambda) del tutorial REST3 (8 repliche 298-450 K; coincide con la formula
# dell'articolo per 16 repliche 298-500 K), riportata in proporzione perche' la replica piu' calda abbia
# kappa_max. Con 16 repliche fino a 480 K la curva pubblicata da' 1.049 all'ultima replica: kappa_max = 1.05
# la riproduce entro 0.001. kappa_max = 1 da' di nuovo REST2.
# Le repliche con kappa = 1 (r0-r2 con 16 repliche fino a 480 K) hanno la topologia di REST2 senza modifiche, quindi r0 resta
# l'Hamiltoniana non scalata. scala.txt ha una colonna in piu' (kappa); 05 e 06 si usano come per REST2.
set -euo pipefail
QUI=$(cd "$(dirname "$0")" && pwd)
[ $# -eq 6 ] || { sed -n 4,6p "$0"; exit 1; }
SIS=$(cd "$1" && pwd); N=$2; TMAX=$3; NS=$4; KMAX=$5; NOME=$6
T0=310.15
for f in equil_s4.gro processed.top; do
  [ -f "$SIS/$f" ] || { echo "manca $SIS/$f"; exit 1; }
done
for c in gmx_mpi plumed; do
  command -v $c > /dev/null || { echo "$c non trovato: prima  source rest2/attiva_gmx2024.sh"; exit 1; }
done
RUN=$SIS/$NOME
[ -e "$RUN" ] && { echo "$RUN esiste gia': scegli un altro nome"; exit 1; }
mkdir -p "$RUN"

python3 "$QUI/03_marca_caldi.py" "$SIS/processed.top" "$RUN/processed_caldi.top" | tee "$RUN/caldi.txt"
NSTEPS=$(python3 -c "print(round($NS * 500000))")
python3 - "$N" "$T0" "$TMAX" "$KMAX" > "$RUN/scala.txt" <<'PY'
import sys
n, t0, tmax, kmax = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
# (lambda, kappa) del tutorial mdlab-um/REST3_tutorial, README.sh
TUT = [(1.000, 1.000), (0.943, 1.000), (0.889, 1.003), (0.838, 1.010),
       (0.790, 1.020), (0.745, 1.027), (0.702, 1.035), (0.662, 1.045)]
def k_tut(l):
    if l >= TUT[1][0]:
        return 1.0
    for (l1, k1), (l2, k2) in zip(TUT, TUT[1:]):
        if l >= l2:
            return k1 + (k2 - k1) * (l - l1) / (l2 - l1)
    (l1, k1), (l2, k2) = TUT[-2:]
    return k2 + (k2 - k1) * (l - l2) / (l2 - l1)     # oltre la tabella: pendenza dell'ultimo tratto
lam = [t0 / (t0 * (tmax / t0) ** (i / (n - 1))) for i in range(n)]
for i, l in enumerate(lam):
    k = 1 + (kmax - 1) * (k_tut(l) - 1) / (k_tut(lam[-1]) - 1)
    print(f"{i} {t0 / l:.2f} {l:.6f} {k:.4f}")
PY

while read -r i T LAM K; do
  d=$RUN/r$i
  mkdir -p "$d"
  if [ "$K" = "1.0000" ]; then
    plumed partial_tempering "$LAM" < "$RUN/processed_caldi.top" > "$d/topol.top"
  else
    plumed partial_tempering "$LAM" < "$RUN/processed_caldi.top" | python3 "$QUI/07_kappa_rest3.py" - "$LAM" "$K" > "$d/topol.top" 2> "$d/kappa.txt"
  fi
  sed -e "s/@NSTEPS@/$NSTEPS/" -e "s/@SEED@/$(( 1000 + 7919 * (i + 1) ))/" "$QUI/rest2.mdp.in" > "$d/rest2.mdp"
  # come in 04: l'unico avviso ammesso e' la carica netta del box (cariche calde x sqrt(lambda))
  if ! gmx_mpi grompp -f "$d/rest2.mdp" -c "$SIS/equil_s4.gro" -p "$d/topol.top" -o "$d/rest2.tpr" \
       -po "$d/mdout.mdp" -maxwarn 1 > "$d/grompp.log" 2>&1; then
    echo "grompp fallito in $d: vedi $d/grompp.log"; exit 1
  fi
  if grep -q "^WARNING" "$d/grompp.log" && ! grep -A3 "^WARNING" "$d/grompp.log" | grep -q "Ewald electrostatics in a system with net charge"; then
    echo "avviso inatteso di grompp in $d: vedi $d/grompp.log"; exit 1
  fi
  : > "$d/plumed.dat"
  printf "replica %2d  T_eff %7.2f K  lambda %s  kappa %s\n" "$i" "$T" "$LAM" "$K"
done < "$RUN/scala.txt"
echo "pronto: $RUN ($N repliche, $NS ns, $NSTEPS passi, kappa_max $KMAX)"
