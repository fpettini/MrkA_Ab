#!/bin/bash
# Protocollo rivisto: minimizzazione (doppia precisione) + 5 stadi di equilibrazione + grompp della produzione.
# Restraint (topol.top): backbone = proteina + N CA C O di NLN/OLT; zuccheri = solo atomi pesanti.
# Le catene laterali di NLN/OLT non sono mai vincolate. Nessun restraint diedro.
set -e

#source /home/f.pettini/gmx_db/bin/GMXRC
#GMX_D=gmx_d        # doppia precisione per le minimizzazioni
#NT_D=4             # thread CPU per gmx_d (con 844k atomi conviene alzarlo)
#MAXWARN=0          # validato con GROMACS 2023: nessun WARNING, solo una NOTE (COM + posres)
#
#### --- minimizzazioni: backbone 1000 -> 400 -> 100, zuccheri liberi ---
#$GMX_D grompp -f min_c1.mdp -c ions.gro -r ions.gro -p topol.top -n index.ndx -o min_c1.tpr -maxwarn $MAXWARN
#$GMX_D mdrun -v -deffnm min_c1 -nt $NT_D
#$GMX_D grompp -f min_c2.mdp -c min_c1.gro -r ions.gro -p topol.top -n index.ndx -o min_c2.tpr -maxwarn $MAXWARN
#$GMX_D mdrun -v -deffnm min_c2 -nt $NT_D
#$GMX_D grompp -f min_c3.mdp -c min_c2.gro -r ions.gro -p topol.top -n index.ndx -o min_c3.tpr -maxwarn $MAXWARN
#$GMX_D mdrun -v -deffnm min_c3 -nt $NT_D
#
#
#source /usr/local/gromacs/bin/GMXRC
#GMX=gmx
#NT=40
#GPU=6
#MAXWARN=0
##REF=min_c3.gro
REF=ions.gro
#
## --- equilibrazione (riferimento dei restraint: min_c3.gro; stato passato con -t .cpt) ---
#$GMX grompp -f equil_s1.mdp -c min_c3.gro -r $REF -p topol.top -n index.ndx -o equil_s1.tpr -maxwarn $MAXWARN
#$GMX mdrun -v -deffnm equil_s1 -nt $NT -gpu_id $GPU
#prev=equil_s1
#for s in equil_s2 equil_s3 equil_s4; do
#  $GMX grompp -f $s.mdp -c $prev.gro -t $prev.cpt -r $REF -p topol.top -n index.ndx -o $s.tpr -maxwarn $MAXWARN
#  $GMX mdrun -v -deffnm $s -nt $NT -gpu_id $GPU
#  prev=$s
#done

# =============================================================================
#  REPLICHE
#
#  Il tronco comune finisce a equil_s4, dove il backbone e' ancora vincolato a
#  100 kJ/mol/nm2: la conformazione proteica e' quindi identica per tutte le
#  repliche e non introduciamo differenze di modello. La diramazione avviene in
#  equil_s5, 20 ns liberi in cui le velocita' sono RIGENERATE con un seme diverso
#  per replica: sono quei 20 ns a decorrelare catene laterali, glicani e solvente.
#  In produzione le velocita' NON si rigenerano (si riprende con -t .cpt), ma
#  ld-seed resta diverso perche' v-rescale e' stocastico per tutti i 500 ns.
#
#  Semi espliciti e non -1: cosi' una replica e' riproducibile e i valori si
#  possono dichiarare in Methods.
# =============================================================================
source /usr/local/gromacs/bin/GMXRC
GMX=gmx
NRUN=3                       # numero di repliche
GPUS=(0 1 2)               # una GPU per replica (si ripete se NRUN > #GPU)
#GPUS=6
NT_REP=80                    # thread per replica: NT totali / repliche in parallelo
PARALLEL=0                   # 1 = le repliche girano insieme, 0 = una dopo l'altra
RUN_PROD=1                   # 0 = prepara solo md.tpr senza lanciare i 500 ns
TEMPLATES=templates
MAXWARN=0
LENGTH=500000

#mdp_from_template() {        # $1 template  $2 uscita  poi CONT GVEL GSEED LSEED
#  sed -e "s|\${CONT}|$3|g" -e "s|\${GVEL}|$4|g" \
#      -e "s|\${GSEED}|$5|g" -e "s|\${LSEED}|$6|g" "$1" > "$2"
#  if grep -q '\${' "$2"; then
#    echo "ERRORE: segnaposto non sostituiti in $2"; grep -n '\${' "$2"; exit 1
#  fi
#}

run_replica() {
  local r=$1
  local dir=run$r
  local gpu=${GPUS[$(( (r-1) % ${#GPUS[@]} ))]}
  local seed=$(( r * 11111 ))          # run1=11111 run2=22222 ...

  #mkdir -p "$dir"
  #
  ## --- equil_s5: diramazione, velocita' rigenerate dal seme ---
  #mdp_from_template "$TEMPLATES/equil_s5.mdp.in" "$dir/equil_s5.mdp" no yes $seed $seed
  #$GMX grompp -f "$dir/equil_s5.mdp" -c equil_s4.gro -r $REF \
  #            -p topol.top -n index.ndx -o "$dir/equil_s5.tpr" -maxwarn $MAXWARN
  #$GMX mdrun -v -deffnm "$dir/equil_s5" -nt $NT_REP -gpu_id $gpu
  #
  ## --- produzione: riprende dal cpt, gen-vel=no, conta solo ld-seed ---
  #mdp_from_template "$TEMPLATES/md.mdp.in" "$dir/run$r.mdp" yes no $seed $seed
  #$GMX grompp -f "$dir/run$r.mdp" -c "$dir/equil_s5.gro" -t "$dir/equil_s5.cpt" \
  #            -p topol.top -n index.ndx -o "$dir/run$r.tpr" -maxwarn $MAXWARN
  #if [ "$RUN_PROD" = "1" ]; then
  #  $GMX mdrun -v -deffnm "$dir/run$r" -nt $NT_REP -gpu_id $gpu
  #fi
  #echo "replica $r completata (GPU $gpu, seme $seed)"
  $GMX convert-tpr -s "$dir/run$r.cpt" -until $LENGTH -o "$dir/run$r.tpr"
  $GMX mdrun -v -deffnm "$dir/run$r" -nt $NT_REP -gpu_id $gpu -s "$dir/run$r.tpr" -cpi "$dir/run$r.cpt" -append
}

for r in $(seq 1 $NRUN); do
  if [ "$PARALLEL" = "1" ]; then run_replica $r & else run_replica $r; fi
done
wait

echo
echo "Fatto. Una cartella per replica:"
for r in $(seq 1 $NRUN); do
  printf "  run%s: seme %s  " $r $(( r * 11111 ))
  [ -f run$r/run$r.tpr ] && printf "run$r.tpr pronto  "
  [ -f run$r/run$r.xtc ] && printf "produzione conclusa"
  echo
done
