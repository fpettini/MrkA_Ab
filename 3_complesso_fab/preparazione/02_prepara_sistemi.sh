#!/usr/bin/env bash
# Prepara i sistemi complesso (MrkA + Fab intero + 7 acque), libero (sola catena A) e fv (MrkA + VL 1-107 + VH 1-117 + 7 acque,
# solo per il test di velocita'). Stesso setup delle apo X-ray: amber14sb, TIP3P, terminali carichi, dodecaedro -d 1.3,
# KCl 0.2 M (francesco, 31 ott), minimizzazioni min_c1-c3 (mdp copiati da xray/). Protonazione: default di pdb2gmx (Lys94 carica); His neutre (PROPKA
# le da' tutte sotto 7) con tautomeri fissati in 02a (ISTIDINE): l'assegnazione automatica di pdb2gmx faceva H204 = HIP. Ponti S-S: automatici di pdb2gmx, contati e confrontati con le Cys presenti.
#
# Acque cristallografiche: le 7 a ponte MrkA-Fab (complesso, fv) + quelle di acque_interne.txt (01b, in tutti i sistemi).
#
# Uso (dalla cartella complesso_fab, GROMACS caricato):
#   bash 02_prepara_sistemi.sh mrka-fab_completo.pdb [sistemi] [acque_interne.txt]
#   es. bash 02_prepara_sistemi.sh mrka-fab_completo.pdb "complesso libero fv" acque_interne.txt
#   senza il terzo argomento: nessuna acqua interna, solo le 7 a ponte nel complesso
# Variabili: GMX (default gmx), NT (thread per le minimizzazioni, default 40),
#   FF (default amber14sb), WATER (default tip3p), SOLV (default spc216.gro; per OPC a 4 siti usare tip4p.gro)
#   es. FF=amber19sb WATER=opc SOLV=tip4p.gro bash 02_prepara_sistemi.sh ...  (nomi da verificare in $GMXDATA/top)
set -euo pipefail
QUI=$(cd "$(dirname "$0")" && pwd)
IN=$(readlink -f "$1"); SISTEMI=${2:-"complesso libero fv"}
INTERNE=${3:+$(readlink -f "$3")}
GMX=${GMX:-gmx}; NT=${NT:-40}; CONC=${CONC:-0.2}; PNAME=${PNAME:-K}; NNAME=${NNAME:-CL}; FF=${FF:-amber14sb}; WATER=${WATER:-tip3p}; SOLV=${SOLV:-spc216.gro}
command -v $GMX > /dev/null || { echo "$GMX non trovato: carica GROMACS (GMXRC)"; exit 1; }
for f in ions min_c1 min_c2 min_c3; do [ -f "$QUI/mdp/$f.mdp" ] || { echo "manca mdp/$f.mdp"; exit 1; }; done
echo "=== $(date "+%F %T") sistemi: $SISTEMI" | tee -a "$QUI/provenienza_02.txt"
sha256sum "$IN" "$0" "$QUI/02a_aiuto.py" $INTERNE | tee -a "$QUI/provenienza_02.txt"
$GMX --version 2>/dev/null | grep -E "GROMACS version|Precision|GPU support" | tee -a "$QUI/provenienza_02.txt"
echo "FF=$FF WATER=$WATER SOLV=$SOLV sale=${PNAME}${NNAME} $CONC M" | tee -a "$QUI/provenienza_02.txt"
python3 "$QUI/02a_aiuto.py" dividi "$IN" "$QUI" $INTERNE

for S in $SISTEMI; do
  cd "$QUI/$S"
  echo "=== $S"
  $GMX pdb2gmx -f proteina.pdb -o proteina.gro -p topol.top -i posre.itp \
       -ff $FF -water $WATER -ignh -merge all > pdb2gmx.log 2>&1 || { echo "pdb2gmx fallito: $S/pdb2gmx.log"; exit 1; }
  # ponti S-S: devono essere uno per ogni coppia di Cys del modello
  NSS=$(grep -c "Linking CYS-" pdb2gmx.log || true); ATT=$(cat ss_attesi.txt)
  echo "ponti S-S: $NSS (attesi $ATT)"; grep "Linking CYS-" pdb2gmx.log || true
  [ "$NSS" = "$ATT" ] || { echo "numero di ponti S-S diverso dall'atteso"; exit 1; }
  # istidine: stati fissati in proteina.pdb da 02a (ISTIDINE, tutte neutre); pdb2gmx non deve riassegnarle ne' fare HIP
  if grep -q "Will use HIS" pdb2gmx.log; then echo "pdb2gmx ha riassegnato delle istidine: $S/pdb2gmx.log"; exit 1; fi
  NHIP=$(awk '/^\[ *atoms/{a=1;next} /^\[/{a=0} a && $4 ~ /HIP$/{print $3}' topol.top | sort -u | wc -l)
  [ "$NHIP" = 0 ] || { echo "istidine cariche (HIP) in topologia: $NHIP"; exit 1; }
  awk '/^\[ *atoms/{a=1;next} /^\[/{a=0} a && $4 ~ /HI[DE]$/{print $3, $4}' topol.top | sort -u | tr '\n' ' ' | sed 's/^/istidine: /;s/$/\n/' | tee -a "$QUI/provenienza_02.txt"
  # restraint con costante da macro, come nelle apo
  sed -i -E 's/^(\s+[0-9]+\s+1)\s+1000\s+1000\s+1000/\1  POSRES_FC  POSRES_FC  POSRES_FC/' posre.itp
  grep -q POSRES_FC posre.itp || { echo "posre.itp non convertito"; exit 1; }
  NCRIS=0
  if [ -f acque.txt ]; then
    python3 "$QUI/02a_aiuto.py" acque proteina.gro acque.txt topol.top conf.gro $WATER; NCRIS=$(wc -l < acque.txt)
  else
    cp proteina.gro conf.gro
  fi
  $GMX editconf -f conf.gro -o box.gro -bt dodecahedron -d 1.3 -c > editconf.log 2>&1
  $GMX solvate -cp box.gro -cs $SOLV -p topol.top -o solv.gro > solvate.log 2>&1
  # -maxwarn 1: l'unico avviso atteso e' la carica netta prima della neutralizzazione (controllato sotto)
  $GMX grompp -f "$QUI/mdp/ions.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 1 > grompp_ions.log 2>&1 \
    || { echo "grompp ions fallito: $S/grompp_ions.log"; exit 1; }
  # (niente grep -q in pipe: con pipefail il SIGPIPE del primo grep fa fallire il controllo)
  NAVV=$(grep -c "^WARNING" grompp_ions.log || true)
  NCAR=$(grep -A3 "^WARNING" grompp_ions.log | grep -cE "non-zero total charge|system with net charge" || true)
  if [ "$NAVV" -gt 1 ] || [ "$NAVV" != "$NCAR" ]; then
    echo "grompp ions: avvisi $NAVV, di cui sulla carica netta $NCAR; vedi $S/grompp_ions.log"; grep -A4 "^WARNING" grompp_ions.log; exit 1
  fi
  # genion sul gruppo SOL; le acque cristallografiche sono SOL anche loro: se genion ne sostituisce una si riprova con un altro seme
  cp topol.top topol_prima_ioni.top
  for SEME in $(seq 1 20); do
    cp topol_prima_ioni.top topol.top
    echo SOL | $GMX genion -s ions.tpr -o ions.gro -p topol.top -pname $PNAME -nname $NNAME -neutral -conc $CONC -seed $SEME > genion.log 2>&1 \
      || { echo "genion fallito: $S/genion.log"; exit 1; }
    python3 "$QUI/02a_aiuto.py" controlla solv.gro ions.gro $NCRIS && break
    [ "$SEME" = 20 ] && { echo "genion tocca sempre un'acqua cristallografica"; exit 1; }
  done
  echo "genion: seme $SEME" >> "$QUI/provenienza_02.txt"
  printf "q\n" | $GMX make_ndx -f ions.gro -o index.ndx > make_ndx.log 2>&1
  for g in Protein Water_and_ions; do grep -q "\[ $g \]" index.ndx || { echo "gruppo $g assente in index.ndx"; exit 1; }; done
  # make_ndx lascia fuori dal gruppo Protein il C-terminale CCYX (Cys214 della catena leggera): gruppo Protein_Other (= Protein | Other) per tc-grps
  python3 "$QUI/02a_aiuto.py" gruppo ions.gro index.ndx
  # minimizzazioni come nelle apo (restraint 1000 su atomi pesanti, steep -> steep -> cg)
  prev=ions
  for m in min_c1 min_c2 min_c3; do
    $GMX grompp -f "$QUI/mdp/$m.mdp" -c $prev.gro -r ions.gro -p topol.top -n index.ndx -o $m.tpr -maxwarn 0 > grompp_$m.log 2>&1 \
      || { echo "grompp $m fallito: $S/grompp_$m.log"; exit 1; }
    $GMX mdrun -deffnm $m -ntmpi 1 -ntomp $NT > mdrun_$m.log 2>&1
    prev=$m
  done
  NAT=$(sed -n 2p min_c3.gro | tr -d ' ')
  CARICA=$(grep -m1 "System has non-zero total charge" grompp_ions.log | awk '{print $NF}' || true)
  printf "%s: acque cristallografiche %s, atomi %s, carica prima degli ioni %s, %s\n" "$S" "$NCRIS" "$NAT" "${CARICA:-0}" "$(grep -E "^($PNAME|$NNAME) " topol.top | tr '\n' ' ')" | tee -a "$QUI/riepilogo_02.txt"
  grep -E "Maximum force|Potential Energy" min_c3.log | tail -2 | tee -a "$QUI/riepilogo_02.txt"
done
echo "pronto: $SISTEMI (vedi riepilogo_02.txt)"
