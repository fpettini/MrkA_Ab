#!/usr/bin/env bash
# Controllo incrociato del port amber19sb.ff (richiesto dal revisore il 31 ott 2026 prima della produzione):
# energia a singolo punto della sola proteina (<sistema>/proteina.gro, uscita di pdb2gmx; sistema = libero o complesso) in un box di 20 nm, senza solvente.
# Si confrontano con OpenMM (04b_confronta_openmm.py) i termini che non dipendono dal cutoff: legami, angoli, diedri propri
# e impropri, CMAP, LJ-14 e Coulomb-14. Gli altri termini non legati NON sono confrontati (GROMACS Verlet non ha "nessun cutoff").
# Uso: bash 04_controllo_ff.sh [libero|complesso]   (dopo 02; default libero). Poi copia <sistema>/vuoto.log, vuoto.pdb e topol.top sul PC.
# Il complesso (31 ott, revisore) copre i residui che MrkA non ha (His, Met) e i terminali del Fab.
set -euo pipefail
QUI=$(cd "$(dirname "$0")" && pwd); GMX=${GMX:-gmx}
SIS=${1:-libero}
cd "$QUI/$SIS"
$GMX editconf -f proteina.gro -o vuoto.gro -box 20 20 20 -c > /dev/null 2>&1
# topologia della sola proteina: topol.top dopo 02 contiene anche acque cristallografiche, solvente e ioni (K, CL; NA per sicurezza)
awk '/^\[ *molecules *\]/{m=1} !(m && $1 ~ /^(SOL|NA|K|CL)$/)' topol.top > vuoto.top
$GMX grompp -f "$QUI/mdp/vuoto.mdp" -c vuoto.gro -p vuoto.top -o vuoto.tpr -maxwarn 0 > grompp_vuoto.log 2>&1 || { echo "grompp fallito: $SIS/grompp_vuoto.log"; exit 1; }
$GMX mdrun -s vuoto.tpr -rerun vuoto.gro -deffnm vuoto -nt 1 -nb cpu > /dev/null 2>&1
$GMX editconf -f vuoto.gro -o vuoto.pdb > /dev/null 2>&1
grep -A8 "Energies (kJ/mol)" vuoto.log | head -9
head -1 "${GMXDATA:-}/top/amber19sb.ff/forcefield.doc" 2>/dev/null; grep -i -m3 -E "sengupta|joung|li.*merz|opc" "${GMXDATA:-}/top/amber19sb.ff/ions_opc.itp" || true
echo "fatto: copia $SIS/vuoto.log, $SIS/vuoto.pdb e $SIS/topol.top sul PC (complesso_fab/$SIS/)"
