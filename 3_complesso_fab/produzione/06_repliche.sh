#!/usr/bin/env bash
# Repliche di complesso e libero dopo 05 (equil_s4 comune). Ogni replica r: equil_s5 (20 ns, velocita' RIGENERATE con
# seme r*11111) -> produzione (velocita' non rigenerate, ld-seed r*11111). Due code separate: le GPU di GPU_COMPLESSO
# lavorano solo sul complesso, quelle di GPU_LIBERO solo sul libero.
# TURNO A BLOCCHI (condizione del revisore, 31 ott 2026): la produzione avanza a tratti uguali per tutte le repliche,
# a rotazione: prima equil_s5 + i primi 10 ns di tutte, poi tutte fino a BLOCCO ns, poi fino a 2*BLOCCO, ... fino a NS_*.
# Cosi' le lunghezze restano simili e il "minimo fra tutte le repliche" (lunghezza confermativa) non spreca dati.
# Una GPU prende il primo tratto in coda la cui replica non e' gia' in corso su un'altra GPU.
#
# Regole del revisore (31 ott 2026), applicate in automatico:
#  - passo uguale per complesso e libero; 4 fs con HMR solo se ZERO avvisi LINCS in equil_s5 e nei primi 10 ns di
#    produzione di OGNI replica. Al primo avviso lo script scrive STOP_LINCS.txt e non parte piu' nessun lavoro nuovo:
#    decide francesco (la regola dice: tutti e due i sistemi tornano a 2 fs). Le repliche gia' in corso continuano.
#  - durata: analisi confermativa alla lunghezza raggiunta da TUTTE le repliche a fine disponibilita', max NS_*.
#  - (revisore, 2 nov 2026 16:00, dopo il test 07: 3 e 4 fs scartati, si usa DT=2) controllo della temperatura in
#    produzione, solo dall'edr: alla fine di ogni tratto (0-10, 10-250, 250-500, ... ns) T media di Protein_Other e di
#    Water_and_ions nel tratto, per replica, in run<r>/controllo_T.txt con l'esito:
#      AVVISO (il sistema prosegue; va elencato nei metodi coi valori): proteina o acqua+ioni fuori da 298 +- 1 K;
#      STOP: proteina sopra 299.5 K, o un gruppo sotto 297 K, o acqua+ioni fuori da 298 +- 1 K, o T della proteina
#            salita di piu' di 0.5 K rispetto al tratto lungo precedente (tratti da 200 ns in su: 10-250, 250-500, ...),
#            o T non leggibili.
#    Con STOP lo script scrive STOP_T_<sistema>.txt e non fa partire tratti nuovi DI QUEL SISTEMA (l'altro prosegue;
#    i tratti in corso finiscono): va segnalato al revisore e a francesco. Per riprendere, dopo la decisione,
#    rinominare STOP_T_<sistema>.txt: le GPU in attesa ripartono da sole entro un minuto.
#
# Uso (dalla cartella complesso_fab):
#   setsid nohup bash 06_repliche.sh > 06.out 2>&1 < /dev/null &            lancia (o riprende dai checkpoint)
#   bash 06_repliche.sh controlla                                            stato: ns raggiunti e avvisi LINCS
# Variabili (default = piano consigliato 31 ott):
#   GPU_COMPLESSO="0 1 2 3 4"  GPU_LIBERO="5 6 7"   (una GPU puo' comparire due volte: due repliche sulla stessa GPU;
#                                                   GPU_LIBERO="" = il libero lo fanno le GPU del complesso, con AIUTA_LIBERO=1)
#   NREP_COMPLESSO=5 NREP_LIBERO=10   NS_COMPLESSO=500 NS_LIBERO=1000   BLOCCO=250 (ns per tratto)
#   DT=4 (oppure 2)   NT=40 thread per replica
#   SCADENZA="AAAA-MM-GG hh:mm"  facoltativa (vuota = nessuna): le repliche si fermano pulite prima di quell'ora e nessun
#                                lavoro nuovo parte dopo; al rilancio (anche con meno GPU) ripartono dal checkpoint.
# Partenza scaglionata (2 nov 2026; il protocollo non cambia, cambia solo quale GPU fa quale tratto):
#   INIZIO_LIBERO="AAAA-MM-GG hh:mm"  le GPU di GPU_LIBERO non prendono lavoro prima di quell'ora (vuota = subito)
#   PRIMA_LIBERO=1     i tratti iniziali del libero (equil_s5 + primi 10 ns) vanno in testa alla coda delle GPU del complesso,
#                      cosi' la regola LINCS e' verificata su ENTRAMBI i sistemi prima di investire GPU-giorni (default 0)
#   AIUTA_LIBERO=1     una GPU del complesso, finita la coda del complesso, prende tratti del libero (default 1)
#   Ogni GPU parte solo quando non ci sono processi di calcolo di altri utenti su di essa (nvidia-smi).
#   TURNI=alternati    (ordine fissato dal revisore il 2 nov 2026, prima dei dati) una coda unica per tutte le GPU:
#                      equil_s5 + primi 10 ns di tutte le repliche (prima il libero, poi il complesso); poi a ogni traguardo
#                      (BLOCCO, 2*BLOCCO, ...) tutte le repliche dei due sistemi, prima il complesso e poi il libero.
#                      Con i default: fase A = 15 repliche fino a 500 ns, fase B = solo il libero fino a 1 us.
#                      Cosi', a qualunque ora si fermino le GPU, il libero non resta indietro rispetto al complesso.
#                      In questo modo PRIMA_LIBERO e AIUTA_LIBERO non servono; GPU_LIBERO (con INIZIO_LIBERO) sono GPU in piu'
#                      sulla stessa coda. Default TURNI=separati (code separate come sopra).
set -euo pipefail
QUI=$(cd "$(dirname "$0")" && pwd)
AZIONE=${1:-esegui}
GMX=${GMX:-gmx}; NT=${NT:-40}; DT=${DT:-4}
NREP_COMPLESSO=${NREP_COMPLESSO:-5}; NREP_LIBERO=${NREP_LIBERO:-10}
NS_COMPLESSO=${NS_COMPLESSO:-500}; NS_LIBERO=${NS_LIBERO:-1000}; BLOCCO=${BLOCCO:-250}
read -r -a GC <<< "${GPU_COMPLESSO-0 1 2 3 4}"; read -r -a GL <<< "${GPU_LIBERO-5 6 7}"   # GPU_LIBERO="" = nessuna GPU solo per il libero
INIZIO_LIBERO=${INIZIO_LIBERO:-}; PRIMA_LIBERO=${PRIMA_LIBERO:-0}; AIUTA_LIBERO=${AIUTA_LIBERO:-1}; TURNI=${TURNI:-separati}
case $TURNI in alternati|separati) ;; *) echo "TURNI deve essere alternati o separati"; exit 1 ;; esac
[ -z "$INIZIO_LIBERO" ] || date -d "$INIZIO_LIBERO" +%s > /dev/null || { echo "INIZIO_LIBERO non e' una data valida"; exit 1; }
STOP="$QUI/STOP_LINCS.txt"
T_RIF=${T_RIF:-298}; T_TOL=${T_TOL:-1}             # avviso: un gruppo fuori da T_RIF +- T_TOL K (acqua+ioni: stop)
T_PROT_MAX=${T_PROT_MAX:-299.5}; T_MIN=${T_MIN:-297}; T_SALITA=${T_SALITA:-0.5}   # stop (revisore, 2 nov 16:00)
case $DT in
  4) DTPS=0.004; HMR=3; NSTX=25000; NSTE=5000; PASSI_NS=250000 ;;
  2) DTPS=0.002; HMR=1; NSTX=50000; NSTE=10000; PASSI_NS=500000 ;;
  *) echo "DT deve essere 2 o 4"; exit 1 ;;
esac

tempo_ns() { if [ -f "$1" ]; then awk '/^ +Step +Time/{getline; t=$2} END{printf "%.1f\n", t/1000}' "$1"; else echo 0; fi; }
avvisi() {   # file [ps massimi]: avvisi LINCS fino a quel tempo
  if [ -f "$1" ]; then awk -v m="${2:-1e30}" '/LINCS WARNING/{t=0; for(i=1;i<NF;i++) if($i=="time") t=$(i+1)+0; if(t<=m) n++} END{print n+0}' "$1"; else echo 0; fi
}
minore() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a+0 < b+0)}'; }
ore() {      # ore per -maxh fino alla SCADENZA (meno 10 minuti); "stop" se non conviene partire
  [ -z "${SCADENZA:-}" ] && { echo ""; return; }
  awk -v s=$(( $(date -d "$SCADENZA" +%s) - $(date +%s) )) 'BEGIN{h=s/3600-0.17; if(h<0.25) print "stop"; else printf "%.2f\n", h}'
}
sostituisci() {   # template uscita CONT GVEL SEME NSTEPS
  sed -e "s|\${CONT}|$3|g" -e "s|\${GVEL}|$4|g" -e "s|\${GSEED}|$5|g" -e "s|\${LSEED}|$5|g" -e "s|\${DT}|$DTPS|g" \
      -e "s|\${NSTEPS}|$6|g" -e "s|\${HMR}|$HMR|g" -e "s|\${NSTX}|$NSTX|g" -e "s|\${NSTE}|$NSTE|g" "$1" > "$2"
  ! grep -q '\${' "$2" || { echo "segnaposto non sostituiti in $2"; return 1; }
}
ferma() { echo "$(date '+%F %T') $1" | tee -a "$STOP"; }
controlla_T() {   # sistema replica traguardo (nella cartella run<r>): T medie del tratto appena finito, dall'edr
  local S=$1 r=$2 fine=$3 ini=0 tp tw e prec="" lungo esito
  [ -s controllo_T.txt ] && ini=$(awk 'END{print $2+0}' controllo_T.txt)
  minore "$ini" "$fine" || return 0          # tratto gia' controllato
  e=energy_T_${ini}_${fine}ns.log
  printf "T-Protein_Other\nT-Water_and_ions\n\n" | $GMX energy -f run$r.edr -b $(( ini * 1000 )) -e $(( fine * 1000 )) -o T_${ini}_${fine}ns.xvg > $e 2>&1 || true
  read -r tp tw <<< "$(awk '$1=="T-Protein_Other"{p=$2} $1=="T-Water_and_ions"{w=$2} END{print p, w}' $e)"
  if [ -z "$tp" ] || [ -z "$tw" ]; then
    echo "$(date '+%F %T') $S run$r, tratto $ini-$fine ns: STOP, T non lette dall'edr ($S/run$r/$e)" | tee -a "$QUI/STOP_T_$S.txt"; return 0
  fi
  [ -s controllo_T.txt ] && prec=$(awk '$2 - $1 >= 200 {p=$3} END{print p}' controllo_T.txt)   # proteina, ultimo tratto lungo
  lungo=$(( fine - ini >= 200 ? 1 : 0 ))
  esito=$(awk -v p="$tp" -v w="$tw" -v t="$T_RIF" -v d="$T_TOL" -v pmax="$T_PROT_MAX" -v tmin="$T_MIN" -v sal="$T_SALITA" \
              -v q="$prec" -v lungo=$lungo 'function a(x){return x<0?-x:x}
    BEGIN{m=""
          if (p > pmax)                         m = m ";proteina>" pmax
          if (p < tmin)                         m = m ";proteina<" tmin
          if (a(w - t) > d)                     m = m ";acqua+ioni-fuori-" t "+-" d
          if (lungo && q != "" && p - q > sal)  m = m sprintf(";proteina-salita-%.2fK", p - q)
          if (m != "") print "STOP(" substr(m, 2) ")"
          else if (a(p - t) > d)                print "AVVISO(proteina-fuori-" t "+-" d ")"
          else                                  print "ok"}')
  echo "$ini $fine $tp $tw $esito" >> controllo_T.txt
  case $esito in
    STOP*)   echo "$(date '+%F %T') $S run$r, tratto $ini-$fine ns: proteina $tp K, acqua+ioni $tw K: $esito" | tee -a "$QUI/STOP_T_$S.txt" ;;
    AVVISO*) echo "$(date '+%F %T') $S run$r, tratto $ini-$fine ns: proteina $tp K, acqua+ioni $tw K: $esito (prosegue)" ;;
  esac
  return 0
}

replica() {       # sistema replica gpu traguardo_ns (fino a dove portare la produzione in questo tratto)
  local S=$1 r=$2 gpu=$3 fine=$4 seme=$(( $2 * 11111 )) ns h CPI n T
  if [ $S = complesso ]; then ns=$NS_COMPLESSO; else ns=$NS_LIBERO; fi
  mkdir -p "$QUI/$S/run$r"; cd "$QUI/$S/run$r"
  until [ -f ../equil_s4.gro ]; do [ -f "$STOP" ] && return 0; sleep 60; done
  # --- equil_s5: 20 ns ---
  if minore "$(tempo_ns equil_s5.log)" 20; then
    if [ ! -f equil_s5.tpr ]; then
      sostituisci "$QUI/mdp/equil_s5.mdp.in" equil_s5.mdp no yes $seme $(( 20 * PASSI_NS )) || return 1
      $GMX grompp -f equil_s5.mdp -c ../equil_s4.gro -r ../min_c3.gro -p ../topol.top -n ../index.ndx -o equil_s5.tpr -maxwarn 0 > grompp_equil_s5.log 2>&1 \
        || { echo "$S run$r: grompp s5 fallito ($S/run$r/grompp_equil_s5.log)"; return 1; }
    fi
    h=$(ore); [ "$h" = stop ] && return 0
    CPI=""; [ -f equil_s5.cpt ] && CPI="-cpi equil_s5.cpt"
    $GMX mdrun -deffnm equil_s5 -ntmpi 1 -ntomp $NT -gpu_id $gpu $CPI ${h:+-maxh $h} >> mdrun_equil_s5.out 2>&1 || { echo "$S run$r: mdrun s5 fallito"; return 1; }
    minore "$(tempo_ns equil_s5.log)" 20 && { echo "$S run$r: s5 fermato a $(tempo_ns equil_s5.log) ns (scadenza?), riprende al rilancio"; return 0; }
  fi
  n=$(avvisi mdrun_equil_s5.out); [ "$n" = 0 ] || { ferma "$S run$r: $n avvisi LINCS in equil_s5"; return 0; }
  # --- produzione fino al traguardo di questo tratto ---
  if [ ! -f run$r.tpr ]; then
    sostituisci "$QUI/mdp/md.mdp.in" run$r.mdp yes no $seme $(( ns * PASSI_NS )) || return 1
    $GMX grompp -f run$r.mdp -c equil_s5.gro -t equil_s5.cpt -p ../topol.top -n ../index.ndx -o run$r.tpr -maxwarn 0 > grompp_run$r.log 2>&1 \
      || { echo "$S run$r: grompp produzione fallito"; return 1; }
  fi
  if minore "$(tempo_ns run$r.log)" $fine; then
    T=run${r}_fino_${fine}ns.tpr
    [ -f $T ] || $GMX convert-tpr -s run$r.tpr -until $(( fine * 1000 )) -o $T > convert_${fine}ns.log 2>&1 || { echo "$S run$r: convert-tpr fallito"; return 1; }
    h=$(ore); [ "$h" = stop ] && return 0
    CPI=""; [ -f run$r.cpt ] && CPI="-cpi run$r.cpt -append"
    $GMX mdrun -s $T -deffnm run$r -ntmpi 1 -ntomp $NT -gpu_id $gpu -cpt 15 $CPI ${h:+-maxh $h} >> mdrun_run$r.out 2>&1 \
      || { echo "$S run$r: mdrun produzione fallito"; return 1; }
    minore "$(tempo_ns run$r.log)" $fine && { echo "$S run$r: fermata a $(tempo_ns run$r.log) ns (scadenza?), riprende al rilancio"; return 0; }
  fi
  n=$(avvisi mdrun_run$r.out 10000); [ "$n" = 0 ] || { ferma "$S run$r: $n avvisi LINCS nei primi 10 ns"; return 0; }
  controlla_T $S $r $fine
  echo "$(date '+%F %T') $S run$r: a $(tempo_ns run$r.log) ns (tratto fino a $fine, GPU $gpu, seme $seme)"
}

prendi() {        # primo tratto in coda la cui replica e' libera; "attendi" se ce ne sono ma tutte occupate
  local q=$1 j="" riga k=0
  until mkdir "$q.lock" 2>/dev/null; do sleep 1; done
  while read -r riga; do
    k=$((k + 1)); set -- $riga
    [ -f "$QUI/STOP_T_$1.txt" ] && continue          # sistema fermato dal controllo della temperatura
    if mkdir "$QUI/$1/.run$2.occupata" 2>/dev/null; then j=$riga; sed -i "${k}d" "$q"; break; fi
  done < "$q"
  [ -z "$j" ] && [ -s "$q" ] && j=attendi
  rmdir "$q.lock"; echo "$j"
}
altri_su_gpu() {  # gpu: numero di processi di calcolo di altri utenti sulla GPU (0 se nvidia-smi manca)
  local p n=0
  command -v nvidia-smi > /dev/null || { echo 0; return; }
  for p in $(nvidia-smi -i "$1" --query-compute-apps=pid --format=csv,noheader 2>/dev/null | tr -d ' '); do
    [ "$(ps -o uid= -p "$p" 2>/dev/null | tr -d ' ')" = "$(id -u)" ] || n=$((n + 1))
  done
  echo $n
}
attendi_gpu() {   # gpu inizio: aspetta l'ora di inizio (se data) e che la GPU sia libera da processi di altri; 1 se STOP o scadenza
  local gpu=$1 inizio=$2 detto=""
  while :; do
    { [ -f "$STOP" ] || [ "$(ore)" = stop ]; } && return 1
    if [ -z "$inizio" ] || [ "$(date +%s)" -ge "$(date -d "$inizio" +%s)" ]; then
      [ "$(altri_su_gpu "$gpu")" = 0 ] && { [ -z "$detto" ] || echo "$(date '+%F %T') GPU $gpu: libera, parto"; return 0; }
    fi
    [ -z "$detto" ] && { echo "$(date '+%F %T') GPU $gpu: attendo ${inizio:+le $inizio e }che non ci siano processi di altri utenti"; detto=1; }
    sleep 60
  done
}
lavora() {        # gpu inizio coda [coda2]: la seconda coda solo quando la prima e' vuota (se la prima ha tratti occupati, attende)
  local gpu=$1 inizio=$2 j q; shift 2
  local code=("$@")      # copia: piu' sotto "set -- $j" sovrascrive gli argomenti
  attendi_gpu "$gpu" "$inizio" || return 0
  while [ ! -f "$STOP" ] && [ "$(ore)" != stop ]; do
    j=""
    for q in "${code[@]}"; do j=$(prendi "$q"); [ -n "$j" ] && break; done
    [ -z "$j" ] && break
    [ "$j" = attendi ] && { sleep 60; continue; }
    set -- $j       # sistema replica traguardo
    ( replica $1 $2 $gpu $3 ) || echo "$(date '+%F %T') GPU $gpu: tratto $j interrotto da un errore, la coda prosegue"
    rmdir "$QUI/$1/.run$2.occupata" 2>/dev/null || true
  done
}

if [ "$AZIONE" = controlla ]; then
  printf "%-10s %-6s %10s %8s %12s %10s\n" sistema replica s5_ns avvisi_s5 prod_ns avvisi_prod
  for S in complesso libero; do
    if [ $S = complesso ]; then N=$NREP_COMPLESSO; else N=$NREP_LIBERO; fi
    for r in $(seq 1 $N); do d="$QUI/$S/run$r"
      printf "%-10s %-6s %10s %8s %12s %10s\n" $S run$r "$(tempo_ns $d/equil_s5.log)" "$(avvisi $d/mdrun_equil_s5.out)" "$(tempo_ns $d/run$r.log)" "$(avvisi $d/mdrun_run$r.out)"
    done
  done
  [ -f "$STOP" ] && { echo; echo "STOP_LINCS.txt presente:"; cat "$STOP"; }
  echo; echo "T media per tratto (K), dall'edr. Avviso: fuori da $T_RIF +- $T_TOL. Stop: proteina > $T_PROT_MAX, gruppo < $T_MIN, acqua+ioni fuori, salita proteina > $T_SALITA:"
  printf "%-10s %-6s %12s %10s %12s  %s\n" sistema replica tratto_ns proteina acqua+ioni esito
  for S in complesso libero; do
    if [ $S = complesso ]; then N=$NREP_COMPLESSO; else N=$NREP_LIBERO; fi
    for r in $(seq 1 $N); do f="$QUI/$S/run$r/controllo_T.txt"
      [ -f "$f" ] && awk -v s=$S -v r=run$r '{printf "%-10s %-6s %12s %10s %12s  %s\n", s, r, $1 "-" $2, $3, $4, $5}' "$f"
    done
  done
  for f in "$QUI"/STOP_T_*.txt; do [ -f "$f" ] && { echo; echo "$(basename "$f") presente:"; cat "$f"; }; done
  exit 0
fi

if pgrep -u "$(id -u)" -f "^([^ ]*/)?gmx(_mpi)? mdrun .*(-deffnm equil_s5|-s run[0-9]+(_10ns)?\.tpr)" > /dev/null; then
  echo "ci sono gia' repliche in corso (mdrun di equil_s5 o run*): non rilancio per non duplicarle"; exit 1
fi
[ -f "$STOP" ] && { echo "c'e' STOP_LINCS.txt: decidere (regola del revisore) e cancellarlo prima di rilanciare"; cat "$STOP"; exit 1; }
for f in "$QUI"/STOP_T_*.txt; do [ -f "$f" ] && { echo "ATTENZIONE: $(basename "$f") presente, quel sistema non prende tratti nuovi:"; cat "$f"; }; done
echo "$(date '+%F %T') DT=$DT fs, complesso ${NREP_COMPLESSO}x${NS_COMPLESSO} ns su GPU ${GC[*]}, libero ${NREP_LIBERO}x${NS_LIBERO} ns su GPU ${GL[*]}" \
     "(inizio ${INIZIO_LIBERO:-subito}), TURNI=$TURNI PRIMA_LIBERO=$PRIMA_LIBERO AIUTA_LIBERO=$AIUTA_LIBERO NT=$NT, scadenza ${SCADENZA:-nessuna}" | tee -a "$QUI/registro_06.txt"
QC="$QUI/coda_complesso.txt"; QL="$QUI/coda_libero.txt"; rm -rf "$QC.lock" "$QL.lock" "$QUI"/*/.run*.occupata
turni() {         # sistema nrep ns: un tratto da 10 ns (s5 + primi 10 ns sorvegliati) per tutte, poi BLOCCO alla volta
  local S=$1 N=$2 ns=$3 t r
  for t in 10 $(seq $BLOCCO $BLOCCO $ns); do for r in $(seq 1 $N); do echo "$S $r $t"; done; done
  [ $(( ns % BLOCCO )) = 0 ] || for r in $(seq 1 $N); do echo "$S $r $ns"; done
}
alternati() {     # coda unica: ordinata per traguardo; a 10 ns prima il libero, agli altri traguardi prima il complesso
  { turni complesso $NREP_COMPLESSO $NS_COMPLESSO; turni libero $NREP_LIBERO $NS_LIBERO; } |
    awk '{o = ($1 == "complesso") ? 0 : 1; if ($3 == 10) o = 1 - o; print $3, o, $2, $1}' | sort -k1,1n -k2,2n -k3,3n | awk '{print $4, $3, $1}'
}
if [ "$TURNI" = alternati ]; then
  alternati > "$QC"; : > "$QL"; CODE_C=("$QC"); CODE_L=("$QC")
else
  turni complesso $NREP_COMPLESSO $NS_COMPLESSO > "$QC"; turni libero $NREP_LIBERO $NS_LIBERO > "$QL"
  if [ "$PRIMA_LIBERO" = 1 ]; then
    { grep -E '^libero [0-9]+ 10$' "$QL" || true; cat "$QC"; } > "$QC.nuova"; mv "$QC.nuova" "$QC"; sed -i -E '/^libero [0-9]+ 10$/d' "$QL"
  fi
  CODE_C=("$QC"); [ "$AIUTA_LIBERO" = 1 ] && CODE_C+=("$QL"); CODE_L=("$QL")
fi
PIDS=()
for g in ${GC[@]+"${GC[@]}"}; do lavora $g "" "${CODE_C[@]}" & PIDS+=($!); sleep 2; done
for g in ${GL[@]+"${GL[@]}"}; do lavora $g "$INIZIO_LIBERO" "${CODE_L[@]}" & PIDS+=($!); sleep 2; done
for p in "${PIDS[@]}"; do wait $p || true; done
echo "$(date '+%F %T') fine" | tee -a "$QUI/registro_06.txt"
bash "$0" controlla
