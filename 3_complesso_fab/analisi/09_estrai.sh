#!/usr/bin/env bash
# Estrae la sola proteina (gruppo Protein_Other di index.ndx) dalle traiettorie di produzione, dall'inizio della
# produzione fino a FINO ns, per tutte le repliche che li hanno gia' raggiunti. Non modifica i file di produzione:
# legge run<r>.xtc anche mentre mdrun lo allunga, perche' si ferma a FINO ns (-e).
#
# Uso (dalla cartella complesso_fab, con GROMACS caricato):
#   FINO=250 bash analisi/09_estrai.sh
# Variabili: GMX (default gmx), FINO (ns di produzione, default 250), SISTEMI (default "complesso libero").
# Il tempo di produzione si conta dal tinit del tpr (il primo fotogramma), quindi s5 e' sempre escluso.
# Scrive in analisi/traiettorie/:
#   proteina_<sistema>.gro               solo proteina, dalla equil_s4.gro del sistema (topologia per le analisi)
#   <sistema>_run<r>_0-<FINO>ns.xtc      solo proteina, un fotogramma ogni 100 ps, tempo come nell'xtc originale
#   provenienza_09.txt                   versione di GROMACS, tinit, comandi, sha256 di ingressi e uscite
# Una replica gia' estratta con lo stesso FINO non si rifa'. Le repliche sotto FINO ns si saltano (con un avviso).
set -euo pipefail
QUI=$(cd "$(dirname "$0")/.." && pwd)
GMX=${GMX:-gmx}; FINO=${FINO:-250}; SISTEMI=${SISTEMI:-"complesso libero"}
OUT="$QUI/analisi/traiettorie"; mkdir -p "$OUT"
PROV="$OUT/provenienza_09.txt"
command -v $GMX > /dev/null || { echo "$GMX non trovato: caricare GROMACS"; exit 1; }
[[ "$FINO" =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "FINO deve essere un numero di ns (es. 250)"; exit 1; }

tempo_ps() { awk '/^ +Step +Time/{getline; t=$2} END{printf "%.3f\n", t+0}' "$1"; }   # ultimo tempo nel log (ps)
{
  echo "=== $(date '+%F %T') 09_estrai.sh sha256 $(sha256sum "$0" | cut -c1-16), FINO=$FINO ns"
  $GMX --version 2>/dev/null | grep -E "GROMACS version|Precision" || true
} >> "$PROV"

for S in $SISTEMI; do
  D="$QUI/$S"
  [ -d "$D" ] || { echo "$S: cartella assente, salto"; continue; }
  for f in "$D/index.ndx" "$D/equil_s4.gro"; do [ -f "$f" ] || { echo "manca $f"; exit 1; }; done
  grep -q "^\[ Protein_Other \]" "$D/index.ndx" || { echo "$S: gruppo Protein_Other assente in index.ndx"; exit 1; }
  for RD in $(ls -d "$D"/run* 2>/dev/null | sort -V); do
    r=$(basename "$RD")
    [ -f "$RD/$r.tpr" ] && [ -f "$RD/$r.xtc" ] && [ -f "$RD/$r.log" ] || { echo "$S $r: tpr, xtc o log assenti, salto"; continue; }
    X="$OUT/${S}_${r}_0-${FINO}ns.xtc"
    if [ -f "$X" ]; then echo "$S $r: gia' estratta ($(basename "$X"))"; continue; fi
    # (awk legge tutto l'output: con pipefail un'uscita anticipata darebbe SIGPIPE a gmx dump)
    TINIT=$($GMX dump -s "$RD/$r.tpr" 2>/dev/null | awk '$1=="tinit"{v=$3; f=1} END{if (f) printf "%.3f\n", v+0}')
    [ -n "$TINIT" ] || { echo "$S $r: tinit non letto dal tpr"; exit 1; }
    FINE_PS=$(awk -v a="$TINIT" -v b="$FINO" 'BEGIN{printf "%.3f\n", a + 1000*b}')
    ORA=$(tempo_ps "$RD/$r.log")
    if awk -v o="$ORA" -v f="$FINE_PS" 'BEGIN{exit !(o + 0.5 < f)}'; then
      echo "$S $r: a $(awk -v o="$ORA" -v t="$TINIT" 'BEGIN{printf "%.1f", (o-t)/1000}') ns di produzione (< $FINO), salto"; continue
    fi
    # topologia della sola proteina (una volta per sistema, stessa per tutte le repliche)
    G="$OUT/proteina_$S.gro"
    if [ ! -f "$G" ]; then
      echo Protein_Other | $GMX trjconv -f "$D/equil_s4.gro" -s "$RD/$r.tpr" -n "$D/index.ndx" -o "$G" > "$OUT/trjconv_gro_$S.log" 2>&1 \
        || { echo "$S: trjconv del gro fallito ($OUT/trjconv_gro_$S.log)"; exit 1; }
      echo "$S: $(basename "$G") <- equil_s4.gro (sha256 $(sha256sum "$D/equil_s4.gro" | cut -c1-16)), index.ndx (sha256 $(sha256sum "$D/index.ndx" | cut -c1-16)): sha256 $(sha256sum "$G" | cut -c1-16), atomi $(sed -n 2p "$G" | tr -d ' ')" | tee -a "$PROV"
    fi
    T="$OUT/.${S}_${r}.tmp.xtc"; rm -f "$T"
    echo "$S $r: estraggo 0-$FINO ns (tinit $TINIT ps, fino a $FINE_PS ps) ..."
    echo Protein_Other | $GMX trjconv -f "$RD/$r.xtc" -s "$RD/$r.tpr" -n "$D/index.ndx" -e "$FINE_PS" -o "$T" > "$OUT/trjconv_${S}_${r}.log" 2>&1 \
      || { echo "$S $r: trjconv fallito ($OUT/trjconv_${S}_${r}.log)"; rm -f "$T"; exit 1; }
    # controllo con gmx check: primo fotogramma, numero di fotogrammi e passo (tabella "Time"); ultimo fotogramma da
    # "Last frame" di gmx check (GROMACS 2023), altrimenti da "Last written" di trjconv (2026, dove gmx check non lo
    # scrive piu'). Deve valere: primo = tinit, ultimo = tinit + FINO, fotogrammi = (ultimo - primo) / passo + 1.
    CK="$OUT/check_${S}_${r}.log"; TL="$OUT/trjconv_${S}_${r}.log"
    $GMX check -f "$T" > "$CK" 2>&1 || true
    PRIMO=$( { tr '\r' '\n' < "$CK" | grep -o 'Reading frame *0 time *[0-9.eE+-]*' || true; } | awk 'NR==1{print $5}')
    read -r NF PASSO <<< "$(awk '$1=="Time" && NF==3 {print $2, $3}' "$CK")"
    ULT=$( { tr '\r' '\n' < "$CK" | grep -o 'Last frame *[0-9]* *time *[0-9.eE+-]*' || true; } | awk 'END{print $5}')
    [ -n "$ULT" ] || ULT=$( { tr '\r' '\n' < "$TL" | grep -o 'Last written: *frame *[0-9]* *time *[0-9.eE+-]*' || true; } | awk 'END{print $6}')
    if ! awk -v p="$PRIMO" -v u="$ULT" -v n="$NF" -v d="$PASSO" -v a="$TINIT" -v f="$FINE_PS" 'BEGIN{
           if (p == "" || u == "" || n == "" || d + 0 <= 0) exit 1
           if (p - a > 0.5 || a - p > 0.5 || u - f > 0.5 || f - u > 0.5) exit 1
           k = (u - p) / d + 1; exit !(k - n < 0.01 && n - k < 0.01) }'; then
      echo "$S $r: controllo dell'estrazione fallito: primo fotogramma ${PRIMO:-?} ps (atteso $TINIT), ultimo ${ULT:-?} ps" \
           "(atteso $FINE_PS), ${NF:-?} fotogrammi ogni ${PASSO:-?} ps (vedi $CK e $TL)"; rm -f "$T"; exit 1
    fi
    mv "$T" "$X"
    echo "$S $r: tinit $TINIT ps, $NF fotogrammi fino a $ULT ps; tpr sha256 $(sha256sum "$RD/$r.tpr" | cut -c1-16); $(basename "$X") sha256 $(sha256sum "$X" | cut -c1-16)" | tee -a "$PROV"
  done
done
echo "fatto: file in analisi/traiettorie/ (provenienza in provenienza_09.txt)"
