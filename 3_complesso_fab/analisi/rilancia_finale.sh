#!/bin/bash
# Analisi finale con UNA sola versione degli script (condizione del revisore, 11 nov 2026: prima di ogni altra
# esecuzione la cartella analisi della workstation va allineata).
#
# Uso, dalla cartella complesso_fab della workstation, dopo aver copiato in analisi/ TUTTI i file di questa cartella:
#     bash analisi/rilancia_finale.sh
#
# 1. controlla le sha256 degli script contro analisi/SHA256SUMS e si ferma se anche uno solo e' diverso;
# 2. mette da parte analisi/risultati (non cancella niente): le uscite precedenti all'allineamento non valgono;
# 3. rilancia 10-18 nell'ordine (10 rifa' da solo le serie che non vengono dagli script attuali; 09 non serve,
#    le traiettorie estratte restano quelle) e si ferma al primo errore.
# 16 chiede GROMACS: se gmx non e' caricato viene saltato e lo si legge nell'ultima riga.
set -euo pipefail
cd "$(dirname "$0")/.."

if ! ( cd analisi && sha256sum -c --quiet SHA256SUMS ); then
    echo "STOP: gli script in analisi/ non sono la versione attesa (vedi le righe FAILED qui sopra)."
    echo "Ricopiare in analisi/ TUTTI i file della cartella analisi_script_attuali e rilanciare."
    exit 1
fi
echo "script: versione attesa (sha256 uguali a analisi/SHA256SUMS)"

if [ -d analisi/risultati ]; then
    vecchi="analisi/risultati_prima_dell_allineamento_$(date +%Y%m%d_%H%M%S)"
    mv analisi/risultati "$vecchi"
    echo "uscite precedenti spostate in $vecchi (non si leggono e non si citano)"
fi

GMX_BIN="${GMX:-gmx}"
saltato16=""
for passo in 10_serie 11_controlli 12_scopo1_contatti 13_scopo2_epitopo 14_scopo3_rmsf 15_figure 16_verifica_immagini \
             17_cerniera 18_figura_articolo; do
    if [ "$passo" = 16_verifica_immagini ] && ! command -v "${GMX_BIN%% *}" >/dev/null 2>&1; then
        saltato16="16_verifica_immagini.py SALTATO: ${GMX_BIN%% *} non trovato (caricare GROMACS e lanciarlo a parte)"
        echo "$saltato16"; continue
    fi
    echo; echo "=================== $passo.py"
    python3 "analisi/$passo.py"
done

echo; echo "FATTO. Ricopiare sul PC le cartelle analisi/cache e analisi/risultati."
[ -n "$saltato16" ] && echo "$saltato16"
exit 0
