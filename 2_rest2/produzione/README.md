# Produzione – REST2

```
source attiva_gmx2024.sh
setsid nohup bash 05_lancia_rest2.sh <cartella_run> ... &     # lancio o ripresa da checkpoint
python3 06_controlla_scambi.py <cartella_run>                 # accettazione e viaggi, anche in corsa
bash 09_allunga_rest2.sh <cartella_run> <ns_totali>           # a run fermo, poi di nuovo 05
```

| Script | Scopo |
|---|---|
| `05_lancia_rest2.sh` | Lancia tutte le repliche in un solo `mpirun` (`-multidir -replex 500 -hrex -notunepme`), una o più repliche per GPU; `setsid` perché Open MPI ignora `nohup`. Riprende da `rest2.cpt`. |
| `06_controlla_scambi.py` | Dal log di r0: tentativi e accettazione per coppia di repliche vicine, viaggi completi 0→N−1→0 di ogni walker. Obiettivo 20–30% su tutte le coppie. |
| `09_allunga_rest2.sh` | Allunga un run finito a una nuova durata totale (nuovi tpr per tutte le repliche, il vecchio resta come `rest2_fino_<ns>ns.tpr`). |
| `COMANDI.md` | Comandi completi usati sulla workstation: compilazione, topologie, test λ = 1, taratura della scala, produzione, NMR, REST3, allungamento a 400 ns. |
