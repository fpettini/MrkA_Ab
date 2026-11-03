# Analisi – REST2

Si analizzano r0 (ensemble non scalato, 310.15 K) e r15 (replica più calda, solo controllo).

```
01_estrai.py ─► cache/<sistema>_r0.npz, _r15.npz
02_controlli.py ─► controlli.txt/.json        (r15 vs r0: core, β, Rg loop, S–S)
03_popolazioni.py ─► popolazioni.txt/.json, riferimenti.txt, distanze_<sistema>.npz
```

| Script | Scopo |
|---|---|
| `01_estrai.py` | Ricompone la proteina attraverso le PBC e salva backbone + SG 43/87 di r0 e r15 per X-ray e NMR. Terzo argomento opzionale: nome del run (per REST3). |
| `02_controlli.py` | Controlli preliminari: RMSD del core, frazione β (DSSP), Rg dei loop caldi, distanza SG–SG; scarto iniziale 20 ns (`--t0`). |
| `03_popolazioni.py` | Analisi principale: distanze di ogni loop dai riferimenti X-ray e NMR (modello 6; controllo = modello più vicino dei 20), popolazioni pX/pN/altro, test di equivalenza (soglia 0.15) fra i run e fra le metà, sensibilità a scarto 120 ns e a soglie 0.75/1.25·D. |
| `provenienza.py` | Scrive la provenienza (sha256, versioni, argomenti) in testa a ogni txt/json e nelle cache. |
| `SPEC.md` | Specifica del piano di analisi approvato (criteri fissati prima dei dati). |
| `LEGGIMI.md` | Tabella passo → script → input/output e comandi per rifare l'analisi. |
