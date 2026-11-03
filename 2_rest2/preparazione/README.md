# Preparazione – REST2

Ordine d'uso:

| # | Script | Scopo |
|---|---|---|
| 0 | `00_compila_gromacs_plumed.sh` | Compila GROMACS 2024.6 con MPI, CUDA 12.8 e patch PLUMED 2.10.2 in un ambiente conda, senza root (GPU Blackwell sm_120). Scrive anche `attiva_gmx2024.sh` per caricare l'ambiente. |
| 1 | `01_muta_nmr.py` | Porta il medoide NMR (9HW9 modello 6) alla sequenza del costrutto X-ray (8 mutazioni: troncamento al CB o innesto della catena laterale X-ray), backbone NMR invariato. |
| 2 | `02_controlla.py` | Controlli sulla struttura NMR mutata prima di pdb2gmx: sequenza, atomi mancanti/spostati, contatti stretti, S–S, χ1, carica, minimizzazione di prova (OpenMM, solvente implicito). |
| – | `eq_rest2_nmr.sh` | Sistema NMR mutato: minimizzazioni min_c1–c3 ed equil_s1–s4 (stessi mdp delle MD classiche, `1_md_classiche_apo/`). Il sistema X-ray parte dall'equil_s4 delle MD classiche. pdb2gmx e solvatazione come in `eq_nmr_ox.sh`. |
| – | `pp.mdp` | mdp per `gmx grompp -pp` → `processed.top` (topologia espansa di ogni sistema). |
| 3 | `03_marca_caldi.py` | Marca con "_" i tipi atomici dei residui caldi in `processed.top` (per `plumed partial_tempering`). |
| 4 | `04_prepara_rest2.sh` | Crea una cartella per replica con topologia scalata (λ = 310.15/T), mdp e tpr; con Tmax = 310.15 fa il test a λ = 1. |
| – | `rest2.mdp.in` | Template mdp di REST2 (stessi parametri della produzione classica, solo proteina ogni 10 ps). |
| 7 | `07_kappa_rest3.py` | Correzione REST3: vdW fra atomi caldi e solvente × κ (`[ nonbond_params ]`). |
| 8 | `08_prepara_rest3.sh` | Come 04 ma per REST3 (scala di κ fino a κ_max = 1.05). Preparato, non usato in produzione. |
