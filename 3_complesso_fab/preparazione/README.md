# Preparazione – complesso MrkA–Fab

| # | Script / file | Scopo |
|---|---|---|
| 1 | `01_ricostruisci_fab.py` | Completa il Fab (PDBFixer: Gln1, Ala9–Gly10, Lys218–Gly223 della pesante, Cys214 della leggera), minimizza solo gli atomi nuovi (OpenMM, GB) e forma il disolfuro L214–H220; scrive un controllo (scarti, SG–SG, ω, chiralità). |
| 1b | `01b_acque_interne.py` | Seleziona le acque cristallografiche sepolte in MrkA (SASA O < 5%, criteri geometrici del revisore) da mettere identiche in complesso e libero. |
| 2 | `02_prepara_sistemi.sh` | Per complesso, libero (e fv, solo test): pdb2gmx (His con tautomeri fissati), acque cristallografiche, dodecaedro `-d 1.3`, solvatazione, KCl 0.2 M, min_c1–c3. FF/acqua via variabili (`FF=amber19sb WATER=opc SOLV=tip4p.gro`). |
| 2a | `02a_aiuto.py` | Funzioni di supporto a 02: divide le catene, aggiunge le acque cristallografiche, controlla che genion non le sostituisca, crea i gruppi index `Protein_Other` e `SOL_bulk`. |
| – | `mdp/min_c1-3.mdp`, `mdp/ions.mdp` | Minimizzazioni (steep, restraint 1000) e mdp per genion. |

Controlli del campo di forze (richiesti dal revisore prima della produzione):

| Script | Scopo |
|---|---|
| `controlli_ff/04_controllo_ff.sh` | Energia a singolo punto della sola proteina nel vuoto con GROMACS (`mdp/vuoto.mdp`). |
| `controlli_ff/04b_confronta_openmm.py` | Confronta i termini legati e 1-4 con OpenMM (amber19-all.xml); criterio: scarto relativo < 1e-4. |
| `controlli_ff/04c_confronta_parametri.py` | Confronta cariche, σ, ε della proteina e i parametri di OPC, K⁺, Cl⁻ fra GROMACS e OpenMM. |
| `densita/box_opc.py` | Densità di OPC pura e OPC + KCl con OpenMM (stesse impostazioni dell'mdp). |
| `densita/densita_s4.py` | Densità globale e dell'acqua di bulk del sistema GROMACS dopo l'equilibrazione NPT, da confrontare con box_opc. |
| `mdp/vuoto.mdp` | mdp per l'energia nel vuoto usata da 04. |
