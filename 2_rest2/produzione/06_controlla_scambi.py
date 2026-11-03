"""Statistiche di scambio di un run REST2, leggendo il log di mdrun (anche a run in corso).

Uso:  python3 rest2/06_controlla_scambi.py <cartella_run>
Stampa per ogni coppia di repliche vicine i tentativi, gli scambi accettati e la probabilita' media,
poi quante volte ogni configurazione ("walker") ha attraversato tutta la scala: un viaggio completo
e' 0 -> N-1 -> 0. Obiettivo: 20-30% di accettazione su tutte le coppie e ogni walker che fa viaggi
completi; una coppia molto sotto le altre e' un collo di bottiglia e la scala va ritoccata.
r0/rest2.xtc e' sempre l'ensemble non scalato (lambda = 1, 310.15 K): le cartelle tengono
l'Hamiltoniana, sono le coordinate a spostarsi.
"""
import os, sys

RUN = sys.argv[1]
scala = [l.split() for l in open(os.path.join(RUN, "scala.txt"))]
N = len(scala)
tentativi, accettati, somma_p = [0] * N, [0] * N, [0.0] * N
walker = list(range(N))                 # walker[posizione] = configurazione che sta li'
ultimo = {w: None for w in range(N)}    # ultimo estremo toccato da ogni walker (0 o N-1)
mezzi = [0] * N                         # traversate complete della scala, in un verso
n_scambi_log, passo = 0, None

for l in open(os.path.join(RUN, "r0", "rest2.log")):
    if l.startswith("Replica exchange at step"):
        passo = int(l.split()[4]); n_scambi_log += 1
    elif l.startswith("Repl ex"):
        t = l[7:].split()
        for k, tok in enumerate(t):
            if tok == "x":
                a = int(t[k - 1]); b = int(t[k + 1])
                accettati[b] += 1
                walker[a], walker[b] = walker[b], walker[a]
        for p, w in enumerate(walker):
            if p in (0, N - 1):
                if ultimo[w] is not None and ultimo[w] != p:
                    mezzi[w] += 1
                ultimo[w] = p
    elif l.startswith("Repl pr"):
        for i in range(1, N):
            campo = l[8 + 5 * (i - 1): 8 + 5 * i].strip()
            if campo:
                tentativi[i] += 1; somma_p[i] += float(campo)

if not n_scambi_log:
    raise SystemExit("nessun tentativo di scambio nel log: il run e' partito con -replex e -hrex?")
print(f"{n_scambi_log} tentativi di scambio, ultimo al passo {passo} ({passo * 0.002 / 1000:.3f} ns)")
print("coppia   T_eff (K)          tentativi  accettati  accettazione  prob. media")
for i in range(1, N):
    acc = accettati[i] / tentativi[i] if tentativi[i] else float("nan")
    pm = somma_p[i] / tentativi[i] if tentativi[i] else float("nan")
    print(f"{i - 1:2d}-{i:<2d}   {float(scala[i - 1][1]):6.1f}-{float(scala[i][1]):6.1f}   {tentativi[i]:9d}  {accettati[i]:9d}  {acc:12.2f}  {pm:11.2f}")
print("viaggi completi 0 -> N-1 -> 0 per walker: " + ", ".join(f"{w}: {mezzi[w] // 2}" for w in range(N)))
print(f"totale viaggi completi: {sum(m // 2 for m in mezzi)}")
