"""Correzione REST3 su una topologia REST2: vdW fra atomi caldi e solvente moltiplicato per kappa.

Uso:  plumed partial_tempering LAMBDA < processed_caldi.top | python3 07_kappa_rest3.py - LAMBDA KAPPA > topol.top
      python3 07_kappa_rest3.py topol_rest2.top LAMBDA KAPPA > topol.top

REST3 (Zhang, Liu, Chen, JCTC 2023, 19, 1602) tiene REST2 com'e' (caldo-caldo x lambda, caldo-freddo x sqrt(lambda),
cariche calde x sqrt(lambda)) e moltiplica per kappa solo il vdW fra soluto caldo e solvente. Nel loro lavoro
tutta la proteina e' calda, quindi "freddo" e "solvente" coincidono; qui sono caldi solo i loop, e kappa va
applicato solo fra atomi caldi e acqua/ioni: il vdW fra loop caldi e core freddo resta x sqrt(lambda) come in REST2.
(Lo script del loro tutorial scala tutte le coppie caldo-freddo: con i soli loop caldi rafforzerebbe anche il
contatto loop-core.)

Si aggiunge una sezione [ nonbond_params ] prima della prima [ moleculetype ], una riga per ogni coppia
(tipo caldo usato, tipo del solvente): sigma = media aritmetica, epsilon = kappa*sqrt(eps_caldo*eps_solvente),
con i sigma ed epsilon dei tipi caldi scritti da partial_tempering (eps gia' x lambda). Con kappa = 1 le righe
coincidono con la regola di combinazione e la topologia equivale a quella REST2.
Solvente = molecole di [ molecules ] senza atomi caldi (qui SOL, NA, CL). Lo script si ferma se la regola di
combinazione non e' la 2 (AMBER), se c'e' gia' una [ nonbond_params ], se un tipo del solvente e' usato anche
dagli atomi freddi della proteina, o se gli epsilon caldi non sono lambda volte quelli freddi (lambda sbagliato).
"""
import sys, math, re

SUFF = "_"

args = sys.argv[1:]
if len(args) != 3:
    raise SystemExit(__doc__.split("\n\n")[1])
IN, LAM, KAPPA = args[0], float(args[1]), float(args[2])
righe = (sys.stdin if IN == "-" else open(IN)).readlines()

# lettura: regola di combinazione, tipi (sigma, eps), tipi usati da ogni molecola, molecole presenti
comb, tipi, tipi_mol, presenti = None, {}, {}, []
sezione, mol, i_mol = None, None, None
for k, l in enumerate(righe):
    s = l.split(";")[0].split()
    if not s:
        continue
    m = re.match(r"\[\s*(\w+)\s*\]", " ".join(s))
    if m:
        sezione = m.group(1)
        if sezione == "nonbond_params":
            raise SystemExit("c'e' gia' una sezione [ nonbond_params ]: questo script non la unisce")
        if sezione == "moleculetype" and i_mol is None:
            i_mol = k
        continue
    if sezione == "defaults":
        comb = int(s[1])
    elif sezione == "atomtypes":
        tipi[s[0]] = (float(s[-2]), float(s[-1]))
    elif sezione == "moleculetype":
        mol = s[0]; tipi_mol[mol] = set()
    elif sezione == "atoms" and len(s) >= 7:
        tipi_mol[mol].add(s[1])
    elif sezione == "molecules" and int(s[1]) > 0:
        presenti.append(s[0])

if comb != 2:
    raise SystemExit(f"regola di combinazione {comb}: lo script e' scritto per la 2 (sigma aritmetico, eps geometrico)")
calde = [m for m in presenti if any(t.endswith(SUFF) for t in tipi_mol[m])]
solvente = [m for m in presenti if m not in calde]
if not calde:
    raise SystemExit("nessun atomo caldo (tipo che finisce con '_'): la topologia non e' marcata")
caldi = sorted({t for m in calde for t in tipi_mol[m] if t.endswith(SUFF)})
freddi_prot = {t for m in calde for t in tipi_mol[m] if not t.endswith(SUFF)}
solv = sorted({t for m in solvente for t in tipi_mol[m]})
comuni = freddi_prot & set(solv)
if comuni:
    raise SystemExit(f"tipi usati sia dal solvente sia da atomi freddi della proteina: {sorted(comuni)}")

# controllo di lambda: eps del tipo caldo = lambda * eps del tipo freddo (partial_tempering scrive 6 cifre)
for t in caldi:
    sc, ec = tipi[t]; sf, ef = tipi[t[:-len(SUFF)]]
    if abs(sc - sf) > 1e-5 * max(sf, 1e-3) or abs(ec - LAM * ef) > 1e-5 * max(ef, 1e-3):
        raise SystemExit(f"tipo {t}: sigma {sc}, eps {ec}; atteso sigma {sf}, eps {LAM * ef:.6g} con lambda {LAM}")

nb = ["[ nonbond_params ]\n",
      f"; REST3: vdW fra atomi caldi e solvente x kappa = {KAPPA:g} (lambda = {LAM:g}); 07_kappa_rest3.py\n",
      "; i  j  func  sigma (nm)  epsilon (kJ/mol)\n"]
for t in caldi:
    for w in solv:
        (si, ei), (sj, ej) = tipi[t], tipi[w]
        nb.append(f"{t:<12} {w:<12} 1  {0.5 * (si + sj):.10g}  {KAPPA * math.sqrt(ei * ej):.10g}\n")
nb.append("\n")
sys.stdout.writelines(righe[:i_mol] + nb + righe[i_mol:])
print(f"REST3 lambda {LAM:g} kappa {KAPPA:g}: {len(caldi)} tipi caldi x {len(solv)} tipi del solvente "
      f"({', '.join(solv)}; molecole {', '.join(solvente)}) = {len(caldi) * len(solv)} coppie", file=sys.stderr)
