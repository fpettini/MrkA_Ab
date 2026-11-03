"""Marca gli atomi "caldi" di REST2 nella topologia espansa (grompp -pp).

Uso:  python3 03_marca_caldi.py processed.top processed_caldi.top [--caldi 59-77,87-101,...]

partial_tempering di PLUMED scala solo gli atomi il cui tipo, nella sezione [ atoms ], finisce con "_".
Qui si aggiunge "_" agli atomi dei residui caldi della proteina (prima [ moleculetype ] del file);
acqua e ioni restano freddi. Di default i residui caldi sono i loop che divergono fra le repliche,
escluso il tratto 44-58 che si affaccia sul solco del filamento donatore.
"""
import sys

CALDI = "59-77,87-101,112-118,124-133,138-147,159-169,180-192"

def intervalli(s):
    res = set()
    for parte in s.split(","):
        a, _, b = parte.partition("-")
        res.update(range(int(a), int(b or a) + 1))
    return res

args = sys.argv[1:]
if "--caldi" in args:
    i = args.index("--caldi"); CALDI = args[i + 1]; del args[i:i + 2]
IN, OUT = args
caldi = intervalli(CALDI)

righe, sezione, molecola, prima = [], None, None, None
n_caldi, q_caldi, res_visti = 0, 0.0, set()
for l in open(IN):
    s = l.split(";")[0].split()
    if s[:1] == ["["]:
        sezione = s[1]
    elif sezione == "moleculetype" and s and molecola is None:
        molecola = prima = s[0]
    elif sezione == "moleculetype" and s:
        molecola = s[0]
    elif sezione == "atoms" and molecola == prima and len(s) >= 8:
        if "_" in s[1]:
            raise SystemExit(f"il tipo {s[1]} contiene gia' '_': il file e' gia' marcato?")
        if int(s[2]) in caldi:
            commento = l[l.index(";"):].rstrip("\n") if ";" in l else ""
            s[1] += "_"
            l = f"{s[0]:>6} {s[1]:>10} {s[2]:>6} {s[3]:>6} {s[4]:>6} {s[5]:>6} {' '.join(f'{x:>10}' for x in s[6:])}   {commento}\n"
            n_caldi += 1; q_caldi += float(s[6]); res_visti.add(int(s[2]))
    righe.append(l)

if not n_caldi:
    raise SystemExit(f"nessun atomo marcato: i residui {CALDI} non sono nella molecola {prima}")
mancanti = sorted(caldi - res_visti)
open(OUT, "w").writelines(righe)
print(f"molecola {prima}: {n_caldi} atomi caldi in {len(res_visti)} residui ({CALDI}), carica {q_caldi:+.3f}")
if mancanti:
    print(f"ATTENZIONE: residui richiesti ma assenti nella topologia: {mancanti}")
