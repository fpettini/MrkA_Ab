#!/usr/bin/env python3
"""Aiuto per 02_prepara_sistemi.sh (solo libreria standard, gira sulla workstation).

  dividi  <mrka-fab_completo.pdb> <cartella> [acque_interne.txt]
          scrive complesso/, libero/, fv/ con proteina.pdb, acque.txt (7 a ponte + interne; libero: solo interne), ss_attesi.txt
  acque   <proteina.gro> <acque.txt> <topol.top> <conf.gro>   aggiunge le 7 acque cristallografiche (TIP3P) dopo la proteina
  controlla <solv.gro> <ions.gro> <n>  esce con 1 se genion ha sostituito una delle n acque cristallografiche
  gruppo  <gro> <index.ndx>  aggiunge il gruppo Protein_Other (make_ndx lascia fuori il C-terminale CCYX) e verifica la copertura
  indice  <solv.gro> <n_acque_cristallo> <solvente.ndx>       gruppo SOL_bulk = acqua di solvatazione (genion non tocca le cristallografiche)
"""
import sys, os, math

FV = {"B": 107, "C": 117}          # fine di VL e VH (VEIK107, VTVSS117)
# Istidine (tutte nel Fab, MrkA non ne ha): neutre perche' PROPKA da' pKa 5.2-6.95 (31 ott 2026). Tautomero dal partner del legame H
# (proposta da approvare dal revisore). Scritte esplicite perche' l'assegnazione automatica di pdb2gmx (rete di legami H) faceva H204 = HISH (+1):
#   L189 HID (ND1-O Asp185 2.82 A); L198 HIE (NE2-O Pro141 2.91; ND1 accettore di N 199/200); H53 HIE (nessun partner);
#   H168 HIE (NE2-OD1 Asn138 L 3.02); H204 HIE (NE2-O Pro151 2.83; ND1 accettore di OG Ser207 2.70; PROPKA 5.95)
ISTIDINE = {("B", 189): "HID", ("B", 198): "HIE", ("C", 53): "HIE", ("C", 168): "HIE", ("C", 204): "HIE"}


def leggi_pdb(p):
    return [l for l in open(p) if l.startswith(("ATOM", "HETATM"))]


def dividi(pdb, out, interne=None):
    righe = leggi_pdb(pdb)
    his = {(l[21], int(l[22:26])) for l in righe if l[17:20] == "HIS"}
    assert his == set(ISTIDINE), f"istidine nel pdb {sorted(his)} diverse da ISTIDINE {sorted(ISTIDINE)}"
    righe = [l[:17] + ISTIDINE[(l[21], int(l[22:26]))] + l[20:] if l[17:20] == "HIS" else l for l in righe]
    prot = [l for l in righe if l[17:20] != "HOH"]
    acque = [l for l in righe if l[17:20] == "HOH"]
    assert len(acque) == 7, f"attese 7 acque, trovate {len(acque)}"
    ponte = [f"{float(l[30:38]) / 10:.4f} {float(l[38:46]) / 10:.4f} {float(l[46:54]) / 10:.4f}  # ponte {l[21]}{l[22:26].strip()}\n" for l in acque]
    dentro = open(interne).readlines() if interne else []
    sistemi = {
        "complesso": [l for l in prot if l[21] in "ABC"],
        "libero": [l for l in prot if l[21] == "A"],
        "fv": [l for l in prot if l[21] == "A" or (l[21] in FV and int(l[22:26]) <= FV[l[21]])],
    }
    for nome, rr in sistemi.items():
        d = os.path.join(out, nome); os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "proteina.pdb"), "w") as f:
            prec = None
            for l in rr:
                if prec and l[21] != prec:
                    f.write("TER\n")
                f.write(l); prec = l[21]
            f.write("TER\nEND\n")
        cys = {(l[21], l[22:26]) for l in rr if l[17:20] == "CYS"}
        open(os.path.join(d, "ss_attesi.txt"), "w").write(f"{len(cys) // 2}\n")
        cris = dentro if nome == "libero" else ponte + dentro
        if cris:
            open(os.path.join(d, "acque.txt"), "w").writelines(cris)
        print(f"{nome}: {len({(l[21], l[22:26]) for l in rr})} residui, Cys {len(cys)}, S-S attesi {len(cys) // 2}, acque cristallografiche {len(cris)}")


def leggi_gro(p):
    r = open(p).read().splitlines()
    return r[0], [x for x in r[2:2 + int(r[1])]], r[2 + int(r[1])]


def acque(gro, txt, top, out, modello="tip3p"):
    tit, at, box = leggi_gro(gro)
    xyz = [(float(a[20:28]), float(a[28:36]), float(a[36:44])) for a in at]
    pesanti = [p for a, p in zip(at, xyz) if not a[10:15].strip().startswith("H")]
    resnr = int(at[-1][0:5])
    nuove = []
    # geometria iniziale (SETTLE e i siti virtuali la riportano a quella del modello): TIP3P 0.9572 A/104.52; OPC 0.8724 A/103.6 + MW a 0.1594 A
    b, ang = (0.08724, math.radians(103.6)) if modello == "opc" else (0.09572, math.radians(104.52))
    for i, riga in enumerate(open(txt)):
        o = tuple(map(float, riga.split("#")[0].split()))
        vic = [p for p in pesanti if sum((p[k] - o[k]) ** 2 for k in range(3)) < 0.4 ** 2]
        c = [sum(p[k] for p in vic) / len(vic) for k in range(3)]
        u = [o[k] - c[k] for k in range(3)]; n = math.sqrt(sum(x * x for x in u)); u = [x / n for x in u]
        t = [1.0, 0.0, 0.0] if abs(u[0]) < 0.9 else [0.0, 1.0, 0.0]          # vettore perpendicolare a u
        w = [u[1] * t[2] - u[2] * t[1], u[2] * t[0] - u[0] * t[2], u[0] * t[1] - u[1] * t[0]]
        n = math.sqrt(sum(x * x for x in w)); w = [x / n for x in w]
        # bisettrice lungo u (verso l'esterno della proteina), idrogeni a +-52.26 gradi
        h1 = [o[k] + b * (math.cos(ang / 2) * u[k] + math.sin(ang / 2) * w[k]) for k in range(3)]
        h2 = [o[k] + b * (math.cos(ang / 2) * u[k] - math.sin(ang / 2) * w[k]) for k in range(3)]
        r = resnr + i + 1
        siti = [("OW", o), ("HW1", h1), ("HW2", h2)]
        if modello == "opc":
            siti.append(("MW", [o[k] + 0.01594 * u[k] for k in range(3)]))
        for nome, p in siti:
            nuove.append(f"{r % 100000:5d}{'SOL':<5s}{nome:>5s}{(len(at) + len(nuove) + 1) % 100000:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}")
    with open(out, "w") as f:
        f.write(f"{tit}\n{len(at) + len(nuove)}\n" + "\n".join(at + nuove) + f"\n{box}\n")
    nsiti = 4 if modello == "opc" else 3
    t = open(top).read().rstrip("\n") + f"\nSOL                 {len(nuove) // nsiti}\n"
    open(top, "w").write(t)
    print(f"aggiunte {len(nuove) // nsiti} acque cristallografiche ({modello}) a {out}")


def controlla(prima, dopo, ncris):
    """esce con 1 se genion ha sostituito una delle acque cristallografiche (le prime ncris SOL)"""
    def ow(g):
        return [a[20:44] for a in leggi_gro(g)[1] if a[5:10].strip() == "SOL" and a[10:15].strip() == "OW"][:int(ncris)]
    if ow(prima) != ow(dopo):
        print("genion ha sostituito un'acqua cristallografica: riprovo con un altro seme"); sys.exit(1)
    print(f"acque cristallografiche intatte dopo genion ({ncris})")


def gruppo(gro, ndx):
    """aggiunge a index.ndx il gruppo Protein_Other (= Protein | Other di make_ndx, come "1 | 12") = tutti gli atomi che non sono SOL/ioni (compreso il C-terminale CCYX,
    che make_ndx NON mette in Protein) e verifica che Protein_Other + Water_and_ions copra tutto il sistema"""
    _, at, _ = leggi_gro(gro)
    solv = {"SOL", "NA", "K", "CL"}
    prot = [i for i, a in enumerate(at, 1) if a[5:10].strip() not in solv]
    testo = open(ndx).read()
    gruppi, nome = {}, None
    for r in testo.splitlines():
        if r.strip().startswith("["):
            nome = r.strip().strip("[] ").strip(); gruppi[nome] = []
        elif nome and r.strip():
            gruppi[nome] += [int(x) for x in r.split()]
    wi = set(gruppi["Water_and_ions"])
    assert set(prot).isdisjoint(wi) and len(prot) + len(wi) == len(at), \
        f"Protein_Other ({len(prot)}) + Water_and_ions ({len(wi)}) != atomi del sistema ({len(at)})"
    atteso = sorted(set(gruppi.get("Protein", [])) | set(gruppi.get("Other", [])))
    assert atteso == prot, "Protein | Other di make_ndx non coincide con gli atomi non solvente"
    if "Protein_Other" in gruppi:
        assert gruppi["Protein_Other"] == prot, "gruppo Protein_Other gia' presente e diverso"
        print("Protein_Other gia' presente"); return
    with open(ndx, "a") as f:
        f.write("\n[ Protein_Other ]\n")
        for k in range(0, len(prot), 15):
            f.write(" ".join(map(str, prot[k:k + 15])) + "\n")
    print(f"gruppo Protein_Other: {len(prot)} atomi (Protein di make_ndx: {len(gruppi.get('Protein', []))}); copertura completa")


def indice(gro, ncris, ndx):
    _, at, _ = leggi_gro(gro)
    sol, visti = [], 0
    prec = None
    for i, a in enumerate(at, 1):
        if a[5:10].strip() != "SOL":
            continue
        if a[10:15].strip() == "OW":
            visti += 1
        if visti > int(ncris):
            sol.append(i)
    with open(ndx, "w") as f:
        f.write("[ SOL_bulk ]\n")
        for k in range(0, len(sol), 15):
            f.write(" ".join(map(str, sol[k:k + 15])) + "\n")
    print(f"SOL_bulk: {len(sol) // 3} acque (escluse le prime {ncris} cristallografiche)")


if __name__ == "__main__":
    {"dividi": dividi, "acque": acque, "controlla": controlla, "gruppo": gruppo, "indice": indice}[sys.argv[1]](*sys.argv[2:])
