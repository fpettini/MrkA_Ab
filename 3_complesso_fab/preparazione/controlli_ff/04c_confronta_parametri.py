#!/usr/bin/env python3
"""Confronto a livello di parametri GROMACS (amber19sb.ff) contro OpenMM (amber19-all.xml), richiesto dal revisore
il 31 ott 2026 per coprire cio' che il controllo energie nel vuoto (04b) non vede: non legati a distanza, acqua, ioni.

Uso: python3 04c_confronta_parametri.py <libero/topol.top> <libero/vuoto.pdb> <cartella ff GROMACS> <uscita.txt>
  cartella ff GROMACS = copia di $GMXDATA/top/amber19sb.ff (servono ffnonbonded.itp, opc.itp, ions_opc.itp)

1. Proteina: per ogni atomo carica, sigma ed epsilon (topol.top + ffnonbonded.itp) contro la NonbondedForce di OpenMM;
   scarto massimo e carica totale per residuo.
2. OPC: cariche, sigma/epsilon dell'O, geometria (O-H, H-H) e posizione del sito M (opc.itp) contro amber19/opc.xml.
3. K+ e Cl-: sigma/epsilon di ions_opc.itp contro il set OPC di OpenMM (Li/Merz 12-6, frcmod.ionslm_126_opc).
Tolleranze: carica 1e-4 e, sigma 1e-5 nm, epsilon relativo 1e-4.
"""
import sys, os, re, math
import xml.etree.ElementTree as ET
import openmm as mm
from openmm import app, unit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module
carica_pdb = import_module("04b_confronta_openmm").carica_pdb

TQ, TS, TE = 1e-4, 1e-5, 1e-4


def sezioni(path):
    """[(sezione, campi)] di un file .itp/.top, senza commenti, ignorando le direttive del preprocessore"""
    sez, out = None, []
    for r in open(path):
        r = r.split(";")[0].strip()
        if not r or r.startswith("#"):
            continue
        m = re.match(r"\[\s*(\S+)\s*\]", r)
        if m:
            sez = m.group(1); continue
        out.append((sez, r.split()))
    return out


def tipi_gromacs(*files):
    """atomtypes: nome -> (sigma nm, epsilon kJ/mol)"""
    t = {}
    for f in files:
        if not os.path.exists(f):
            continue
        for s, c in sezioni(f):
            if s == "atomtypes":
                # nome [at.num] massa carica ptype sigma epsilon: sigma ed epsilon sono gli ultimi due campi
                t[c[0]] = (float(c[-2]), float(c[-1]))
    return t


def main(top, pdb, ffdir, out):
    righe, ok = [], True
    tipi = tipi_gromacs(os.path.join(ffdir, "ffnonbonded.itp"), os.path.join(ffdir, "ions_opc.itp"),
                        os.path.join(ffdir, "opc.itp"))
    # --- 1. proteina ---
    atomi = [c for s, c in sezioni(top) if s == "atoms"]
    pdbf = carica_pdb(pdb)
    ff = app.ForceField("amber19-all.xml")
    sistema = ff.createSystem(pdbf.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    nb = [f for f in sistema.getForces() if isinstance(f, mm.NonbondedForce)][0]
    oatomi = list(pdbf.topology.atoms())
    assert len(atomi) == len(oatomi) == nb.getNumParticles(), f"atomi: top {len(atomi)}, pdb {len(oatomi)}"
    dq = ds = de = 0.0; peggio = {}; nomi_diversi = 0; senza_tipo = set()
    qres_g, qres_o = {}, {}; nres, prec = 0, None
    for c, a in zip(atomi, oatomi):
        tipo, resnr, resn, nome, q = c[1], int(c[2]), c[3], c[4], float(c[6])
        nres += (resnr, resn) != prec; prec = (resnr, resn)   # con -merge all le catene riusano i numeri: conto i residui in ordine
        qo, so, eo = nb.getParticleParameters(a.index)
        qo = qo.value_in_unit(unit.elementary_charge); so = so.value_in_unit(unit.nanometer)
        eo = eo.value_in_unit(unit.kilojoule_per_mole)
        nome_o = {"OC1": "O", "OC2": "OXT"}.get(nome, nome)
        if resn.endswith("ILE") and nome == "CD":
            nome_o = "CD1"                     # stesso atomo: GROMACS chiama CD il CD1 dell'isoleucina
        nomi_diversi += nome_o != a.name and not (nome_o.startswith("H") and a.name.startswith("H"))
        qres_g[(nres, resnr, resn)] = qres_g.get((nres, resnr, resn), 0) + q
        qres_o[(nres, resnr, resn)] = qres_o.get((nres, resnr, resn), 0) + qo
        if abs(q - qo) > dq:
            dq = abs(q - qo); peggio["carica"] = f"{resn}{resnr} {nome}: {q:+.5f} contro {qo:+.6f}"
        if tipo not in tipi:
            senza_tipo.add(tipo); continue
        sg, eg = tipi[tipo]
        if abs(sg - so) > ds and eo > 0:
            ds = abs(sg - so); peggio["sigma"] = f"{resn}{resnr} {nome} ({tipo}): {sg:.6f} contro {so:.6f} nm"
        rel = abs(eg - eo) / max(eo, 1e-12) if eo > 0 or eg > 0 else 0.0
        if rel > de:
            de = rel; peggio["epsilon"] = f"{resn}{resnr} {nome} ({tipo}): {eg:.6f} contro {eo:.6f} kJ/mol"
    dres = max(abs(qres_g[k] - qres_o[k]) for k in qres_g)
    nonintere = [f"{k[2]}{k[1]} {v:+.4f}" for k, v in qres_g.items() if abs(v - round(v)) > 1e-3]
    righe += ["== 1. Proteina (%s: %d atomi, %d residui, %d catene in OpenMM)" % (top, len(atomi), len(qres_g), pdbf.topology.getNumChains()),
              f"carica: scarto massimo per atomo {dq:.2e} e  ({peggio.get('carica', '-')})",
              f"sigma: scarto massimo {ds:.2e} nm  ({peggio.get('sigma', '-')})",
              f"epsilon: scarto relativo massimo {de:.2e}  ({peggio.get('epsilon', '-')})",
              f"carica per residuo: scarto massimo {dres:.2e} e; carica totale GROMACS {sum(qres_g.values()):+.4f}, OpenMM {sum(qres_o.values()):+.4f}",
              f"residui con carica GROMACS non intera (>1e-3): {len(nonintere)} {' '.join(nonintere[:10])}",
              f"nomi d'atomo diversi nello stesso ordine (idrogeni esclusi): {nomi_diversi}"]
    if senza_tipo:
        righe.append(f"ATTENZIONE: tipi senza sigma/epsilon in ffnonbonded.itp: {sorted(senza_tipo)}"); ok = False
    ok &= dq < TQ and ds < TS and de < TE and nomi_diversi == 0

    # --- 2. OPC e 3. ioni, dal file xml di OpenMM ---
    xml = os.path.join(os.path.dirname(app.__file__), "data", "amber19", "opc.xml")
    radice = ET.parse(xml).getroot()
    nbx = {e.get("type"): e for e in radice.iter("Atom") if e.get("sigma") is not None}
    resx = {r.get("name"): r for r in radice.iter("Residue")}
    righe.append("== 2. OPC (opc.itp contro amber19/opc.xml)")
    w = resx.get("HOH")
    qx = {a.get("name"): (a.get("type"),) for a in w.iter("Atom")}
    q_o = {e.get("type"): float(e.get("charge")) for e in radice.iter("Atom") if e.get("charge") is not None}
    sez = sezioni(os.path.join(ffdir, "opc.itp"))
    at = [c for s, c in sez if s == "atoms"]
    mappa = {"OW": "O", "HW1": "H1", "HW2": "H2", "MW": "M"}
    for c in at:
        nome = c[4]; tx = qx[mappa[nome]][0]; qg = float(c[6])
        d = abs(qg - q_o[tx]); ok &= d < TQ
        righe.append(f"carica {nome}: {qg:+.6f} contro {q_o[tx]:+.6f}  scarto {d:.1e}")
        if c[1] in tipi and float(nbx[tx].get("epsilon")) > 0:
            sg, eg = tipi[c[1]]; sx, ex = float(nbx[tx].get("sigma")), float(nbx[tx].get("epsilon"))
            okl = abs(sg - sx) < TS and abs(eg - ex) / ex < TE; ok &= okl
            righe.append(f"LJ {nome} ({c[1]}): sigma {sg:.6f}/{sx:.6f} nm, epsilon {eg:.6f}/{ex:.6f} kJ/mol {'ok' if okl else 'NON TORNA'}")
    sett = [c for s, c in sez if s == "settles"]
    vs = [c for s, c in sez if s and s.startswith("virtual_sites") or s == "dummies3"]
    if sett:
        doh, dhh = float(sett[0][2]), float(sett[0][3])
        cx = {frozenset((b.get("atomName1"), b.get("atomName2"))): float(b.get("length")) for b in radice.iter("Bond") if b.get("length")}
        righe.append(f"settles: O-H {doh:.8f} nm, H-H {dhh:.8f} nm")
        # vincoli di OpenMM (rigidWater): O-H dal legame, H-H dall'angolo HW-OW-HW del xml
        lx = [float(b.get("length")) for b in radice.iter("Bond") if b.get("class1") and b.get("length")]
        ax = [float(a.get("angle")) for a in radice.iter("Angle") if a.get("class1") and a.get("angle")]
        if lx and ax:
            hh = 2 * lx[0] * math.sin(ax[0] / 2)
            okg = abs(doh - lx[0]) < 1e-7 and abs(dhh - hh) < 1e-7; ok &= okg
            righe.append(f"  opc.xml: O-H {lx[0]:.8f} nm, H-H {hh:.8f} nm (da angolo {math.degrees(ax[0]):.4f} gradi) {'ok' if okg else 'NON TORNA'}")
    if vs:
        righe.append(f"sito virtuale (opc.itp): {' '.join(vs[0])}")
    for v in w.iter("VirtualSite"):
        righe.append(f"  opc.xml VirtualSite: {v.attrib}")
        if vs:
            a_g = float(vs[0][-2]); b_g = float(vs[0][-1])
            wx = [float(v.get(k)) for k in ("weight1", "weight2", "weight3")]
            dw = max(abs(1 - a_g - b_g - wx[0]), abs(a_g - wx[1]), abs(b_g - wx[2]))
            # spostamento di M: dw * |(H1-O) + (H2-O)| = dw * 2 * dOH * cos(HOH/2)
            dm = dw * 2 * 0.08724331 * math.cos(math.radians(103.6 / 2))
            righe.append(f"  pesi di M: GROMACS {1 - a_g - b_g:.8f}/{a_g:.8f}/{b_g:.8f}, OpenMM {wx[0]:.8f}/{wx[1]:.8f}/{wx[2]:.8f};"
                         f" scarto {dw:.1e}, spostamento di M {dm:.1e} nm")
            ok &= dm < 1e-5
    righe.append("== 3. Ioni K+ e Cl- (ions_opc.itp contro ionslm_126_opc di OpenMM)")
    for nome_g, nome_x in (("K+_OPC", "ionslm_126_opc-K+"), ("CL-_OPC", "ionslm_126_opc-Cl-")):
        cand = [k for k in tipi if k.upper() == nome_g]
        e = nbx.get(nome_x)
        if not cand or e is None:
            righe.append(f"{nome_g}: tipo non trovato (GROMACS {cand}, OpenMM {nome_x in nbx})"); ok = False; continue
        sx, ex = float(e.get("sigma")), float(e.get("epsilon"))
        for k in cand:
            sg, eg = tipi[k]; okl = abs(sg - sx) < TS and abs(eg - ex) / ex < TE; ok &= okl
            righe.append(f"{k}: sigma {sg:.6f}/{sx:.6f} nm, epsilon {eg:.6f}/{ex:.6f} kJ/mol {'ok' if okl else 'NON TORNA'}")
    righe.append("ESITO: " + ("tutto entro le tolleranze" if ok else "ALMENO UN PARAMETRO NON TORNA: vedi sopra"))
    open(out, "w").write("\n".join(righe) + "\n")
    print("\n".join(righe))


if __name__ == "__main__":
    main(*sys.argv[1:5])
