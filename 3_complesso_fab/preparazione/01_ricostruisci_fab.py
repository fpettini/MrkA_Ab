#!/usr/bin/env python3
"""Ricostruisce i residui mancanti del Fab 03F20 nel complesso MrkA-Fab.

Uso: python3 01_ricostruisci_fab.py <mrka-fab_xray_full.pdb> <cartella_uscita>

Ingresso: pdb_last/mrka-fab_xray_full.pdb di francesco (A = MrkA wt 36-217 con loop 116-117 gia' ricostruito,
B = leggera 1-213, C = pesante 2-217 senza 9-10, E/S = 7 acque d'interfaccia).
Sequenze fornite da francesco il 31 ott 2026 (pesante fino a Gly223, tag His escluso).

Passi:
 1. tiene solo l'altloc A (Ser19 C: occupanze 0.50/0.50), toglie gli idrogeni (pdb2gmx -ignh li rigenera);
    le 7 acque (catene E, S) sono tolte durante la modellazione (stanno all'epitopo, >19 A dai residui nuovi)
    e rimesse alla fine con coordinate identiche;
 2. PDBFixer aggiunge C Gln1, Ala9, Gly10, Lys218-Gly223 e B Cys214;
 3. minimizzazione OpenMM (amber14 + GB-neck2) con TUTTI gli atomi gia' presenti congelati (massa 0) e un richiamo
    legame SG(B214)-SG(C220) del campo di forze (template CYX) per formare il disolfuro intercatena; solo minimizzazione (max 5000 iterazioni), cutoff non periodico 1.2 nm;
 4. scrive il modello senza idrogeni e un controllo (scarto degli atomi originali, SG-SG, contatti, omega, chiralita').
Gli omega dei legami peptidici che toccano un residuo nuovo sono tenuti trans con k(1+cos omega), k 500 kJ/mol.
"""
import sys, os, json, hashlib
import numpy as np
import gemmi
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer
from pdbfixer.pdbfixer import Sequence

PESANTE = ("QVQLQQWGAGPVKPSETLSLTCAVYDGSFTTYYWTWIRQPPGRGLEWIGEVMHSGSTFYNPSLQSRVTISVDTSKNQFSLKLASVTSADTAVYYCAR"
           "GTRERYFDIWGRGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGT"
           "QTYICNVNHKPSNTKVDKRVEPKSCDKG")
LEGGERA = ("AIQMTQSPSTLSASVGDTVTITCRASQSINSWLAWYQQKPGNVPKLLIQKASNLESGVPSRFSGSGSGTEFTLTINNLQPDDFATYDCQQYRTFPWTF"
           "GQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTH"
           "QGLSSPVTKSFNRGEC")
TRE = {v.one_letter_code.upper(): k for k, v in [(n, gemmi.find_tabulated_residue(n)) for n in
       ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]]}
PRIMO = {"A": 36, "B": 1, "C": 1}   # numero del primo residuo di ogni catena nella sequenza completa


BACKBONE = {"N", "H", "H1", "H2", "H3", "CA", "C", "O", "OXT"}


def chiralita(xyz):
    n, ca, c, cb = (xyz[k] for k in ("N", "CA", "C", "CB"))
    return float(np.dot(n - ca, np.cross(c - ca, cb - ca)))   # > 0 per gli L-amminoacidi (verificato sui 562 originali)


def correggi_chiralita(top, pos, residui):
    """Riflette la catena laterale (CB in poi e HA) sul piano N-CA-C dei residui D; ritorna i residui corretti."""
    pos = np.array(pos)
    corretti = []
    for r in top.residues():
        if (r.chain.id, int(r.id)) not in residui or r.name == "GLY":
            continue
        idx = {a.name: a.index for a in r.atoms()}
        if chiralita({k: pos[idx[k]] for k in ("N", "CA", "C", "CB")}) > 0:
            continue
        n, ca, c = pos[idx["N"]], pos[idx["CA"]], pos[idx["C"]]
        nrm = np.cross(n - ca, c - ca); nrm /= np.linalg.norm(nrm)
        for name, i in idx.items():
            if name not in BACKBONE:
                pos[i] = pos[i] - 2 * np.dot(pos[i] - ca, nrm) * nrm
        corretti.append(f"{r.chain.id}{r.id}{r.name}")
    return pos, corretti


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main(inp, out):
    os.makedirs(out, exist_ok=True)
    st = gemmi.read_structure(inp)
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    m = st[0]
    seqA = [r.name for r in m["A"]]
    assert [r.seqid.num for r in m["A"]] == list(range(36, 218)), "catena A non e' 36-217 continua"
    acque = gemmi.Structure(); acque.cell = st.cell; acque.add_model(gemmi.Model("1"))
    for c in m:
        if all(r.name == "HOH" for r in c):
            acque[0].add_chain(c)
    st_prot = st.clone()
    for n in [c.name for c in st_prot[0] if all(r.name == "HOH" for r in c)]:
        st_prot[0].remove_chain(n)
    pulito = os.path.join(out, "00_ingresso_pulito.pdb")
    st_prot.write_pdb(pulito)
    orig = {(c.name, r.seqid.num, a.name): np.array(a.pos.tolist()) / 10 for c in m for r in c for a in r}

    fx = PDBFixer(filename=pulito)
    fx.sequences = [Sequence("A", seqA), Sequence("B", [TRE[x] for x in LEGGERA]), Sequence("C", [TRE[x] for x in PESANTE])]
    fx.findMissingResidues()
    catene = list(fx.topology.chains())
    mancanti = {f"{catene[k[0]].id}:{k[1]}": v for k, v in fx.missingResidues.items()}
    print("residui mancanti trovati:", mancanti)
    attesi = {"B:213": ["CYS"], "C:0": ["GLN"], "C:7": ["ALA", "GLY"], "C:214": ["LYS", "SER", "CYS", "ASP", "LYS", "GLY"]}
    assert mancanti == attesi, "residui mancanti diversi da quelli attesi"
    fx.findNonstandardResidues(); assert not fx.nonstandardResidues
    fx.findMissingAtoms()
    print("atomi mancanti in residui esistenti:", {f"{r.chain.id}{r.id}{r.name}": [a.name for a in v] for r, v in fx.missingAtoms.items()})
    fx.addMissingAtoms(seed=20261031)
    fx.addMissingHydrogens(7.4)

    # rinumera: catena A da 36, B e C da 1, acque come in ingresso
    for ch in fx.topology.chains():
        if ch.id in PRIMO:
            for i, r in enumerate(ch.residues()):
                r.id = str(PRIMO[ch.id] + i)

    # residui nuovi = (catena, numero) assenti nell'ingresso
    origres = {(k[0], k[1]) for k in orig}
    nuovi = {(r.chain.id, int(r.id), r.name) for r in fx.topology.residues() if (r.chain.id, int(r.id)) not in origres}
    # disolfuro intercatena B Cys214 - C Cys220: tolgo HG e aggiungo il legame SG-SG (template CYX)
    md = app.Modeller(fx.topology, fx.positions)
    PONTE = {("B", 214), ("C", 220)}
    md.delete([a for a in md.topology.atoms() if a.name == "HG" and (a.residue.chain.id, int(a.residue.id)) in PONTE])
    sgat = {a.residue.chain.id: a for a in md.topology.atoms() if a.name == "SG" and (a.residue.chain.id, int(a.residue.id)) in PONTE}
    md.topology.addBond(sgat["B"], sgat["C"])
    top, pos = md.topology, md.positions
    sg = {k: v.index for k, v in sgat.items()}
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=1.2 * unit.nanometer, constraints=None)
    for a in top.atoms():
        if (a.residue.chain.id, int(a.residue.id), a.name) in orig:
            system.setParticleMass(a.index, 0)   # atomi pesanti del cristallo congelati
    # omega trans sui legami peptidici che toccano un residuo nuovo (PDBFixer puo' costruirli cis)
    om = mm.CustomTorsionForce("k_om*(1+cos(theta))"); om.addGlobalParameter("k_om", 500.0)
    for ch in top.chains():
        rs = [r for r in ch.residues() if r.name != "HOH"]
        for r1, r2 in zip(rs, rs[1:]):
            if (ch.id, int(r1.id), r1.name) in nuovi or (ch.id, int(r2.id), r2.name) in nuovi:
                at1 = {x.name: x.index for x in r1.atoms()}; at2 = {x.name: x.index for x in r2.atoms()}
                om.addTorsion(at1["CA"], at1["C"], at2["N"], at2["CA"], [])
    system.addForce(om)
    print("omega vincolati trans:", om.getNumTorsions())
    integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 5 / unit.picosecond, 0.5 * unit.femtosecond)
    integ.setRandomNumberSeed(20261031)
    sim = app.Simulation(top, system, integ, mm.Platform.getPlatformByName("CPU"))
    # PDBFixer costruisce D alcuni residui terminali aggiunti (verificato 31 ott): correggo prima e dopo la minimizzazione
    nuovi_rn = {(c, n) for c, n, _ in nuovi}
    p0, corr = correggi_chiralita(top, pos.value_in_unit(unit.nanometer), nuovi_rn)
    print("chiralita' corretta prima della minimizzazione:", corr)
    sim.context.setPositions(p0 * unit.nanometer)
    def stato(tag):
        s = sim.context.getState(getEnergy=True, getPositions=True)
        p = s.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        print(f"{tag}: E {s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.0f} kJ/mol, SG-SG {np.linalg.norm(p[sg['B']]-p[sg['C']])*10:.2f} A", flush=True)
        return p
    stato("dopo PDBFixer")
    # solo minimizzazione (scelta di francesco, 31 ott: il resto in GROMACS)
    sim.minimizeEnergy(maxIterations=5000); p = stato("minimizzato")
    for giro in range(3):
        p, corr = correggi_chiralita(top, p, nuovi_rn)
        print(f"chiralita' corretta dopo la minimizzazione (giro {giro + 1}):", corr)
        if not corr:
            break
        sim.context.setPositions(p * unit.nanometer)
        sim.minimizeEnergy(maxIterations=5000); p = stato(f"minimizzato di nuovo ({giro + 1})")
    else:
        sys.exit("chiralita' non stabile dopo 3 giri")

    # scrive senza idrogeni
    with open(os.path.join(out, "tmp_conH.pdb"), "w") as f:
        app.PDBFile.writeFile(top, p * unit.nanometer, f, keepIds=True)
    st2 = gemmi.read_structure(os.path.join(out, "tmp_conH.pdb")); st2.remove_hydrogens()
    os.remove(os.path.join(out, "tmp_conH.pdb"))
    for ch in st2[0]:
        for r in ch:
            if r.name == "CYS" and (ch.name, r.seqid.num) in {("B", 214), ("C", 220)}:
                pass
    for c in acque[0]:
        st2[0].add_chain(c)
    st2.setup_entities()
    fin = os.path.join(out, "mrka-fab_completo.pdb")
    st2.write_pdb(fin)

    # controlli
    m2 = st2[0]; ctl = {}
    dev = [np.linalg.norm(np.array(a.pos.tolist()) / 10 - orig[(c.name, r.seqid.num, a.name)]) * 10
           for c in m2 for r in c for a in r if (c.name, r.seqid.num, a.name) in orig]
    ctl["atomi_originali"] = len(dev); ctl["scarto_max_originali_A"] = round(max(dev), 4)
    ctl["residui_nuovi"] = sorted({f"{c}{n}{rn}" for c, n, rn in nuovi if rn != "HOH"})
    sgB = m2["B"]["214"][0]["SG"][0].pos; sgC = m2["C"]["220"][0]["SG"][0].pos
    ctl["SG_B214_SG_C220_A"] = round(sgB.dist(sgC), 3)
    ns = gemmi.NeighborSearch(m2, gemmi.UnitCell(), 5).populate()
    corti = []
    for c in m2:
        for r in c:
            if (c.name, r.seqid.num) not in {(x[0], x[1]) for x in nuovi}: continue
            for a in r:
                for mk in ns.find_atoms(a.pos, "\0", radius=3.0):
                    cra = mk.to_cra(m2)
                    if cra.chain.name == c.name and abs(cra.residue.seqid.num - r.seqid.num) <= 1: continue
                    d = cra.atom.pos.dist(a.pos)
                    if d < 3.0 and not ({a.name, cra.atom.name} == {"SG"}):
                        corti.append(f"{c.name}{r.seqid.num}{r.name}:{a.name}-{cra.chain.name}{cra.residue.seqid.num}{cra.residue.name}:{cra.atom.name} {d:.2f}")
    ctl["contatti_<3A_atomi_nuovi"] = sorted(set(corti))
    om = []
    for c in m2:
        rs = [r for r in c if r.name != "HOH"]
        for r1, r2 in zip(rs, rs[1:]):
            if r2.seqid.num != r1.seqid.num + 1: continue
            w = np.degrees(gemmi.calculate_dihedral(r1["CA"][0].pos, r1["C"][0].pos, r2["N"][0].pos, r2["CA"][0].pos))
            cn = r1["C"][0].pos.dist(r2["N"][0].pos)
            if abs(abs(w) - 180) > 30 or not 1.25 < cn < 1.40:
                om.append(f"{c.name}{r1.seqid.num}-{r2.seqid.num} omega {w:.0f} C-N {cn:.2f}")
    ctl["legami_peptidici_anomali"] = om
    def chir(r):
        return chiralita({x: np.array(r[x][0].pos.tolist()) for x in ("N", "CA", "C", "CB")})
    segno_orig = [np.sign(chir(r)) for c in m2 for r in c if r.name not in ("GLY", "HOH") and (c.name, r.seqid.num) in origres]
    ctl["chiralita_originali_segni"] = {str(k): int(v) for k, v in zip(*np.unique(segno_orig, return_counts=True))}
    ctl["chiralita_nuovi_segno"] = {f"{c.name}{r.seqid.num}{r.name}": int(np.sign(chir(r))) for c in m2 for r in c
                                    if r.name not in ("GLY", "HOH") and (c.name, r.seqid.num) not in origres}
    seqs = {c.name: gemmi.one_letter_code([r.name for r in c if r.name != "HOH"]) for c in m2 if c.name in "ABC"}
    ctl["sequenza_B_uguale"] = seqs["B"] == LEGGERA; ctl["sequenza_C_uguale"] = seqs["C"] == PESANTE
    ctl["numerazione"] = {c.name: [c[0].seqid.num, c[len(c) - 1].seqid.num, len(c)] for c in m2}
    ctl["sha256_ingresso"] = sha(inp); ctl["sha256_uscita"] = sha(fin)
    ctl["versioni"] = {"openmm": mm.__version__, "gemmi": gemmi.__version__}
    json.dump(ctl, open(os.path.join(out, "controllo_ricostruzione.json"), "w"), indent=1)
    print(json.dumps(ctl, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:3])
