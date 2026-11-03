#!/usr/bin/env python3
"""Confronta i termini legati e 1-4 di GROMACS (amber19sb.ff) con OpenMM (amber19-all.xml = ff19SB) sulla stessa struttura.

Uso: python3 04b_confronta_openmm.py <vuoto.log> <vuoto.pdb> <uscita.txt>

Criterio del revisore (31 ott 2026): scarto relativo < 1e-4 per ogni termine; se un termine non torna, non si parte.
Termini: Bond, Angle, diedri (propri + impropri insieme, OpenMM li tiene in una sola forza), CMAP, LJ-14, Coulomb-14.
"""
import sys, re
import openmm as mm
from openmm import app, unit

AA = {"ALA","ARG","ASN","ASP","ASH","CYS","CYX","GLN","GLU","GLH","GLY","HID","HIE","HIP","HIS","ILE","LEU","LYS","LYN",
      "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}


def energie_gromacs(log):
    righe = open(log).read().splitlines()
    i = next(k for k, l in enumerate(righe) if "Energies (kJ/mol)" in l)
    e = {}
    k = i + 1
    while k + 1 < len(righe) and righe[k].strip():
        nomi = [righe[k][j:j + 15].strip() for j in range(0, len(righe[k]), 15)]
        val = [righe[k + 1][j:j + 15].strip() for j in range(0, len(righe[k + 1]), 15)]
        e.update({n: float(v) for n, v in zip(nomi, val) if n})
        k += 2
    return e


def carica_pdb(p):
    """pdb di editconf -> pdb per OpenMM: toglie i prefissi N/C dei terminali e, poiche' con -merge all editconf non scrive
    la catena, apre una catena nuova (A, B, C...) a ogni residuo N-terminale, cosi' OpenMM non lega fra loro catene diverse."""
    testo, catene, ultimo = [], -1, None
    for l in open(p):
        if l.startswith(("ATOM", "HETATM")):
            rn = l[17:21].strip()
            res = (l[22:27], rn)
            if res != ultimo and len(rn) == 4 and rn[0] == "N" and rn[1:] in AA:
                catene += 1
            ultimo = res
            if len(rn) == 4 and rn[0] in "NC" and rn[1:] in AA:
                rn = rn[1:]
            an = {"OC1": "O", "OC2": "OXT"}.get(l[12:16].strip(), l[12:16].strip())
            an = f" {an:<3s}" if len(an) < 4 else an
            l = l[:12] + an + l[16:17] + f"{rn:>3s} " + "ABCDEFGHIJ"[max(catene, 0)] + l[22:]
        testo.append(l)
    tmp = p + ".openmm.pdb"
    open(tmp, "w").writelines(testo)
    return app.PDBFile(tmp)


def main(log, pdb, out):
    eg = energie_gromacs(log)
    pdbf = carica_pdb(pdb)
    ff = app.ForceField("amber19-all.xml")
    sistema = ff.createSystem(pdbf.topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    gruppi = {}
    for i, f in enumerate(sistema.getForces()):
        f.setForceGroup(i); gruppi[type(f).__name__] = i
    nb = sistema.getForces()[gruppi["NonbondedForce"]]
    for i in range(nb.getNumParticles()):
        q, s, e = nb.getParticleParameters(i); nb.setParticleParameters(i, 0, s, 0)   # solo eccezioni = coppie 1-4
    eccez = [nb.getExceptionParameters(i) for i in range(nb.getNumExceptions())]

    def energia(gruppo):
        ctx = mm.Context(sistema, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
        ctx.setPositions(pdbf.positions)
        return ctx.getState(getEnergy=True, groups={gruppo}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    eo = {"Bond": energia(gruppi["HarmonicBondForce"]), "Angle": energia(gruppi["HarmonicAngleForce"]),
          "Diedri": energia(gruppi["PeriodicTorsionForce"]), "CMAP": energia(gruppi["CMAPTorsionForce"])}
    for i, (a, b, qq, s, e) in enumerate(eccez):
        nb.setExceptionParameters(i, a, b, 0, s, e)
    eo["LJ-14"] = energia(gruppi["NonbondedForce"])
    for i, (a, b, qq, s, e) in enumerate(eccez):
        nb.setExceptionParameters(i, a, b, qq, s, 0)
    eo["Coulomb-14"] = energia(gruppi["NonbondedForce"])

    diedri_g = sum(v for k, v in eg.items() if "Dih." in k and "CMAP" not in k)
    cmap_g = sum(v for k, v in eg.items() if "CMAP" in k)
    confronto = [("Bond", eg.get("Bond")), ("Angle", eg.get("Angle")), ("Diedri", diedri_g), ("CMAP", cmap_g),
                 ("LJ-14", eg.get("LJ-14")), ("Coulomb-14", eg.get("Coulomb-14"))]
    righe, ok = [f"termini GROMACS letti: {sorted(eg)}", f"{'termine':<12}{'GROMACS':>16}{'OpenMM':>16}{'scarto rel':>12}"], True
    for nome, g in confronto:
        o = eo[nome]
        rel = abs(g - o) / max(abs(o), 1e-9) if g is not None else float("inf")
        ok &= rel < 1e-4
        righe.append(f"{nome:<12}{g if g is not None else float('nan'):16.4f}{o:16.4f}{rel:12.2e}{'' if rel < 1e-4 else '  <-- NON TORNA'}")
    righe.append("ESITO: " + ("tutti i termini entro 1e-4" if ok else "NON si parte: almeno un termine fuori tolleranza"))
    righe.append("Non confrontati: LJ e Coulomb a distanza (dipendono dal cutoff).")
    open(out, "w").write("\n".join(righe) + "\n")
    print("\n".join(righe))


if __name__ == "__main__":
    main(*sys.argv[1:4])
