"""Densita' di OPC pura e OPC + KCl con i parametri di OpenMM amber19/opc.xml (= ions_opc/opc.itp di GROMACS, verificati in 04c),
impostazioni come il nostro mdp: PME, cutoff 1.2 nm, correzione di dispersione, 298 K, 1 bar, acqua rigida, 2 fs.
Uso: python box_opc.py <lato_nm> <n_coppie_KCl> <ps_equil> <ps_prod> <uscita.txt>"""
import sys, os, time, numpy as np
import openmm as mm, openmm.app as app, openmm.unit as u
lato, ncop, peq, ppr, out = float(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
ff = app.ForceField('amber19/opc.xml')
m = app.Modeller(app.Topology(), [])
m.addSolvent(ff, model='tip4pew', boxSize=mm.Vec3(lato, lato, lato) * u.nanometer, neutralize=False)
# sostituisce acque casuali con K e Cl (come genion)
rng = np.random.default_rng(7)
acque = [r for r in m.topology.residues() if r.name == 'HOH']
scelte = rng.choice(len(acque), 2 * ncop, replace=False)
pos = m.getPositions(); nuove_ioni = []
for k, i in enumerate(scelte):
    o = [a for a in acque[i].atoms() if a.element is not None and a.element.symbol == 'O'][0]
    nuove_ioni.append(('K' if k < ncop else 'CL', pos[o.index]))
m.delete([acque[i] for i in scelte])
for nome, p in nuove_ioni:
    t = app.Topology(); c = t.addChain(); r = t.addResidue(nome, c)
    t.addAtom(nome, app.Element.getBySymbol('K' if nome == 'K' else 'Cl'), r)
    m.add(t, [p])
nw = sum(1 for r in m.topology.residues() if r.name == 'HOH')
sis = ff.createSystem(m.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.2 * u.nanometer, constraints=app.HBonds, rigidWater=True)
nbf = [f for f in sis.getForces() if isinstance(f, mm.NonbondedForce)][0]
assert nbf.getUseDispersionCorrection()
sis.addForce(mm.MonteCarloBarostat(1 * u.bar, 298 * u.kelvin, 25))
massa = sum(sis.getParticleMass(i).value_in_unit(u.dalton) for i in range(sis.getNumParticles()))
mw = nw * (16.0 + 2 * 1.008)
integ = mm.LangevinMiddleIntegrator(298 * u.kelvin, 1 / u.picosecond, 0.002 * u.picoseconds)
piat = mm.Platform.getPlatformByName('CPU')
ctx = mm.Context(sis, integ, piat, {'Threads': os.environ.get('NTHREADS', '4')})
ctx.setPositions(m.getPositions()); ctx.setPeriodicBoxVectors(*m.topology.getPeriodicBoxVectors())
ctx.applyConstraints(1e-6); mm.LocalEnergyMinimizer.minimize(ctx, 10, 500)
ctx.setVelocitiesToTemperature(298 * u.kelvin, 3)
t0 = time.time(); integ.step(int(peq / 0.002)); veq = (time.time() - t0)
NA = 6.02214076e23; V = []
nblk = int(ppr / 0.2)
for i in range(nblk):
    integ.step(100)
    V.append(ctx.getState().getPeriodicBoxVolume().value_in_unit(u.nanometer ** 3))
V = np.array(V); rho = massa / NA / (V * 1e-21); rhow = mw / NA / (V * 1e-21)
# errore standard con medie a blocchi (10 blocchi)
b = rho[: len(rho) // 10 * 10].reshape(10, -1).mean(1); bw = rhow[: len(rhow) // 10 * 10].reshape(10, -1).mean(1)
mol = ncop / (nw * 18.016e-3)
riga = (f"lato {lato} nm, {nw} acque, {ncop} coppie KCl ({mol:.3f} mol/kg), equil {peq} ps ({peq/1000/(veq/86400):.1f} ns/g), prod {ppr} ps: "
        f"densita' soluzione {rho.mean():.4f} +- {b.std(ddof=1)/np.sqrt(10):.4f} g/cm3, massa d'acqua per volume {rhow.mean():.4f} +- {bw.std(ddof=1)/np.sqrt(10):.4f}")
print(riga); open(out, 'a').write(riga + "\n")
