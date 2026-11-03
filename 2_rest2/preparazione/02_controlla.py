"""Controlli sulla struttura NMR mutata prima di pdb2gmx.

Uso: python rest2/02_controlla.py <nmr_originale.pdb> <nmr_mutato.pdb> <xray/conf.gro>
Stampa: sequenza, atomi mancanti, atomi spostati, contatti stretti, ponte S-S, chi1 dei residui
innestati, carica netta e una minimizzazione di prova con amber14 (OpenMM, solvente implicito).
"""
import sys, warnings, numpy as np
warnings.filterwarnings("ignore")
ORIG, MUT, XRAY = sys.argv[1:4]
MUTATI = [86, 106, 128, 143, 145, 160, 167, 182]

def leggi_pdb(path):
    A = []
    for l in open(path):
        if l.startswith("ATOM"):
            A.append((int(l[22:26]), l[17:20].strip(), l[12:16].strip(),
                      np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return A

def leggi_gro(path):
    r = open(path).read().splitlines()
    return [(int(l[0:5]), l[5:10].strip(), l[10:15].strip(), 10 * np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])]))
            for l in r[2:2 + int(r[1])]]

def diedro(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 /= np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))

o, mu, x = leggi_pdb(ORIG), leggi_pdb(MUT), [a for a in leggi_gro(XRAY) if not a[2].startswith("H")]
seq = lambda A: {k: n for k, n, _, _ in A}
so, sm, sx = seq(o), seq(mu), seq(x)
alias = {"CYX": "CYS"}
diff = [k for k in sx if alias.get(sm[k], sm[k]) != alias.get(sx[k], sx[k])]
print(f"1. Sequenza: {len(sm)} residui ({min(sm)}-{max(sm)}); differenze dalla sequenza X-ray: {diff or 'nessuna'}")

term = {"O1": "O", "O2": "OXT"}
nomi = lambda A, k: {term.get(a, a) if k == 202 else a for kk, _, a, _ in A if kk == k}
mancanti = {k: (nomi(x, k) - nomi(mu, k), nomi(mu, k) - nomi(x, k)) for k in sx}
mancanti = {k: v for k, v in mancanti.items() if v[0] or v[1]}
print(f"2. Atomi pesanti rispetto all'X-ray (mancanti, in piu'): {mancanti or 'identici in tutti i residui'}")

po = {(k, a): p for k, _, a, p in o}
spost = [np.linalg.norm(p - po[(k, a)]) for k, _, a, p in mu if (k, a) in po]
nuovi = [(k, a) for k, _, a, p in mu if (k, a) not in po or k in (145, 167, 182) and a not in ("N", "CA", "C", "O")]
print(f"3. Atomi conservati dal file NMR: {len(spost)}, spostamento massimo {max(spost):.3f} A; atomi nuovi o innestati: {len(nuovi)}")

def contatti(A, soglia=3.0):
    P = np.array([a[3] for a in A]); R = np.array([a[0] for a in A]); N = [a[2] for a in A]
    d = np.linalg.norm(P[:, None] - P[None], axis=2)
    out = []
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            if R[j] - R[i] <= 1:          # stesso residuo o legame peptidico
                continue
            if {N[i], N[j]} == {"SG"}:    # ponte disolfuro
                continue
            if d[i, j] < soglia:
                out.append((R[i], N[i], R[j], N[j], d[i, j]))
    return out

co, cm = contatti(o), contatti(mu)
cmm = [c for c in cm if c[0] in MUTATI or c[2] in MUTATI]
print(f"4. Contatti fra atomi pesanti < 3.0 A (residui non adiacenti): originale {len(co)}, mutata {len(cm)}; "
      f"con un residuo mutato: {len(cmm)}")
for c in sorted(cmm, key=lambda c: c[4]):
    tipo = "legame H/ponte salino possibile" if c[1][0] in "NO" and c[3][0] in "NO" else "contatto sterico"
    print(f"      {c[0]}{c[1]} - {c[2]}{c[3]}  {c[4]:.2f} A  ({tipo})")
print(f"   distanza minima fra atomi pesanti nella mutata: {min(c[4] for c in cm):.2f} A (originale {min(c[4] for c in co):.2f} A)")

sg = [p for k, _, a, p in mu if a == "SG"]
print(f"5. SG43-SG87: {np.linalg.norm(sg[0] - sg[1]):.3f} A; nomi dei residui 43/87: {sm[43]}/{sm[87]}")

get = lambda A, k, a: next(p for kk, _, aa, p in A if kk == k and aa == a)
for k, g in ((145, "CG1"), (167, "CG"), (182, "OG1")):
    c1m = diedro(*(get(mu, k, a) for a in ("N", "CA", "CB", g)))
    c1x = diedro(*(get(x, k, a) for a in ("N", "CA", "CB", g)))
    phi_m = diedro(get(mu, k - 1, "C"), get(mu, k, "N"), get(mu, k, "CA"), get(mu, k, "C"))
    phi_x = diedro(get(x, k - 1, "C"), get(x, k, "N"), get(x, k, "CA"), get(x, k, "C"))
    print(f"6. {k}{sm[k]}: chi1 {c1m:7.1f} (X-ray {c1x:7.1f}); phi {phi_m:7.1f} (X-ray {phi_x:7.1f})")

# minimizzazione di prova con OpenMM (non sostituisce pdb2gmx: serve solo a vedere se la struttura regge)
import openmm as mm, openmm.app as app, openmm.unit as u
def con_oxt(path):
    """Copia temporanea con l'OXT del C-terminale (pdb2gmx lo aggiunge da se', OpenMM no)."""
    righe = [l for l in open(path) if l.startswith("ATOM")]
    ultimo = int(righe[-1][22:26])
    p = {l[12:16].strip(): np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) for l in righe if int(l[22:26]) == ultimo}
    uo, uca = [(p[n] - p["C"]) / np.linalg.norm(p[n] - p["C"]) for n in ("O", "CA")]
    oxt = p["C"] - 1.25 * (uo + uca) / np.linalg.norm(uo + uca)
    l = righe[-1]
    righe.append(f"{l[:6]}{int(l[6:11]) + 1:5d}  OXT{l[16:30]}{oxt[0]:8.3f}{oxt[1]:8.3f}{oxt[2]:8.3f}{l[54:76]} O\n")
    tmp = path + ".oxt.pdb"
    open(tmp, "w").writelines(righe + ["TER\n", "END\n"])
    return tmp

def prova(path):
    pdb = app.PDBFile(con_oxt(path))
    for r in pdb.topology.residues():
        if r.name == "CYX":
            r.name = "CYS"
    top = pdb.topology
    top.createDisulfideBonds(pdb.positions)
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    mod = app.Modeller(top, pdb.positions); mod.addHydrogens(ff, pH=7.0)
    s = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
    q = sum(s.getForce(i).getParticleParameters(j)[0].value_in_unit(u.elementary_charge)
            for i in range(s.getNumForces()) if isinstance(s.getForce(i), mm.NonbondedForce)
            for j in range(s.getNumParticles()))
    ctx = mm.Context(s, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("CPU"), {"Threads": "4"})
    ctx.setPositions(mod.positions)
    e0 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    mm.LocalEnergyMinimizer.minimize(ctx, 100.0, 2000)
    st = ctx.getState(getEnergy=True, getPositions=True)
    e1 = st.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    P0 = np.array(mod.positions.value_in_unit(u.angstrom)); P1 = st.getPositions(asNumpy=True).value_in_unit(u.angstrom)
    pes = [a.index for a in mod.topology.atoms() if a.element.symbol != "H"]
    res = np.array([int(a.residue.id) for a in mod.topology.atoms()])[pes]
    d = np.linalg.norm(P1[pes] - P0[pes], axis=1)
    nss = sum(1 for b in mod.topology.bonds() if b[0].name == "SG" and b[1].name == "SG")
    return q, e0, e1, d, res, nss
for nome, path in (("originale", ORIG), ("mutata", MUT)):
    q, e0, e1, d, res, nss = prova(path)
    peggio = sorted({(r, round(float(d[res == r].max()), 2)) for r in set(res)}, key=lambda t: -t[1])[:4]
    dm = {k: round(float(d[res == k].max()), 2) for k in MUTATI}
    print(f"7. OpenMM amber14 + GBn2, {nome}: carica {q:+.0f}, ponti S-S {nss}, energia {e0:.0f} -> {e1:.0f} kJ/mol, "
          f"spostamento atomi pesanti medio {d.mean():.2f} A max {d.max():.2f} A")
    print(f"      residui piu' spostati {peggio}; residui mutati {dm}")
