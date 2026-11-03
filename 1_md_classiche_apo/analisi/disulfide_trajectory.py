#!/usr/bin/env python3
"""
disulfide_trajectory.py -- monitora i ponti disolfuro lungo una traiettoria.

I ponti sono quelli DICHIARATI nel topologico (legami SG-SG in [ bonds ], #include
risolti), non quelli dedotti per distanza: cosi' un ponte che si allunga o si rompe
resta in tabella invece di scomparire dal conteggio.

Le distanze usano l'immagine minima, quindi funzionano anche su frame in cui il
soluto e' avvolto sul bordo della cella (non serve passare da trjconv -pbc).

Legge solo i 2N atomi SG, non tutta la traiettoria, quindi gira in pochi secondi
anche su file da decine di GB.

Uso:
  python disulfide_trajectory.py equil_s5.xtc -s equil_s5.gro -p topol.top
  python disulfide_trajectory.py run1.xtc -s min_c3.gro -p topol.top --csv ss.csv
  python disulfide_trajectory.py *.xtc -s min_c3.gro -p topol.top --stride 5
"""
import argparse, csv, glob, os, re, sys
import numpy as np


def parse_top(path, _seen=None):
    """-> {moleculetype: (n_atomi, [(i,j) 0-based], [nomi atomo])}, [(nome, copie)]"""
    _seen = _seen or set()
    rp = os.path.realpath(path)
    if rp in _seen:
        return {}, []
    _seen.add(rp)
    mols, order, cur, name = {}, [], None, None
    base, in_mol = os.path.dirname(rp), False
    for raw in open(rp, encoding="utf-8", errors="ignore"):
        s = raw.split(";")[0].strip()
        if not s:
            continue
        m = re.match(r'#include\s+"([^"]+)"', s)
        if m:
            cand = os.path.join(base, m.group(1))
            if os.path.exists(cand):
                a, b = parse_top(cand, _seen)
                mols.update(a); order += b
            continue
        if s.startswith("#"):
            continue
        if s.startswith("[") and s.endswith("]"):
            cur = s.strip("[] ").strip().lower()
            in_mol = (cur == "molecules")
            if cur == "moleculetype":
                name = None
            continue
        p = s.split()
        if cur == "moleculetype" and name is None:
            name = p[0]; mols[name] = [0, [], []]
        elif cur == "atoms" and name and len(p) >= 5 and p[0].isdigit():
            mols[name][0] += 1; mols[name][2].append((p[3], int(p[2]), p[4]))
        elif cur == "bonds" and name and len(p) >= 2 and p[0].isdigit() and p[1].isdigit():
            mols[name][1].append((int(p[0]) - 1, int(p[1]) - 1))
        elif in_mol and len(p) >= 2:
            order.append((p[0], int(p[1])))
    return mols, order


def ss_pairs(top, molname=None):
    mols, order = parse_top(top)
    if not order:
        sys.exit(f"{top}: [ molecules ] non trovata")
    sol = molname or order[0][0]
    nat, bonds, atoms = mols[sol]
    sg = {i for i, (rn, ri, an) in enumerate(atoms) if an == "SG"}
    pairs = sorted({tuple(sorted(b)) for b in bonds if b[0] in sg and b[1] in sg})
    return pairs, atoms, sol, nat


def labels(pairs, atoms, ref_pdb, offsets):
    """etichette biologiche dal PDB di riferimento, se disponibile"""
    resid = {i: atoms[i][1] for p in pairs for i in p}
    lab = {}
    if ref_pdb and os.path.exists(ref_pdb):
        per = {}
        for l in open(ref_pdb, errors="ignore"):
            if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() not in ("NME", "ACE"):
                per.setdefault(l[21], set()).add(int(l[22:26]))
        for ch, rr in per.items():
            if ch not in offsets:
                continue
            rr = sorted(rr); run = [rr[0]]
            for r in rr[1:]:
                if r == run[-1] + 1: run.append(r)
                else: break
            for k, r in enumerate(run):
                lab[r] = f"{ch}{offsets[ch] + k}"
    def nm(i):
        return lab.get(resid[i], str(resid[i]))
    return [f"{nm(i)}-{nm(j)}" for i, j in pairs]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("traj", nargs="+")
    ap.add_argument("-s", "--struct", required=True, help="gro/pdb per la topologia mdtraj")
    ap.add_argument("-p", "--top", default="topol.top")
    ap.add_argument("--ref-pdb", default="../default_capped.pdb")
    ap.add_argument("--offsets", default="A:15,B:15,C:15,D:1,E:1")
    ap.add_argument("--break", dest="brk", type=float, default=2.5,
                    help="oltre questa distanza (A) il ponte e' considerato rotto")
    ap.add_argument("--warn", type=float, default=2.3, help="soglia di allerta (A)")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--chunk", type=int, default=200)
    ap.add_argument("--csv")
    a = ap.parse_args()
    import mdtraj as md

    pairs, atoms, sol, nat = ss_pairs(a.top)
    offs = {k: int(v) for k, v in (t.split(":") for t in a.offsets.split(","))}
    labs = labels(pairs, atoms, a.ref_pdb, offs)
    idx = sorted({i for p in pairs for i in p})
    remap = {g: k for k, g in enumerate(idx)}
    local = np.array([[remap[i], remap[j]] for i, j in pairs])
    print(f"{a.top}: molecola '{sol}' ({nat} atomi), {len(pairs)} ponti dichiarati, "
          f"{len(idx)} atomi SG da leggere")

    files = [f for pat in a.traj for f in (sorted(glob.glob(pat)) or [pat]) if os.path.exists(f)]
    if not files:
        sys.exit("nessuna traiettoria trovata")

    for traj in files:
        D, T = [], []
        for blk in md.iterload(traj, top=a.struct, chunk=a.chunk,
                               atom_indices=idx, stride=a.stride):
            D.append(md.compute_distances(blk, local, periodic=True) * 10.0)  # nm -> A
            T.append(blk.time)
        if not D:
            print(f"\n{traj}: nessun frame letto"); continue
        D = np.concatenate(D); T = np.concatenate(T)
        n = len(T)
        print(f"\n########## {os.path.basename(traj)}: {n} frame, "
              f"t = {T[0]/1000:.1f} - {T[-1]/1000:.1f} ns")
        mx = D.max(1)
        nb = (D > a.brk).any(1).sum(); nw = (D > a.warn).any(1).sum()
        print(f"  distanza massima su tutti i ponti e tutti i frame: {D.max():.2f} A "
              f"(ponte {labs[int(np.unravel_index(D.argmax(), D.shape)[1])]}, "
              f"t = {T[np.unravel_index(D.argmax(), D.shape)[0]]/1000:.2f} ns)")
        print(f"  frame con almeno un ponte > {a.warn} A: {nw}   > {a.brk} A (rotto): {nb}")
        if nb:
            print(f"\n  FRAME CON PONTI ROTTI (> {a.brk} A):")
            for f in np.where((D > a.brk).any(1))[0]:
                who = [f"{labs[k]} {D[f,k]:.2f}" for k in np.where(D[f] > a.brk)[0]]
                print(f"    frame {f:5d}  t = {T[f]/1000:8.3f} ns   " + ";  ".join(who))
        else:
            print(f"  nessun ponte rotto: tutti i {len(pairs)} disolfuri restano "
                  f"integri per tutta la traiettoria")
        sus = np.argsort(-D.max(0))[:6]
        print(f"\n  {'ponte':<14}{'media':>8}{'dev.st':>8}{'min':>7}{'max':>7}"
              f"{'frame>'+str(a.warn):>12}")
        for k in sus:
            print(f"  {labs[k]:<14}{D[:,k].mean():8.2f}{D[:,k].std():8.3f}"
                  f"{D[:,k].min():7.2f}{D[:,k].max():7.2f}{(D[:,k]>a.warn).sum():12d}")
        print(f"  (i 6 ponti con la massima escursione; media su tutti: "
              f"{D.mean():.3f} +/- {D.std():.3f} A)")
        if a.csv:
            out = a.csv if len(files) == 1 else \
                os.path.splitext(a.csv)[0] + "_" + os.path.basename(traj).split(".")[0] + ".csv"
            with open(out, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["frame", "time_ps"] + labs)
                for f in range(n):
                    w.writerow([f, f"{T[f]:.1f}"] + [f"{v:.3f}" for v in D[f]])
            print(f"  scritto {out}")


if __name__ == "__main__":
    main()
