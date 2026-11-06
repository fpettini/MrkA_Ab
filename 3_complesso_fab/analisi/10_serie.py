"""Serie per fotogramma di ogni replica (complesso e libero), dalle traiettorie della sola proteina di 09_estrai.sh.

Uso (dalla cartella complesso_fab):
    python3 analisi/10_serie.py [--trj analisi/traiettorie] [--cache analisi/cache] [--pdb mrka-fab_completo.pdb] [--rifai]

Per ogni <sistema>_run<r>_0-<FINO>ns.xtc (se ce ne sono piu' d'uno per replica, il FINO piu' alto):
 1. legge proteina_<sistema>.gro (mai il tpr), riconosce le catene (A; nel complesso anche L e H) e le confronta
    con mrka-fab_completo.pdb; controllo all'avvio A87 CYX, A94 LYS (complesso anche L214 CCYX, H220 CYX);
 2. ricompone le catene spezzate dalle condizioni periodiche (ogni N all'immagine minima del C precedente, ogni
    altro atomo all'immagine minima dell'N del suo residuo) e, nel complesso, mette L accanto ad A e H accanto a L
    (centri dei Ca all'immagine minima); controlli: legame C-N massimo, ponti S-S del cristallo, e ogni 25
    fotogrammi che nessuna immagine del Fab sia piu' vicina ad A di quella scelta;
 3. per fotogramma: RMSD dei 43 Ca del core (fit sui 43 Ca, riferimento = cristallo completato); RMSD del backbone
    dell'epitopo dopo lo stesso fit; Ca della catena A dopo lo stesso fit (per l'RMSF); distanza minima dalle
    immagini periodiche (tutti gli atomi della proteina, come gmx mindist -pi);
    solo complesso: coppie di residui A-Fab con atomi pesanti entro 5.0 A e la loro distanza minima; distanza
    minima A-Fab fra atomi pesanti (controllo del riassemblaggio, nan se > 5 A); distanza fra il centro dei Ca
    dell'Fv (VL 1-107 + VH 1-117, ricostruiti esclusi) e il centro dei 43 Ca del core di A (revisore, 5 nov);
    RMSD dei Ca del Fab e dell'Fv (fit propri, residui ricostruiti esclusi); omega delle proline cis del cristallo
    (L8, L95, L141, H151, H153); RMSD del backbone dell'epitopo dopo fit sull'epitopo stesso (misura secondaria
    dello scopo 2, aggiunta il 6 nov 2026: comune_fab.rmsd_epitopo_proprio). Nel libero la secondaria NON si
    calcola qui: la calcola 13 solo a FINALE per lunghezze (condizione del revisore, 6 nov).
Scrive CACHE/serie_<sistema>_run<r>.npz con la provenienza (sha256 di script, gro, xtc e pdb). Una serie gia'
presente si rifa' se e' cambiata la traiettoria, questo script o comune_fab.py.
Tempo: t_ns = tempo dall'inizio della produzione (primo fotogramma dell'xtc), s5 escluso.
"""
import os, sys, re, json, time, argparse, warnings
import numpy as np
warnings.filterwarnings("ignore")
import MDAnalysis as mda
from scipy.spatial import cKDTree
import comune_fab as C

ap = argparse.ArgumentParser()
ap.add_argument("--trj", default="analisi/traiettorie")
ap.add_argument("--cache", default="analisi/cache")
ap.add_argument("--pdb", default="mrka-fab_completo.pdb")
ap.add_argument("--rifai", action="store_true", help="rifa' anche le serie gia' presenti")
ap.add_argument("--blocco", type=int, default=200, help="fotogrammi letti insieme (memoria)")
arg = ap.parse_args()
os.makedirs(arg.cache, exist_ok=True)
BASI = {"TRJ": arg.trj, "CACHE": arg.cache, "RADICE": "."}
CUT_MAX = max(C.CONTATTI)
PEP_MAX, SS_MAX = 2.0, 2.5


def trova():
    best = {}
    for f in os.listdir(arg.trj):
        m = re.fullmatch(r"(complesso|libero)_(run\d+)_0-(\d+(?:\.\d+)?)ns\.xtc", f)
        if m:
            k = (m.group(1), m.group(2))
            if k not in best or float(m.group(3)) > best[k][0]:
                best[k] = (float(m.group(3)), os.path.join(arg.trj, f))
    return {k: v[1] for k, v in sorted(best.items())}


def elabora(s, r, xtc):
    gro = os.path.join(arg.trj, f"proteina_{s}.gro")
    out = os.path.join(arg.cache, f"serie_{s}_{r}.npz")
    sha_xtc = C.sha256(xtc)
    if os.path.exists(out) and not arg.rifai:
        vecchio = json.loads(str(np.load(out, allow_pickle=True)["provenienza"]))
        if (sha_xtc in vecchio["input"].values() and vecchio["sha256_script"] == C.sha256(__file__)
                and vecchio["sha256_comune_fab"] == C.sha256(C.__file__)):
            print(f"   {s} {r}: serie gia' fatta da questa traiettoria con questi script, salto (--rifai per rifarla)")
            return
    top = C.Topologia(gro, s, arg.pdb)
    xyz, _ = C.leggi_pdb(arg.pdb)
    u = mda.Universe(gro, xtc)
    if u.trajectory.n_atoms != top.n:
        raise SystemExit(f"{s} {r}: atomi nell'xtc {u.trajectory.n_atoms}, nel gro {top.n}")
    rif = lambda c, k, n: xyz[c][(k, n)]
    # gruppi di atomi
    resA = [k for k, _, _ in top.res["A"]]
    i_core = np.array([top.atomo("A", k, "CA") for k in C.CORE43]); r_core = np.array([rif("A", k, "CA") for k in C.CORE43])
    i_epi, r_epi = C.atomi_epitopo(top, arg.pdb)          # stessi atomi per primaria e secondaria dello scopo 2
    i_caA = np.array([top.atomo("A", k, "CA") for k in resA])
    pep = [(top.atomo(c, a[0], "C"), top.atomo(c, b[0], "N")) for c in C.CATENE[s] for a, b in zip(top.res[c][:-1], top.res[c][1:])]
    pep = np.array(pep)
    # ponti S-S del cristallo completato (SG-SG < 2.5 A), cercati anche fra catene diverse
    sg = [(c, k) for c in C.CATENE[s] for (k, n) in xyz[c] if n == "SG"]
    ss = [(a, b) for i, a in enumerate(sg) for b in sg[i + 1:] if np.linalg.norm(rif(*a, "SG") - rif(*b, "SG")) < SS_MAX]
    i_ss = np.array([(top.atomo(*a, "SG"), top.atomo(*b, "SG")) for a, b in ss])
    piano, gruppi = C.piano_riassemblaggio(top)
    cplx = s == "complesso"
    if cplx:
        resL = [k for k, _, _ in top.res["L"]]; resH = [k for k, _, _ in top.res["H"]]
        fab = [("L", k) for k in resL if k not in C.RICOSTRUITI["L"]] + [("H", k) for k in resH if k not in C.RICOSTRUITI["H"]]
        i_fab = np.array([top.atomo(c, k, "CA") for c, k in fab]); r_fab = np.array([rif(c, k, "CA") for c, k in fab])
        fv = [(c, k) for c, k in fab if C.FV[c][0] <= k <= C.FV[c][1]]
        i_fv = np.array([top.atomo(c, k, "CA") for c, k in fv]); r_fv = np.array([rif(c, k, "CA") for c, k in fv])
        d0 = np.linalg.norm(r_core.mean(0) - r_fv.mean(0))          # centri: 43 Ca del core e Ca dell'Fv
        omq = np.array([(top.atomo(c, k - 1, "CA"), top.atomo(c, k - 1, "C"), top.atomo(c, k, "N"), top.atomo(c, k, "CA")) for c, k in C.PRO_CIS])
        for c, k in C.PRO_CIS:
            nome = [rn for kk, rn, _ in top.res[c] if kk == k][0]
            if C.uno(nome) != "P":
                raise SystemExit(f"{c}{k} non e' una prolina ({nome})")
        hA = top.indici("A", solo_pesanti=True); hF = np.concatenate([top.indici("L", solo_pesanti=True), top.indici("H", solo_pesanti=True)])
        rA = top.resid[hA]
        cF = (top.catena[hF] == "H").astype(np.int8)        # 0 = L, 1 = H
        kF = top.resid[hF]
    # lettura a blocchi
    nf = u.trajectory.n_frames
    T = np.empty(nf); RC = np.empty(nf); RE = np.empty(nf); PI = np.empty(nf); PEP = np.empty(nf); SSd = np.empty(nf)
    CA = np.empty((nf, len(i_caA), 3), np.float32)
    if cplx:
        RF = np.empty(nf); RV = np.empty(nf); DC = np.empty(nf); RP = np.empty(nf); DM = np.full(nf, np.nan); OM = np.empty((nf, len(C.PRO_CIS)))
        cf, ca_, cc, ck, cd = [], [], [], [], []
        ass_ctrl, ass_err = 0, 0
    vicini = np.array([v for v in C._VICINI if np.any(v != 0)])
    t0 = time.time()
    for a in range(0, nf, arg.blocco):
        b = min(nf, a + arg.blocco)
        X = np.empty((b - a, top.n, 3)); M = np.empty((b - a, 3, 3))
        for j, ts in enumerate(u.trajectory[a:b]):
            X[j] = ts.positions; M[j] = ts.triclinic_dimensions; T[a + j] = ts.time
        Y = C.riassembla(X, M, piano, gruppi)
        PEP[a:b] = np.linalg.norm(Y[:, pep[:, 0]] - Y[:, pep[:, 1]], axis=-1).max(1)
        SSd[a:b] = np.linalg.norm(Y[:, i_ss[:, 0]] - Y[:, i_ss[:, 1]], axis=-1).max(1) if len(i_ss) else 0.0
        R, cP, cQ = C.kabsch(Y[:, i_core], r_core)
        RC[a:b] = C.rmsd(C.applica(Y[:, i_core], R, cP, cQ), r_core)
        RE[a:b] = C.rmsd(C.applica(Y[:, i_epi], R, cP, cQ), r_epi)
        CA[a:b] = C.applica(Y[:, i_caA], R, cP, cQ).astype(np.float32)
        if cplx:
            RP[a:b] = C.rmsd_epitopo_proprio(Y, i_epi, r_epi)     # scopo 2, secondaria: solo complesso
            R2, c2, q2 = C.kabsch(Y[:, i_fab], r_fab); RF[a:b] = C.rmsd(C.applica(Y[:, i_fab], R2, c2, q2), r_fab)
            R3, c3, q3 = C.kabsch(Y[:, i_fv], r_fv); RV[a:b] = C.rmsd(C.applica(Y[:, i_fv], R3, c3, q3), r_fv)
            DC[a:b] = np.linalg.norm(Y[:, i_core].mean(1) - Y[:, i_fv].mean(1), axis=-1)
            OM[a:b] = C.diedro(*(Y[:, omq[:, q]] for q in range(4)))
        for j in range(b - a):
            f = a + j
            PI[f] = C.mindist_immagini(Y[j], M[j])
            if cplx:
                xa, xf = Y[j, hA], Y[j, hF]
                tf = cKDTree(xf)
                sp = cKDTree(xa).sparse_distance_matrix(tf, CUT_MAX, output_type="ndarray")
                if len(sp):
                    DM[f] = sp["v"].min()
                    ii, jj, vv = sp["i"], sp["j"], sp["v"]
                    key = (rA[ii].astype(np.int64) * 2 + cF[jj]) * 1000 + kF[jj]
                    o = np.lexsort((vv, key))
                    key, vv, ii, jj = key[o], vv[o], ii[o], jj[o]
                    prima = np.r_[True, key[1:] != key[:-1]]
                    cf.append(np.full(prima.sum(), f, np.int32)); ca_.append(rA[ii[prima]].astype(np.int16))
                    cc.append(cF[jj[prima]]); ck.append(kF[jj[prima]].astype(np.int16)); cd.append(vv[prima].astype(np.float32))
                if f % 25 == 0:          # nessuna immagine del Fab piu' vicina ad A di quella scelta
                    ass_ctrl += 1
                    d_scelta = tf.query(xa, k=1)[0].min()
                    for v in vicini:
                        sh = v @ M[j]
                        if cKDTree(xf + sh).query(xa, k=1, distance_upper_bound=d_scelta)[0].min() < d_scelta - 1e-3:
                            ass_err += 1; break
        print(f"\r   {s} {r}: {b}/{nf} fotogrammi ({time.time() - t0:.0f} s)", end="", flush=True)
    print()
    # tempo e controlli
    if np.any(np.diff(T) <= 0):
        raise SystemExit(f"{s} {r}: tempi non crescenti nell'xtc")
    dt = np.diff(T)
    if dt.max() - dt.min() > 0.5:
        raise SystemExit(f"{s} {r}: passo fra fotogrammi non costante ({dt.min():.1f}-{dt.max():.1f} ps)")
    tns = np.round((T - T[0]) / 1000.0, 4)
    if PEP.max() >= PEP_MAX or (len(i_ss) and SSd.max() >= SS_MAX):
        raise SystemExit(f"{s} {r}: ricomposizione non riuscita (C-N max {PEP.max():.2f} A, S-S max {SSd.max():.2f} A)")
    if cplx and ass_err:
        raise SystemExit(f"{s} {r}: in {ass_err}/{ass_ctrl} fotogrammi controllati un'immagine del Fab e' piu' vicina ad A di quella scelta")
    P = C.prov(__file__, [gro, xtc, arg.pdb], {"sistema": s, "replica": r, "trj": arg.trj, "pdb": arg.pdb}, BASI)
    P["controllo_avvio"] = top.controllo_avvio
    dati = dict(t_ns=tns, tinit_ps=T[0], dt_ps=float(np.median(dt)), rmsd_core=RC, rmsd_epi=RE, ca_fit=CA,
                ca_resid=np.array(resA), mindist_pi=PI, pep_max=PEP, ss_max=SSd, n_ss=len(i_ss),
                provenienza=np.array(json.dumps(P)))
    msg = ""
    if cplx:
        dati.update(rmsd_epi_proprio=RP, rmsd_fab=RF, rmsd_fv=RV, dcentri=DC, dcentri_cristallo=d0, dmin_afab=DM, omega=OM,
                    pro_cis=np.array([f"{c}{k}" for c, k in C.PRO_CIS]),
                    c_f=np.concatenate(cf) if cf else np.zeros(0, np.int32), c_a=np.concatenate(ca_) if ca_ else np.zeros(0, np.int16),
                    c_c=np.concatenate(cc) if cc else np.zeros(0, np.int8), c_k=np.concatenate(ck) if ck else np.zeros(0, np.int16),
                    c_d=np.concatenate(cd) if cd else np.zeros(0, np.float32), c_cut=CUT_MAX,
                    assemblaggio=np.array([ass_ctrl, ass_err]))
        dm = "> 5" if np.all(np.isnan(DM)) else f"{np.nanmax(DM):.2f}" + (" (> 5 A in alcuni)" if np.any(np.isnan(DM)) else "")
        msg = (f", centri Fv-core {DC.min():.1f}-{DC.max():.1f} A (cristallo {d0:.1f}), assemblaggio ok in {ass_ctrl} controlli, "
               f"distanza minima A-Fab al massimo {dm} A")
    np.savez(out, **dati)
    pi = np.nanmin(PI) if np.any(~np.isnan(PI)) else float("nan")
    print(f"   {s} {r}: {nf} fotogrammi, tinit {T[0]:.0f} ps, 0-{tns[-1]:.1f} ns ogni {np.median(dt):.0f} ps; "
          f"C-N max {PEP.max():.2f} A, S-S max {SSd.max():.2f} A ({len(i_ss)} ponti); RMSD core {RC.mean():.2f} A; "
          f"immagini min {pi:.1f} A{msg}")


if __name__ == "__main__":
    lista = trova()
    if not lista:
        raise SystemExit(f"nessuna traiettoria in {arg.trj} (lanciare prima 09_estrai.sh)")
    for (s, r), xtc in lista.items():
        elabora(s, r, xtc)
