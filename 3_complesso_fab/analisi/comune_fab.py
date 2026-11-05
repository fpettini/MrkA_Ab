"""Definizioni e funzioni comuni alle analisi del complesso MrkA-Fab e della MrkA libera.

Le definizioni (residui, soglie, scarti, conteggi) vengono dal documento "MrkA-Fab: criteri di analisi"
(https://claude.ai/artifact/DPDeZEeE17qTTZxD3ukXdX), chiuso il 2 nov 2026 prima di aprire le traiettorie.
Le scelte operative che il documento non fissava sono marcate "NON NEL DOC" (rev. 5 nov): il revisore le ha decise
il 5 nov 2026, prima di aprire le traiettorie, e sono nel registro del documento.

Numerazione: wt della catena A (costrutto 36-217), Fab L 1-214 e H 1-223, come in mrka-fab_completo.pdb
(catene A, B = L, C = H). Nei gro i numeri si ripetono fra le catene e non c'e' l'identificativo di catena:
le catene si riconoscono dal salto di numerazione e dalla sequenza, e si seleziona sempre per catena.
Unita': distanze in A, tempi in ns (il tempo 0 e' l'inizio della produzione, s5 escluso).
"""
import os, json, math, hashlib, platform, datetime, itertools
import numpy as np
from importlib import metadata

VERSIONE = "analisi_fab 2026-11-05b"
PACCHETTI = ["numpy", "scipy", "MDAnalysis", "matplotlib"]

# ---------------- definizioni fissate nel documento dei criteri ----------------
CORE43 = [44, 45, 46, 47, 78, 79, 80, 81, 82, 83, 84, 85, 88, 94, 102, 103, 104, 105, 107, 108, 122, 134,
          135, 136, 137, 138, 148, 162, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 195, 196, 197, 198, 199]
EPITOPO = [42, 86, 88, 90, 91, 92, 93, 94, 96, 167, 168]
BACKBONE = ("N", "CA", "C", "O")
SEGMENTI = {"87-101": list(range(87, 102)), "138-147": list(range(138, 148))}
PRO_CIS = [("L", 8), ("L", 95), ("L", 141), ("H", 151), ("H", 153)]       # proline cis del cristallo
SCARTI = (50.0, 100.0)                     # ns: principale e sensibilita'
CONTATTI = (3.5, 4.0, 5.0)                 # A: principale 4.0, sensibilita' 3.5 e 5.0
CONTATTO_PRINCIPALE = 4.0
DISEGNO = {"complesso": (5, 500.0), "libero": (10, 1000.0)}   # repliche, durata massima (ns)
MIN_VALIDE = {"complesso": 4, "libero": 8}  # sotto: "non concludibile"
RMSD_CORE_MAX, RMSD_CORE_FRAZ = 2.0, 0.5   # replica non valida: RMSD 43 Ca > 2.0 A in > 50% dei fotogrammi
PI_MIN, PI_FRAZ = 12.0, 0.05               # mindist -pi < 1.2 nm: segnalata se succede, non valida in > 5%
CONTATTI_PERSI_FRAZ = 0.5                  # complesso: non valida se perde > 50% dei contatti del cristallo
OCC_PERSO = 0.5                            # NON NEL DOC (rev. 5 nov): perso = occupanza < 0.5 nella finestra (come lo scopo 1)
CENTRI_MAX = 5.0                           # complesso: non valida se la distanza dei centri cresce > 5 A
CENTRI_FRAZ = 0.5                          # NON NEL DOC (rev. 5 nov): ... in > 50% dei fotogrammi della finestra
# NON NEL DOC (rev. 5 nov): centri = Ca dell'Fv (VL 1-107 + VH 1-117, ricostruiti esclusi) e i 43 Ca del core di A
AFAB_MAX = 4.0                             # NON NEL DOC (rev. 5 nov): riassemblaggio, distanza minima A-Fab < 4 A
PI_TOLL = 0.005                            # NON NEL DOC (rev. 5 nov): accordo con gmx mindist -pi (A)
SOGLIA_EPI_FLAG = 2.0                      # scopo 2: si segnala una soglia sopra 2 A
ISTERESI_NS = 1.0                          # scopo 2: uscita/rientro = almeno 1 ns di fila
CICLI_FRAZ = 0.5                           # scopo 2: cicli in almeno 5/10 repliche; con meno valide ceil(0.5 x valide) (rev. 5 nov)
# NON NEL DOC (rev. 5 nov): residui ricostruiti, esclusi dall'RMSD del Fab e dell'Fv e dai centri; Fv = VL 1-107 + VH 1-117
RICOSTRUITI = {"A": [116, 117], "L": [214], "H": [1, 9, 10, 218, 219, 220, 221, 222, 223]}
FV = {"L": (1, 107), "H": (1, 117)}
BLOCCO_SCOPO2_NS = 50.0                    # NON NEL DOC (rev. 5 nov): blocchi del decorso dello scopo 2 (descrittivo)
PI_CUT = 25.0                              # oltre 25 A la distanza dalle immagini si registra come "> 25"

UNO = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V', 'CYX': 'C', 'HIE': 'H', 'HID': 'H', 'HIP': 'H', 'LYN': 'K', 'ASH': 'D', 'GLH': 'E',
       'CYM': 'C'}
CATENE = {"complesso": ["A", "L", "H"], "libero": ["A"]}
PDB_CATENA = {"A": "A", "B": "L", "C": "H"}           # mrka-fab_completo.pdb -> nomi delle analisi


def uno(n):
    if n in UNO:
        return UNO[n]
    if len(n) == 4 and n[0] in "NC" and n[1:] in UNO:
        return UNO[n[1:]]
    return None


def pesante(nome):
    return not nome.startswith("H")


# ---------------- provenienza (stessa regola di analisi_rest2/provenienza.py) ----------------
def sha256(path, blocco=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(blocco), b""):
            h.update(b)
    return h.hexdigest()


def etichetta(path, basi):
    a = os.path.abspath(path)
    for nome, d in (basi or {}).items():
        d = os.path.abspath(d)
        if a == d:
            return f"<{nome}>"
        if a.startswith(d + os.sep):
            return f"<{nome}>/" + os.path.relpath(a, d).replace(os.sep, "/")
    return os.path.basename(a)


def prov(script, inputs, argomenti, basi=None):
    pk = {}
    for p in PACCHETTI:
        try:
            pk[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            pk[p] = None
    return {"versione": VERSIONE, "script": os.path.basename(script), "sha256_script": sha256(script),
            "sha256_comune_fab": sha256(os.path.abspath(__file__)),
            "data": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
            "python": platform.python_version(), "pacchetti": pk,
            "argomenti": {k: (etichetta(v, basi) if isinstance(v, str) and os.path.exists(v) else v)
                          for k, v in argomenti.items()},
            "input": {etichetta(i, basi): sha256(i) for i in inputs if os.path.exists(i)}}


def testo_prov(p):
    r = [f"# {p['script']} ({p['versione']}), sha256 {p['sha256_script'][:16]}, comune_fab.py {p['sha256_comune_fab'][:16]}",
         f"# eseguito {p['data']}, Python {p['python']}, " + ", ".join(f"{k} {v}" for k, v in p["pacchetti"].items()),
         f"# argomenti: {p['argomenti']}"]
    r += [f"# input {k}: sha256 {v[:16]}" for k, v in p["input"].items()]
    return r


# ---------------- struttura di riferimento (cristallo completato) ----------------
def leggi_pdb(path):
    """mrka-fab_completo.pdb -> {catena: {(resid, nome): xyz}}, {catena: {resid: resname}}; solo ATOM, niente H."""
    xyz, rn = {}, {}
    for l in open(path):
        if not l.startswith("ATOM") or l[16] not in " A":
            continue
        c = PDB_CATENA.get(l[21])
        nome = l[12:16].strip()
        el = l[76:78].strip()
        if c is None or el == "H" or nome.startswith("H"):
            continue
        k = int(l[22:26])
        xyz.setdefault(c, {})[(k, nome)] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        rn.setdefault(c, {})[k] = l[17:20]
    return xyz, rn


def contatti_struttura(xa, ra, xf, rf, cut):
    """Coppie di residui (resA, catena Fab, resF) con almeno una coppia di atomi pesanti entro cut (<=)."""
    from scipy.spatial import cKDTree
    t = cKDTree(xf)
    out = set()
    for i, js in enumerate(t.query_ball_point(xa, cut)):
        for j in js:
            out.add((int(ra[i]), rf[j][0], int(rf[j][1])))
    return out


def contatti_cristallo(pdb, cut):
    xyz, _ = leggi_pdb(pdb)
    ka = [k for k in xyz["A"] if pesante(k[1])]
    kf = [(c, k) for c in ("L", "H") for k in xyz[c] if pesante(k[1])]
    xa = np.array([xyz["A"][k] for k in ka]); ra = [k[0] for k in ka]
    xf = np.array([xyz[c][k] for c, k in kf]); rf = [(c, k[0]) for c, k in kf]
    return contatti_struttura(xa, ra, xf, rf, cut)


# ---------------- topologia dal gro della sola proteina ----------------
class Topologia:
    """Atomi della proteina dal gro (gruppo Protein_Other), catene riconosciute e controllate contro il pdb."""

    def __init__(self, gro, sistema, pdb):
        import MDAnalysis as mda
        if sistema not in CATENE:
            raise SystemExit(f"sistema sconosciuto: {sistema}")
        self.sistema = sistema
        u = mda.Universe(gro)                      # mai il tpr: MDAnalysis rinumera i residui
        self.n = u.atoms.n_atoms
        self.nomi = np.array(u.atoms.names)
        self.resnomi = np.array(u.atoms.resnames)
        self.resid = np.array(u.atoms.resids)
        # catene: nuovo inizio quando il numero non e' il precedente + 1 (come 08_verifica_numerazione.py)
        res = []
        for r in u.residues:
            res.append((int(r.resid), r.resname, r.atoms.indices))
        for k, r in enumerate(res):
            if uno(r[1]) is None:
                raise SystemExit(f"{gro}: residuo non proteico {r[1]}{r[0]} nel gruppo della proteina")
        catene, cur = [], []
        for r in res:
            if cur and r[0] != cur[-1][0] + 1:
                catene.append(cur); cur = []
            cur.append(r)
        catene.append(cur)
        attese = CATENE[sistema]
        if len(catene) != len(attese):
            raise SystemExit(f"{gro}: {len(catene)} catene, attese {len(attese)} ({sistema})")
        _, rnp = leggi_pdb(pdb)
        self.catena = np.empty(self.n, dtype="<U1")
        self.res = {}                               # catena -> lista (resid, resname, indici)
        for nome, c in zip(attese, catene):
            seq = "".join(uno(r[1]) for r in c)
            ref = rnp[nome]
            seqp = "".join(UNO[ref[k]] for k in sorted(ref))
            if seq != seqp or c[0][0] != min(ref) or c[-1][0] != max(ref):
                raise SystemExit(f"{gro}: la catena {nome} ({c[0][0]}-{c[-1][0]}) non coincide con mrka-fab_completo.pdb")
            self.res[nome] = c
            for r in c:
                self.catena[r[2]] = nome
        # controllo all'avvio deciso dal revisore il 2 nov 2026
        attesi = [("A", 87, "CYX"), ("A", 94, "LYS")]
        if sistema == "complesso":
            attesi += [("L", 214, "CCYX"), ("H", 220, "CYX")]
        for c, k, n in attesi:
            trovato = [r[1] for r in self.res[c] if r[0] == k]
            if trovato != [n]:
                raise SystemExit(f"{gro}: controllo all'avvio fallito, {c}{k} = {trovato}, atteso {n}")
        self.controllo_avvio = ", ".join(f"{c}{k} {n}" for c, k, n in attesi)

    def atomo(self, catena, k, nome):
        for r in self.res[catena]:
            if r[0] == k:
                w = [i for i in r[2] if self.nomi[i] == nome]
                if len(w) == 1:
                    return int(w[0])
        raise SystemExit(f"atomo {catena}{k} {nome} non trovato")

    def indici(self, catena, residui=None, nomi=None, solo_pesanti=False):
        out = []
        for r in self.res[catena]:
            if residui is not None and r[0] not in residui:
                continue
            for i in r[2]:
                if nomi is not None and self.nomi[i] not in nomi:
                    continue
                if solo_pesanti and not pesante(self.nomi[i]):
                    continue
                out.append(int(i))
        return np.array(out, dtype=int)

    def piano_ricomposizione(self):
        """Per ogni residuo: (C del residuo precedente o -1, N, altri atomi). Ogni N si mette all'immagine minima
        rispetto al C precedente, ogni altro atomo rispetto all'N del suo residuo (vettori < 10 A)."""
        piano = []
        for c in CATENE[self.sistema]:
            cprec = -1
            for k, rn, idx in self.res[c]:
                n = [i for i in idx if self.nomi[i] == "N"]
                cc = [i for i in idx if self.nomi[i] == "C"]
                if len(n) != 1 or len(cc) != 1:
                    raise SystemExit(f"residuo {c}{k}{rn} senza N o C unici")
                altri = np.array([i for i in idx if i != n[0]], dtype=int)
                piano.append((cprec, int(n[0]), altri))
                cprec = int(cc[0])
        return piano


# ---------------- condizioni periodiche ----------------
def immagine_minima(v, M, Minv):
    """v (F,k,3), M (F,3,3) con i vettori di cella per righe: toglie i vettori di cella interi (vettori corti)."""
    f = np.einsum("fki,fij->fkj", v, Minv)
    return v - np.einsum("fki,fij->fkj", np.round(f), M)


_VICINI = np.array(list(itertools.product((-1, 0, 1), repeat=3)), dtype=float)


def immagine_minima_vera(v, M, Minv):
    """v (F,3): immagine di lunghezza minima, cercata fra le 27 celle vicine dopo l'arrotondamento."""
    f = np.einsum("fi,fij->fj", v, Minv)
    base = v - np.einsum("fi,fij->fj", np.round(f), M)
    cand = base[:, None, :] + np.einsum("ni,fij->fnj", _VICINI, M)
    k = np.argmin((cand ** 2).sum(-1), axis=1)
    return cand[np.arange(len(v)), k]


def ricomponi(X, M, Minv, piano):
    Y = X.copy()
    for cprec, n, altri in piano:
        if cprec >= 0:
            Y[:, n] = Y[:, cprec] + immagine_minima(X[:, n][:, None] - Y[:, cprec][:, None], M, Minv)[:, 0]
        if len(altri):
            Y[:, altri] = Y[:, n][:, None] + immagine_minima(X[:, altri] - Y[:, n][:, None], M, Minv)
    return Y


def assembla(Y, M, Minv, gruppi):
    """gruppi: [(atomi, ca, atomi_ancora, ca_ancora)]: sposta ogni gruppo di un vettore di cella intero perche' il
    centro dei suoi Ca sia all'immagine minima rispetto al centro dei Ca dell'ancora (gia' sistemata)."""
    for atomi, ca, _, ca_anc in gruppi:
        d = Y[:, ca].mean(1) - Y[:, ca_anc].mean(1)
        s = immagine_minima_vera(d, M, Minv) - d
        Y[:, atomi] += s[:, None, :]
    return Y


def piano_riassemblaggio(top):
    """Piano per ricomponi() e gruppi per assembla(): nel complesso L accanto ad A e H accanto a L (10 e 16)."""
    gruppi = []
    if top.sistema == "complesso":
        gruppi = [(top.indici("L"), top.indici("L", nomi=("CA",)), None, top.indici("A", nomi=("CA",))),
                  (top.indici("H"), top.indici("H", nomi=("CA",)), None, top.indici("L", nomi=("CA",)))]
    return top.piano_ricomposizione(), gruppi


def riassembla(X, M, piano, gruppi):
    """X (F,N,3), M (F,3,3): catene ricomposte e, nel complesso, riassemblate."""
    Minv = np.linalg.inv(M)
    Y = ricomponi(X, M, Minv, piano)
    return assembla(Y, M, Minv, gruppi) if gruppi else Y


def mindist_immagini(X, M, cut=PI_CUT):
    """Distanza minima fra la proteina (X, N x 3, ricomposta) e le sue immagini periodiche, come gmx mindist -pi
    (tutti gli atomi del gruppo). Restituisce nan se nessuna coppia e' entro cut."""
    from scipy.spatial import cKDTree
    c = X.mean(0)
    R = np.sqrt(((X - c) ** 2).sum(1)).max()
    best = np.inf
    for n in itertools.product(range(-2, 3), repeat=3):
        if n == (0, 0, 0) or next(x for x in n if x != 0) < 0:      # n e -n danno le stesse distanze
            continue
        T = np.array(n, dtype=float) @ M
        L = np.linalg.norm(T)
        lim = min(cut, best)
        if L > 2 * R + lim:
            continue
        u = T / L
        p = X @ u
        I = np.where(p >= p.min() + L - lim)[0]
        J = np.where(p <= p.max() - L + lim)[0]
        if len(I) == 0 or len(J) == 0:
            continue
        d, _ = cKDTree(X[I]).query(X[J] + T, k=1, distance_upper_bound=lim)
        m = d.min()
        if m < best:
            best = m
    return best if best <= cut else np.nan


# ---------------- sovrapposizione ----------------
def kabsch(P, Q):
    """P (F,k,3) mobili, Q (k,3) riferimento: rotazioni R (F,3,3) e centri, tali che (P - cP) @ R ~ Q - cQ."""
    cP = P.mean(1)
    cQ = Q.mean(0)
    H = np.einsum("fki,kj->fij", P - cP[:, None], Q - cQ)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("fij,fjk->fik", U, Vt)))
    D = np.ones((len(P), 3)); D[:, 2] = d
    R = np.einsum("fij,fj,fjk->fik", U, D, Vt)
    return R, cP, cQ


def applica(X, R, cP, cQ):
    return np.einsum("fki,fij->fkj", X - cP[:, None], R) + cQ


def rmsd(A, B):
    return np.sqrt(((A - B) ** 2).sum(-1).mean(-1))


def diedro(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


# ---------------- finestre, conteggi, stato dell'analisi ----------------
def soglia_repliche(n_valide):
    return math.ceil(0.8 * n_valide - 1e-9)


def soglia_cicli(n_valide):
    """scopo 2: repliche con almeno un ciclo uscita-rientro; 5 con 10 valide, ceil(0.5 x valide) con meno (rev. 5 nov)"""
    return math.ceil(CICLI_FRAZ * n_valide - 1e-9)


def finestra(t, scarto, fine):
    return (t >= scarto - 1e-6) & (t <= fine + 1e-6)


def leggi_cache(cache):
    """{(sistema, replica): dati npz} per tutte le serie in CACHE."""
    out = {}
    if not os.path.isdir(cache):
        raise SystemExit(f"cartella delle serie assente: {cache} (lanciare prima 10_serie.py)")
    for f in sorted(os.listdir(cache)):
        if f.startswith("serie_") and f.endswith(".npz"):
            s, r = f[6:-4].split("_")
            out[(s, r)] = dict(np.load(os.path.join(cache, f), allow_pickle=True))
            out[(s, r)]["_file"] = os.path.join(cache, f)
    if not out:
        raise SystemExit(f"nessuna serie in {cache}")
    return out


def ordina_repliche(rep):
    return sorted(rep, key=lambda r: int(r[3:]))


def durate(serie):
    """Durate confermative (criteri, sezione 4): T3 = minimo fra le repliche del complesso (max 500 ns), usata
    per gli scopi 1 e 3 e per la soglia dello scopo 2; T2 = minimo fra le repliche del libero (max 1000 ns)."""
    lung = {}
    for (s, r), d in serie.items():
        lung.setdefault(s, {})[r] = float(d["t_ns"][-1])
    T3 = min(min(lung.get("complesso", {"-": 0.0}).values()), DISEGNO["complesso"][1])
    T2 = min(min(lung.get("libero", {"-": 0.0}).values()), DISEGNO["libero"][1])
    return T3, T2, lung


def stato_analisi(lung, T3, T2, finale=False):
    """Analisi FINALE per uno scopo solo con tutte le repliche del disegno alla durata massima, oppure con --finale
    (riduzione del disegno decisa da francesco prima di guardare le traiettorie). Altrimenti INTERMEDIA."""
    nc = len(lung.get("complesso", {})); nl = len(lung.get("libero", {}))
    c_ok = nc == DISEGNO["complesso"][0] and T3 >= DISEGNO["complesso"][1] - 1e-6
    l3 = nl == DISEGNO["libero"][0] and min(lung.get("libero", {"-": 0}).values()) >= T3 - 1e-6
    l2 = nl == DISEGNO["libero"][0] and T2 >= DISEGNO["libero"][1] - 1e-6
    f = lambda ok: "FINALE" if (ok or finale) else "INTERMEDIA"
    return {"scopo1": f(c_ok), "scopo2": f(c_ok and l2), "scopo3": f(c_ok and l3)}


def leggi_controlli(out):
    f = os.path.join(out, "controlli.json")
    if not os.path.exists(f):
        raise SystemExit(f"{f} assente: lanciare prima 11_controlli.py")
    return json.load(open(f))


def valide(ctrl, sistema, finestra_, scarto):
    """Repliche valide per sistema, finestra e scarto."""
    return ctrl["riepilogo"][f"{sistema}|{finestra_}|{int(scarto)}"]


def avviso_intermedia(stato):
    if stato == "FINALE":
        return []
    return ["# ANALISI INTERMEDIA, solo descrittiva (condizioni del revisore, 5 nov 2026): numeri per replica; le regole",
            "# di conteggio (ceil 0.8, 4/5, 8/10, cicli, 'non concludibile') e gli esiti non si applicano; niente test;",
            "# nessuna frase del paper da questi numeri; nessuna replica si rilancia o si sostituisce per quello che si vede."]
