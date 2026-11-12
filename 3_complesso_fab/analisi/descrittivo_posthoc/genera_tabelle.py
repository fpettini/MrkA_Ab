"""Tabelle per il testo dell'articolo (EN) e per il report (IT), generate SOLO dai file dei risultati registrati.
Nessun ricalcolo dalle traiettorie: ogni numero e' una foglia dei json (o, per la distanza nel cristallo, del pdb).
Le tabelle S2 (continua 2), S4 (continua) e S6 vengono dai json dei descrittivi (d9, d10, d11), se la cartella e' data.
Uso: python3 genera_tabelle.py <cartella risultati> <pdb> <cartella uscita> [cartella dei json d9-d11] [cartella degli script di analisi]"""
import json, re, sys, os, hashlib
import numpy as np
sys.path.insert(0, sys.argv[5] if len(sys.argv) > 5 else os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import comune_fab as C

RIS, PDB, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
DESC = sys.argv[4] if len(sys.argv) > 4 else None
L = lambda f: json.load(open(os.path.join(RIS, f)))
sha = lambda f: hashlib.sha256(open(os.path.join(RIS, f), "rb").read()).hexdigest()[:16]
s1, s2, s3, ct, cer = L("scopo1_contatti.json"), L("scopo2_epitopo.json"), L("scopo3_rmsf.json"), L("controlli.json"), L("diagnostica_cerniera.json")
assert s1["stato"] == s2["stato"] == s3["stato"] == cer["stato"] == "FINALE"
xyz, rn = C.leggi_pdb(PDB)
tre = lambda r: r.capitalize()
CR = ["run1", "run2", "run3", "run4", "run5"]; FR = [f"run{i}" for i in range(1, 11)]
ES = {"persistente": "persistent", "perso": "lost", "senza esito": "no outcome", "nuovo": "new",
      "popolata": "populated", "intermedia": "intermediate", "non popolata": "not populated",
      "non distinguibile dalla memoria della partenza": "indistinguishable from memory of the start"}

def dmin(a, c, k):
    A = [v for (r, n), v in xyz["A"].items() if r == a and C.pesante(n)]
    B = [v for (r, n), v in xyz[c].items() if r == k and C.pesante(n)]
    return min(float(np.linalg.norm(x - y)) for x in A for y in B)

def chiave(p):
    a, c, k = p.split("|"); return int(a), c, int(k)

def nomeEN(p):
    a, c, k = chiave(p); return f"{tre(rn['A'][a])}{a}", f"{c} {tre(rn[c][k])}{k}"

def md(righe, intest):
    out = ["| " + " | ".join(intest) + " |", "|" + "|".join(["---"] * len(intest)) + "|"]
    out += ["| " + " | ".join(str(x) for x in r) + " |" for r in righe]
    return "\n".join(out)

EN, IT, F = [], [], []          # F = fatti (numeri per il testo)

# ---------------- Tabella S1: contatti del cristallo a 4.0 A
a, b = s1["risultati"]["4.0"]["50"]["coppie"], s1["risultati"]["4.0"]["100"]["coppie"]
cris = sorted([p for p in a if a[p]["cristallo"]], key=chiave)
assert len(cris) == 27
righe = []; righeIT = []
for p in cris:
    assert a[p]["esito"] == b[p]["esito"], p
    ra, rf = nomeEN(p); x, c, k = chiave(p)
    fr = a[p]["frazioni"]
    righe.append([ra, rf, f"{dmin(x, c, k):.2f}"] + [f"{v:.2f}" for v in fr] + [ES[a[p]["esito"]]])
    righeIT.append([f"{ra} – {rf}"] + [f"{v:.2f}" for v in fr] + [a[p]["esito"]])
EN.append("**Table S1. Occupancy of the 27 crystal contacts in the five complex replicas.** Fraction of frames (50–500 ns) with any two heavy atoms of the pair within 4.0 Å. Crystal distance: minimum heavy-atom distance in the completed model. Outcomes are identical with the 100-ns discard. H, heavy chain; L, light chain.\n\n"
          + md(righe, ["MrkA", "Fab", "Crystal distance (Å)", "Replica 1", "2", "3", "4", "5", "Outcome"]))
IT.append("**Contatti del cristallo a 4.0 Å: frazione dei fotogrammi per replica (50–500 ns) ed esito.** Gli esiti sono gli stessi con scarto 100 ns.\n\n"
          + md(righeIT, ["Coppia", "run1", "run2", "run3", "run4", "run5", "Esito"]))
n50 = {p: sum(v >= 0.5 for v in a[p]["frazioni"]) for p in cris}; n100 = {p: sum(v >= 0.5 for v in b[p]["frazioni"]) for p in cris}
pers = [p for p in cris if a[p]["esito"] == "persistente"]
minimo = [p for p in pers if n50[p] == 4]; scende = [p for p in pers if n50[p] == 5 and n100[p] == 4]
assert all(n100[p] >= 4 for p in pers)
coppia = lambda p: "–".join(nomeEN(p))
basso = lambda p, d: min(v for v in d[p]["frazioni"] if v >= 0.5)
NOTA_S1 = ("Persistent by the minimum count (4 of 5 replicas at or above 0.5; in brackets the lowest of the four, with the 50- and 100-ns discards): " + ", ".join(f"{coppia(p)} ({basso(p, a):.2f}, {basso(p, b):.2f})" for p in minimo)
           + ". Persistent in 5 replicas with the 50-ns discard and in 4 with the 100-ns discard: " + ", ".join(coppia(p) for p in scende) + ".")
F.append("scopo1 4.0 A: " + NOTA_S1)
EN[-1] = EN[-1].replace("H, heavy chain; L, light chain.", NOTA_S1 + " H, heavy chain; L, light chain.")
assert NOTA_S1 in EN[-1]
nc = sorted([p for p in a if not a[p]["cristallo"] and max(a[p]["frazioni"]) >= 0.5], key=chiave)
righe = []
for p in nc:
    ra, rf = nomeEN(p); x, c, k = chiave(p)
    righe.append([ra, rf, f"{dmin(x, c, k):.2f}"] + [f"{v:.2f}" for v in a[p]["frazioni"]])
EN.append("**Table S1 (continued). Pairs not in contact in the crystal at 4.0 Å with occupancy ≥ 0.5 in at least one replica.** None meets the rule for a new contact (≥ 0.5 in at least 4 of 5 replicas).\n\n"
          + md(righe, ["MrkA", "Fab", "Crystal distance (Å)", "Replica 1", "2", "3", "4", "5"]))
from collections import Counter
for cut in ("3.5", "4.0", "5.0"):
    A_, B_ = s1["risultati"][cut]["50"]["coppie"], s1["risultati"][cut]["100"]["coppie"]
    fin = Counter()
    for p in set(A_) | set(B_):
        ea = A_[p]["esito"] if p in A_ else "senza esito"; eb = B_[p]["esito"] if p in B_ else "senza esito"
        cr = (A_.get(p) or B_.get(p))["cristallo"]
        e = ea if ea == eb else "non robusto"
        fin[("cristallo" if cr else "non cristallo", e)] += 1
    F.append(f"scopo1 {cut} A: coppie nel cristallo {sum(v for (c_, e), v in fin.items() if c_ == 'cristallo')}; " + "; ".join(f"{c_} {e} {v}" for (c_, e), v in sorted(fin.items())))
sat = [p for p in cris if min(a[p]["frazioni"]) >= 0.88]
F.append("scopo1 4.0 A: coppie con frazione >= 0.88 in ogni replica: " + ", ".join("–".join(nomeEN(p)) + f" (min {min(a[p]['frazioni']):.2f})" for p in sat))
perres = {}
for p in cris:
    perres.setdefault(chiave(p)[0], []).append(a[p]["esito"])
F.append("scopo1 4.0 A per residuo persistenti/totali: " + ", ".join(f"{r} {v.count('persistente')}/{len(v)}" for r, v in sorted(perres.items())))

# ---------------- Tabella S2: scopo 2
def conteggi(P):
    out = {}
    for nome in P["soglie"]:
        fr = [P["libero"][r][nome]["frazione"] for r in FR]; cy = [P["libero"][r][nome]["cicli"] for r in FR]
        out[nome] = (sum(f >= 0.5 for f in fr), sum(f <= 0.05 for f in fr), sum(c >= 1 for c in cy))
    return out
NOMI = {"P95": "95th percentile, pooled (criterion)", "P90": "90th percentile, pooled", "P99": "99th percentile, pooled",
        "P95_replica_min": "smallest per-replica 95th percentile", "P95_replica_max": "largest per-replica 95th percentile"}
righe = []
for mis, ramo in (("Core-fitted (registered)", s2), ("Self-fitted (post hoc)", s2["secondaria"])):
    for sc in ("50", "100"):
        P = ramo["per_scarto"][sc]; cn = conteggi(P)
        assert P["k"] == 8 and P["k_cicli"] == 5
        for nome in ("P95", "P90", "P99", "P95_replica_min", "P95_replica_max"):
            g, l, c = cn[nome]
            et = NOMI[nome] if ramo is s2 else NOMI[nome].replace("(criterion)", "(same definition as the criterion)")
            righe.append([mis, sc, et, f"{P['soglie'][nome]:.2f}", g, l, c, ES[P["esiti"][nome]]])
        F.append(f"scopo2 {mis} scarto {sc}: " + "; ".join(f"{n} soglia {P['soglie'][n]:.3f} ge0.5 {cn[n][0]} le0.05 {cn[n][1]} cicli {cn[n][2]} esito {P['esiti'][n]}" for n in P["soglie"]))
        F.append(f"scopo2 {mis} scarto {sc}: P95 per replica del complesso " + ", ".join(f"{P['soglia_per_replica'][r]:.2f}" for r in CR) + "; RMSD medio complesso " + ", ".join(f"{P['complesso_media'][r]:.2f}" for r in CR))
EN.append("**Table S2. Bound envelope of the epitope and counts over the ten free replicas.** Envelope: percentile of the epitope backbone RMSD from the crystal in the complex. Counts: free replicas with a fraction of frames inside the envelope ≥ 0.5 (8 required for \"populated\"), ≤ 0.05 (8 required for \"not populated\"), and with at least one exit and re-entry (5 required). The pooled 95th percentile is the envelope of the criteria; the other four definitions are the sensitivity checks fixed with them. Free MrkA window: discard to 1000 ns. Self-fitted values pending confirmation on the workstation.\n\n"
          + md(righe, ["Measure", "Discard (ns)", "Envelope definition", "Envelope (Å)", "≥ 0.5", "≤ 0.05", "Exit and re-entry", "Outcome"]))
righe = []; righeIT = []
P1, P2 = s2["per_scarto"]["50"], s2["secondaria"]["per_scarto"]["50"]
for i, r in enumerate(FR, 1):
    x, y = P1["libero"][r], P2["libero"][r]
    righe.append([i, f"{x['rmsd_media']:.2f}", f"{x['P95']['frazione']:.2f}", x["P95"]["uscite"], x["P95"]["rientri"],
                  f"{y['rmsd_media']:.2f}", f"{y['P95']['frazione']:.2f}", y["P95"]["uscite"], y["P95"]["rientri"]])
    righeIT.append([r, f"{x['P95']['frazione']:.2f}", x["P95"]["rientri"], f"{y['P95']['frazione']:.2f}", y["P95"]["rientri"]])
EN.append("**Table S2 (continued). Free replicas, 50–1000 ns, pooled 95th-percentile envelope** (core-fitted 3.88 Å; self-fitted 1.54 Å). Mean RMSD in Å; fraction of frames inside the envelope; exits and re-entries counted with a persistence of at least 1 ns.\n\n"
          + md(righe, ["Free replica", "Core-fitted: mean RMSD", "fraction inside", "exits", "re-entries", "Self-fitted: mean RMSD", "fraction inside", "exits", "re-entries"]))
IT.append("**Scopo 2, repliche del libero (50–1000 ns, soglia P95):** frazione dei fotogrammi dentro l'involucro e rientri (cicli).\n\n"
          + md(righeIT, ["Replica", "Primaria: frazione (3.88 Å)", "cicli", "Secondaria: frazione (1.54 Å)", "cicli"]))
for nome, P in (("primaria", P1), ("secondaria", P2)):
    F.append(f"scopo2 {nome} scarto 50 P95: frazioni " + ", ".join(f"{P['libero'][r]['P95']['frazione']:.2f}" for r in FR) + "; rientri " + ", ".join(str(P['libero'][r]['P95']['rientri']) for r in FR)
             + "; uscite " + ", ".join(str(P['libero'][r]['P95']['uscite']) for r in FR))
    for r in ("run1", "run2", "run8"):
        F.append(f"scopo2 {nome} blocchi {r}: " + " ".join(f"{b[2]:.2f}" for b in P["libero"][r]["P95"]["blocchi"]))

# ---------------- Tabella S3: validita'
righe = []
for sis, n in (("complesso", 5), ("libero", 10)):
    for i in range(1, n + 1):
        R = ct["repliche"][f"{sis}_run{i}"]
        if sis == "complesso":
            x, y = R["50"]["complesso"], R["100"]["complesso"]
            assert x["esito"] == y["esito"] == "valida"
            righe.append([f"Complex {i}", "50–500", f"{x['rmsd_core_media']:.2f}", f"{x['rmsd_core_max']:.2f}", f"{x['immagini_min_A']:.1f}",
                          f"{x['Q_media']:.2f}", f"{x['contatti_persi']:.2f} ({y['contatti_persi']:.2f})", f"{x['centri_crescita_media']:.2f}", f"{x['centri_crescita_max']:.2f}", "valid"])
        else:
            for fin, lab in (("scopo3", "50–500"), ("scopo2", "50–1000")):
                x, y = R["50"][fin], R["100"][fin]
                assert x["esito"] == y["esito"] == "valida"
                righe.append([f"Free {i}", lab, f"{x['rmsd_core_media']:.2f}", f"{x['rmsd_core_max']:.2f}", f"{x['immagini_min_A']:.1f}", "–", "–", "–", "–", "valid"])
EN.append("**Table S3. Validity of the replicas (50-ns discard).** Core RMSD: Cα RMSD of the 43 β-strand residues of MrkA from the crystal (Å); no frame exceeded 2 Å. Image distance: minimum distance between periodic images (Å; limit 12 Å). Q: mean fraction of the 27 crystal contacts present. Lost: fraction of the crystal contacts with occupancy < 0.5 (limit 0.5; in parentheses, 100-ns discard). Fv–core: increase of the distance between the Fv and MrkA core centers relative to the crystal (Å; limit 5 Å in more than half of the frames). All replicas were also valid with the 100-ns discard.\n\n"
          + md(righe, ["Replica", "Window (ns)", "Core RMSD, mean", "max", "Image distance, min", "Q", "Lost", "Fv–core, mean", "max", "Outcome"]))
cc = [ct["repliche"][f"complesso_run{i}"]["50"]["complesso"] for i in range(1, 6)]
l3 = [ct["repliche"][f"libero_run{i}"]["50"]["scopo3"] for i in range(1, 11)]; l2 = [ct["repliche"][f"libero_run{i}"]["50"]["scopo2"] for i in range(1, 11)]
rng = lambda v, f="{:.2f}": f.format(min(v)) + "–" + f.format(max(v))
F.append("validita' scarto 50: RMSD core media complesso " + rng([x["rmsd_core_media"] for x in cc]) + "; libero 50-500 " + rng([x["rmsd_core_media"] for x in l3]) + "; libero 50-1000 " + rng([x["rmsd_core_media"] for x in l2])
         + "; RMSD core max complesso " + rng([x["rmsd_core_max"] for x in cc]) + "; libero " + rng([x["rmsd_core_max"] for x in l2])
         + "; immagini min complesso " + rng([x["immagini_min_A"] for x in cc], "{:.1f}") + "; libero " + rng([x["immagini_min_A"] for x in l2], "{:.1f}")
         + "; Q " + rng([x["Q_media"] for x in cc]) + "; persi " + rng([x["contatti_persi"] for x in cc])
         + "; persi scarto 100 " + rng([ct["repliche"][f"complesso_run{i}"]["100"]["complesso"]["contatti_persi"] for i in range(1, 6)])
         + "; crescita centri media " + rng([x["centri_crescita_media"] for x in cc]) + "; max " + rng([x["centri_crescita_max"] for x in cc])
         + "; riassemblaggio max " + rng([ct["riassemblaggio"][f"complesso_run{i}"]["max_A"] for i in range(1, 6)]))

# ---------------- Tabella S4: RMSF di segmento
righe = []
for seg in ("87-101", "138-147"):
    for sc in ("50", "100"):
        P = s3["per_scarto"][sc]; cs = [P["segmenti"][seg]["complesso"][r] for r in CR]; ls = [P["segmenti"][seg]["libero"][r] for r in FR]
        t = P["test"][seg]
        righe.append([seg.replace("-", "–"), sc, ", ".join(f"{v:.2f}" for v in cs), f"{np.mean(cs):.2f} ± {np.std(cs, ddof=1):.2f}",
                      ", ".join(f"{v:.2f}" for v in ls), f"{np.mean(ls):.2f} ± {np.std(ls, ddof=1):.2f}",
                      f"{'one-sided' if t['alternativa'] == 'less' else 'two-sided'}: t = {t['welch_t']:.2f}, p = {t['welch_p']:.2f}, Holm p = {t['welch_p_holm']:.2f}",
                      f"U = {t['mw_U']:.0f}, p = {t['mw_p']:.2f}, Holm p = {t['mw_p_holm']:.2f}"])
        F.append(f"scopo3 {seg} scarto {sc}: complesso {np.mean(cs):.3f} +- {np.std(cs, ddof=1):.3f} (min {min(cs):.2f}, max {max(cs):.2f}); libero {np.mean(ls):.3f} +- {np.std(ls, ddof=1):.3f} (min {min(ls):.2f}, max {max(ls):.2f}); diff {np.mean(cs) - np.mean(ls):.3f}; "
                 f"welch t {t['welch_t']:.3f} p {t['welch_p']:.4f} holm {t['welch_p_holm']:.4f}; MW U {t['mw_U']:.0f} p {t['mw_p']:.4f} holm {t['mw_p_holm']:.4f}; differenza {t['differenza']}")
EN.append("**Table S4. Mean Cα RMSF of the two segments per replica (Å).** Window: discard to 500 ns in both systems. Complex replicas 1–5 and free replicas 1–10, in order; mean ± SD over replicas. Welch tests: one-sided (complex < free) for 87–101, two-sided for 138–147; Holm correction over the two segments. Mann–Whitney: exact test, same alternatives.\n\n"
          + md(righe, ["Segment", "Discard (ns)", "Complex, per replica", "Complex, mean ± SD", "Free MrkA, per replica", "Free, mean ± SD", "Welch", "Mann–Whitney"]))
P = s3["per_scarto"]["50"]; res = P["resid"]
fuori = sorted({res[i] for k, v in P["rmsf"].items() for i, x in enumerate(v) if x > 4.6})
mx = max((x, res[i], k) for k, v in P["rmsf"].items() for i, x in enumerate(v) if res[i] > 40)
F.append(f"scopo3 RMSF per residuo scarto 50: residui con RMSF > 4.6 A in almeno una replica: {fuori}; massimo fuori da 36-40: {mx[0]:.2f} A al residuo {mx[1]} ({mx[2]})")

# ---------------- Tabella S5: orientazione dell'Fv (17)
righe = []
for i, r in enumerate(CR, 1):
    x = cer["repliche"][r]; f = x["finestra"]; p = x["primo"]
    righe.append([i, f"{p['rot_fv']:.0f}", f"{f['rot_fv']['media']:.0f} ({f['rot_fv']['min']:.0f}–{f['rot_fv']['max']:.0f})", f"{f['tip_fv']['media']:.2f}", f"{f['tip_core']['media']:.2f}",
                  f"{f['epi_epi']['media']:.2f}", f"{f['epi_core']['media']:.2f}", f"{x['corr_epi_core_rot']:.2f}"])
e = cer["equil_s4"]
EN.append("**Table S5. Orientation of the Fv relative to the MrkA core in the complex (descriptive; pending confirmation on the workstation).** Rotation of the Fv (221 Cα atoms, built residues excluded) after superposition on the MrkA core, relative to the crystal: first production frame, and mean (minimum–maximum) over 50–500 ns, one frame per ns. RMSD from the crystal (Å, mean over 50–500 ns) of the hairpin tip (Cα of residues 90–93) after fitting on the Fv or on the core, and of the epitope backbone after fitting on itself or on the core. r: Pearson correlation between the core-fitted epitope RMSD and the rotation over the frames of each replica (one per ns, 50–500 ns); descriptive, because frames within a replica are not independent. "
          f"At the end of the restrained equilibration the rotation was {e['rot_fv']:.1f}° and the core-fitted epitope RMSD {e['epi_core']:.2f} Å.\n\n"
          + md(righe, ["Complex replica", "Rotation, first frame (°)", "Rotation, 50–500 ns (°)", "Tip, Fv fit", "Tip, core fit", "Epitope, self fit", "Epitope, core fit", "r"]))
fin = lambda k, q="media": [cer["repliche"][r]["finestra"][k][q] for r in CR]
F.append("cerniera: rot media " + ", ".join(f"{v:.1f}" for v in fin("rot_fv")) + f"; rot min {min(fin('rot_fv', 'min')):.1f} max {max(fin('rot_fv', 'max')):.1f}; rot primo fotogramma " + ", ".join(f"{cer['repliche'][r]['primo']['rot_fv']:.1f}" for r in CR)
         + f"; rot fine s4 {e['rot_fv']:.2f}; tip|Fv " + ", ".join(f"{v:.2f}" for v in fin("tip_fv")) + "; tip|core " + ", ".join(f"{v:.2f}" for v in fin("tip_core"))
         + "; epi|epi " + ", ".join(f"{v:.2f}" for v in fin("epi_epi")) + "; epi|core " + ", ".join(f"{v:.2f}" for v in fin("epi_core")) + "; tip|Fv max " + ", ".join(f"{v:.2f}" for v in fin("tip_fv", "max"))
         + "; trasl media " + ", ".join(f"{v:.1f}" for v in fin("trasl_fv")))
pr = {int(k): [cer["repliche"][r]["per_residuo_finestra"][k] for r in CR] for k in cer["repliche"]["run1"]["per_residuo_finestra"]}
F.append("cerniera: scostamento Ca dopo fit sul core (50-fine) per residuo, run1-5: " + "; ".join(f"{k} " + ", ".join(f"{v:.1f}" for v in vs) for k, vs in sorted(pr.items())))
tip = [v for k in (90, 91, 92, 93) for v in pr[k]]
F.append(f"cerniera: Ca 90-93 dopo fit sul core, minimo {min(tip):.1f} massimo {max(tip):.1f} A (medie per replica e residuo)")

# ---------------- Tabelle dai descrittivi (d9 interfaccia, d10 permanenza, d11 RMSF oltre 500 ns)
if DESC:
    d9, d10, d11 = (json.load(open(os.path.join(DESC, f))) for f in ("d9_interfaccia.json", "d10_persistenza_scan.json", "d11_verifiche.json"))
    # S2 (continua 2): repliche con cicli in funzione della permanenza richiesta
    perm = ("0.5", "1.0", "1.5", "2.0", "2.5", "3.0", "4.0", "5.0")
    righe = []
    for mis, k in (("Core-fitted (registered)", "prim"), ("Self-fitted (post hoc)", "sec")):
        for sc in ("50", "100"):
            x = d10[f"{k}_{sc}"]
            assert x["repliche_con_cicli"]["1.0"] == conteggi((s2 if k == "prim" else s2["secondaria"])["per_scarto"][sc])["P95"][2]
            righe.append([mis, sc, x["dentro"]] + [x["repliche_con_cicli"][p] for p in perm])
    def cambi(k):
        v = d10[k]["repliche_con_cicli"]; ks = sorted(v, key=float); out = []
        for i in range(1, len(ks)):
            if v[ks[i]] != v[ks[i - 1]]: out.append(f"{ks[i]} ns: {v[ks[i]]}")
        return ", ".join(out)
    F.append("permanenza (post hoc): primaria scarto 50, cambi del conteggio: " + cambi("prim_50") + " | scarto 100: " + cambi("prim_100"))
    F.append("permanenza (post hoc): secondaria scarto 50: " + cambi("sec_50") + " | scarto 100: " + cambi("sec_100"))
    F.append("permanenza (post hoc): tratto di fila fuori piu' lungo, primaria scarto 50: " + ", ".join(f"{r} {v:.1f}" for r, v in d10["prim_50"]["fuori_di_fila_max_ns"].items())
             + " | scarto 100: " + ", ".join(f"{r} {v:.1f}" for r, v in d10["prim_100"]["fuori_di_fila_max_ns"].items()))
    F.append("permanenza (post hoc): ultimo rientro (ns), primaria scarto 50: " + ", ".join(f"{r} {v:.1f}" for r, v in d10["prim_50"]["ultimo_rientro_ns"].items() if v is not None)
             + f"; cicli nella sola finestra 500-1000 ns: primaria {d10['prim_50']['cicli_500_1000']}, secondaria {d10['sec_50']['cicli_500_1000']}")
    EN.append("**Table S2 (continued 2). Free replicas with at least one exit and re-entry against the persistence required for an exit or a re-entry (post hoc, 12 Nov 2026).** "
              "Pooled 95th-percentile envelope. The registered persistence is 1.0 ns. Inside: replicas with a fraction of frames inside the envelope ≥ 0.5, which does not depend on the persistence. "
              "With 8 or more replicas inside, 5 or more replicas with a cycle read as populated and fewer than 5 as indistinguishable from memory of the starting structure. Self-fitted values pending confirmation on the workstation.\n\n"
              + md(righe, ["Measure", "Discard (ns)", "Inside"] + [p + " ns" + (" (registered)" if p == "1.0" else "") for p in perm]))
    # S4 (continua): RMSF del libero oltre la durata confermativa
    righe = []
    for seg in ("87-101", "138-147"):
        for fin in ("50-500", "550-1000", "50-1000"):
            v = d11["rmsf_libero"][seg][fin]
            if fin == "50-500":
                assert max(abs(x - s3["per_scarto"]["50"]["segmenti"][seg]["libero"][r]) for x, r in zip(v, FR)) < 1e-6
            righe.append([seg.replace("-", "–"), fin.replace("-", "–"), ", ".join(f"{x:.2f}" for x in v), f"{np.mean(v):.2f} ± {np.std(v, ddof=1):.2f}", f"{min(v):.2f}–{max(v):.2f}"])
            F.append(f"RMSF libero {seg} {fin}: media {np.mean(v):.3f} DS {np.std(v, ddof=1):.3f} min {min(v):.2f} max {max(v):.2f}")
    EN.append("**Table S4 (continued). Mean Cα RMSF of the two segments in free MrkA beyond the first 500 ns (Å).** A descriptive analysis planned in the criteria, without a test: the complex has no data beyond 500 ns. "
              "The window 550–1000 ns has the same length as the window of the comparison (50–500 ns). Free replicas 1–10, in order; mean ± SD and range over replicas.\n\n"
              + md(righe, ["Segment", "Window (ns)", "Free MrkA, per replica", "Mean ± SD", "Range"]))
    F.append("rilevabile (A): " + ", ".join(f"{k} {v:.3f}" for k, v in d11["rilevabile"].items()))
    F.append("RMSF Ca 36-41, media sulle 15 repliche: " + ", ".join(f"{k} {v:.2f}" for k, v in d11["rmsf_36_41_media"].items()) + f"; residuo 40 min-max {d11['rmsf_40_minmax'][0]:.2f}-{d11['rmsf_40_minmax'][1]:.2f}")
    # S6: interfaccia
    S = d9["scarti"]["50"]; cr = d9["cristallo"]; D = d9["definizioni"]
    righe = []
    hb = sorted(((k, v) for k, v in S["hb"].items() if max(v) >= 0.10), key=lambda x: -np.mean(x[1]))
    for k, v in hb:
        ra, rf = nomeEN(k); dc = cr["polari"].get(k)
        at = sorted(((kk, vv) for kk, vv in S["hb_atomi"].items() if kk.startswith(k + "|")), key=lambda x: -np.mean(x[1]))[0][0].split("|")
        atomi = f"{at[3]} → {at[4]}" if at[5] == "A" else f"{at[4]} ← {at[3]}"
        righe.append(["hydrogen bond", ra, rf, atomi, f"{dc:.2f}" if dc is not None else "> 4.5"] + [f"{x:.2f}" for x in v])
    sb = sorted(((k, v) for k, v in S["sb"].items() if max(v) >= 0.10), key=lambda x: -np.mean(x[1]))
    for k, v in sb:
        ra, rf = nomeEN(k); dc = cr["salini"].get(k)
        righe.append(["salt bridge", ra, rf, "–", f"{dc:.2f}" if dc is not None else "–"] + [f"{x:.2f}" for x in v])
    bs = [S["bsa"][r]["media"] for r in CR]
    righe.append(["buried surface (Å²)", "–", "–", "–", f"{cr['bsa']:.0f}"] + [f"{x:.0f}" for x in bs])
    EN.append("**Table S6. Hydrogen bonds, salt bridges and buried surface at the MrkA–Fab interface per complex replica (50–500 ns).** A descriptive analysis planned in the criteria and run on 12 Nov 2026, after the outcomes; "
              "its geometric definitions were not in the criteria and were fixed in the script before the numbers were seen. No outcome depends on it. "
              f"Hydrogen bond: donor N or O with a bound hydrogen, acceptor any O or an N without bound hydrogen, donor–acceptor distance ≤ {D['D_HB']} Å and hydrogen–donor–acceptor angle ≤ {D['ANG_HB']:.0f}°; "
              "occupancy per residue pair, with the most frequent atom pair (MrkA atom first; the arrow points from donor to acceptor). "
              f"Salt bridge: Lys NZ or Arg NE, NH1 or NH2 within {D['D_SB']} Å of a carboxylate oxygen of Asp or Glu. Pairs with an occupancy of at least 0.10 in at least one replica are listed. "
              "Crystal: minimum N/O–N/O distance of the pair in the completed model, which has no hydrogens (Å; for salt bridges, the N–O distance of the charged groups). "
              f"Buried surface: Shrake–Rupley, probe {D['sonda']} Å, heavy atoms with Bondi radii, sum over the two partners, one frame per ns; mean per replica.\n\n"
              + md(righe, ["Quantity", "MrkA", "Fab", "Atoms", "Crystal", "Replica 1", "2", "3", "4", "5"]))
    F.append(f"interfaccia: BSA cristallo {cr['bsa']:.0f} A^2; per replica " + ", ".join(f"{x:.0f}" for x in bs) + "; rapporto " + ", ".join(f"{x / cr['bsa']:.2f}" for x in bs)
             + "; DS sui fotogrammi " + ", ".join(f"{S['bsa'][r]['ds']:.0f}" for r in CR))
    F.append("interfaccia: legami H con occupanza >= 0.5 in >= 4 repliche: " + ", ".join(coppia(k) + " (" + "-".join(f"{f(v):.2f}" for f in (min, max)) + ")" for k, v in hb if sum(x >= 0.5 for x in v) >= 4))
    F.append("interfaccia: ponti salini: " + "; ".join(coppia(k) + " " + ", ".join(f"{x:.2f}" for x in v) for k, v in sb))
    F.append("interfaccia: coppie di residui con legame H per fotogramma, somma delle occupanze per replica: " + ", ".join(f"{sum(v[i] for v in S['hb'].values()):.1f}" for i in range(5)))

prov = "sha256 (16 cifre) dei file letti: " + ", ".join(f"{f} {sha(f)}" for f in ("scopo1_contatti.json", "scopo2_epitopo.json", "scopo3_rmsf.json", "controlli.json", "diagnostica_cerniera.json")) \
       + f"; pdb {hashlib.sha256(open(PDB, 'rb').read()).hexdigest()[:16]}"
os.makedirs(OUT, exist_ok=True)
# nelle tabelle il segno meno e' quello tipografico (U+2212), non il trattino; fatti.txt resta in ASCII
meno = lambda t: re.sub(r"(?<=[ (=|])-(?=\d)", "\u2212", t)
def _ordine(t):
    """Ordina le tabelle per numero (S1, S1 continued, S2, S2 continued, S2 continued 2, ...), comunque siano state aggiunte."""
    import re
    m = re.match(r"\*\*Table S(\d+)(?: \(continued(?: (\d+))?\))?\.", t)
    assert m, t[:60]
    return (int(m.group(1)), 0 if "(continued" not in t[:40] else int(m.group(2) or 1))
EN.sort(key=_ordine)
open(os.path.join(OUT, "tabelle_EN.md"), "w").write(meno("\n\n".join(EN)) + "\n")
open(os.path.join(OUT, "tabelle_IT.md"), "w").write(meno("\n\n".join(IT)) + "\n")
open(os.path.join(OUT, "fatti.txt"), "w").write(prov + "\n\n" + "\n".join(F) + "\n")
print(prov); print(); print("\n".join(F))
