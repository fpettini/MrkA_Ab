"""Scopo 1: contatti dell'epitopo in soluzione (criteri, sezione "Scopo 1: contatti dell'epitopo in soluzione").

Uso (dalla cartella complesso_fab, dopo 11_controlli.py):
    python3 analisi/12_scopo1_contatti.py [--cache analisi/cache] [--out analisi/risultati] [--pdb mrka-fab_completo.pdb]

Contatto = coppia di residui MrkA-Fab con almeno una coppia di atomi pesanti entro 4.0 A (sensibilita' 3.5 e 5.0).
Per ogni replica valida del complesso, finestra dallo scarto (50 ns; 100 di sensibilita') a T3: frazione dei
fotogrammi in contatto. Con N repliche valide e k = ceil(0.8 N):
  persistente = nel cristallo, frazione >= 0.5 in almeno k repliche;  perso = nel cristallo, < 0.5 in almeno k;
  nuovo = non nel cristallo, >= 0.5 in almeno k;  le altre coppie restano senza esito (solo descrittive).
Un esito vale solo se e' lo stesso con scarto 50 e 100 ns, altrimenti "non robusto". Con meno di 4 repliche
valide: "non concludibile". Riferimento: contatti del cristallo completato con lo stesso criterio (a 4.0 A sono
gli stessi del cristallo originale: 27 coppie, verificato il 5 nov 2026).
Analisi INTERMEDIA (condizioni del revisore, 5 nov): solo le frazioni per replica; nessun esito (persistente, perso,
nuovo, hotspot) e nessuna soglia di repliche.
Scrive OUT/scopo1_contatti.txt, .json e .csv (tutte le coppie, tutte le distanze e tutti gli scarti).
"""
import os, json, csv, argparse
import numpy as np
import comune_fab as C

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="analisi/cache")
ap.add_argument("--out", default="analisi/risultati")
ap.add_argument("--pdb", default="mrka-fab_completo.pdb")
arg = ap.parse_args()

ctrl = C.leggi_controlli(arg.out)
serie = C.leggi_cache(arg.cache)
T3 = ctrl["T3_ns"]; stato = ctrl["stato"]["scopo1"]
FIN = stato == "FINALE"
NO = "-"                      # nessun esito nell'analisi intermedia
_, rn = C.leggi_pdb(arg.pdb)
nome = lambda c, k: f"{c}{k}{rn[c][k]}"


def frazioni(d, w, cut):
    """{(resA, catena, resF): frazione dei fotogrammi della finestra con distanza minima <= cut}"""
    fr = np.where(w)[0]
    sel = np.isin(d["c_f"], fr) & (d["c_d"] <= cut)
    chiave = (d["c_a"][sel].astype(np.int64) * 2 + d["c_c"][sel]) * 1000 + d["c_k"][sel]
    chiavi, conta = np.unique(chiave, return_counts=True)
    return {(int(x // 2000), "LH"[int(x // 1000 % 2)], int(x % 1000)): float(n) / len(fr) for x, n in zip(chiavi, conta)}


def esito(cr, f, k, n):
    if not FIN:
        return NO
    if n < C.MIN_VALIDE["complesso"]:
        return "non concludibile"
    ge = sum(x >= 0.5 for x in f); lt = len(f) - ge
    if cr:
        return "persistente" if ge >= k else ("perso" if lt >= k else "senza esito")
    return "nuovo" if ge >= k else "senza esito"


ris = {}
for cut in C.CONTATTI:
    cris = C.contatti_cristallo(arg.pdb, cut)
    ris[cut] = {}
    for sc in C.SCARTI:
        v = C.valide(ctrl, "complesso", "complesso", sc)
        reps = v["valide"]; n = len(reps); k = (C.soglia_repliche(n) if n else 0) if FIN else None
        fr = {r: frazioni(serie[("complesso", r)], C.finestra(serie[("complesso", r)]["t_ns"], sc, T3), cut) for r in reps}
        coppie = sorted(set(cris).union(*[set(x) for x in fr.values()]))
        ris[cut][sc] = {"repliche": reps, "k": k, "coppie": {}}
        for p in coppie:
            f = [fr[r].get(p, 0.0) for r in reps]
            ris[cut][sc]["coppie"][p] = {"cristallo": p in cris, "frazioni": f, "n_ge_0.5": int(sum(x >= 0.5 for x in f)),
                                         "esito": esito(p in cris, f, k, n)}


def finale(cut, p):
    if not FIN:
        return NO
    e = [ris[cut][sc]["coppie"][p]["esito"] if p in ris[cut][sc]["coppie"] else "senza esito" for sc in C.SCARTI]
    return e[0] if e[0] == e[1] else f"non robusto ({e[0]} a 50 ns, {e[1]} a 100 ns)"


P = C.prov(__file__, [v["_file"] for (s, _), v in serie.items() if s == "complesso"] + [arg.pdb, os.path.join(arg.out, "controlli.json")],
           {"cache": arg.cache, "out": arg.out, "pdb": arg.pdb}, {"CACHE": arg.cache, "OUT": arg.out})
L = C.testo_prov(P) + C.avviso_intermedia(stato)
L += ["", f"Scopo 1 ({stato}): finestra dallo scarto a T3 = {T3:.1f} ns; contatto = atomi pesanti entro la distanza (<=)."]
for cut in C.CONTATTI:
    r50 = ris[cut][50.0]
    L.append("")
    L.append(f"=== {cut} A{' (principale)' if cut == C.CONTATTO_PRINCIPALE else ' (sensibilita)'}: repliche valide {r50['repliche']}"
             + (f", soglia {r50['k']}/{len(r50['repliche'])}" if FIN else ""))
    tutte = sorted(set(ris[cut][50.0]["coppie"]) | set(ris[cut][100.0]["coppie"]))
    fin = {p: finale(cut, p) for p in tutte}
    cont = {}
    for p, e in fin.items():
        cont[e.split(" (")[0]] = cont.get(e.split(" (")[0], 0) + 1
    ncr = sum(1 for p in tutte if ris[cut][50.0]["coppie"].get(p, {}).get("cristallo"))
    if FIN:
        L.append(f"Coppie nel cristallo: {ncr}. Esiti: " + ", ".join(f"{k} {v}" for k, v in sorted(cont.items())))
    else:
        L.append(f"Coppie nel cristallo: {ncr}. Per replica (scarto 50 ns), contatti del cristallo con frazione >= 0.5: "
                 + ", ".join(f"{r} {sum(1 for p, c in r50['coppie'].items() if c['cristallo'] and c['frazioni'][i] >= 0.5)}"
                             for i, r in enumerate(r50["repliche"])))
    if cut != C.CONTATTO_PRINCIPALE:
        continue
    L.append("Contatti del cristallo, frazione per replica (scarto 50 ns)" + (" ed esito:" if FIN else ":"))
    L.append(f"  {'coppia':28s}" + "".join(f"{r:>7s}" for r in r50["repliche"]) + ("   esito" if FIN else ""))
    for p in tutte:
        c = r50["coppie"].get(p)
        if c and c["cristallo"]:
            L.append(f"  {nome('A', p[0]) + ' - ' + nome(p[1], p[2]):28s}" + "".join(f"{x:7.2f}" for x in c["frazioni"]) + (f"   {fin[p]}" if FIN else ""))
    L.append("Coppie non nel cristallo con frazione >= 0.5 in almeno una replica (scarto 50 ns):")
    nuovi = [p for p in tutte if p in r50["coppie"] and not r50["coppie"][p]["cristallo"] and r50["coppie"][p]["n_ge_0.5"] > 0]
    for p in nuovi:
        c = r50["coppie"][p]
        L.append(f"  {nome('A', p[0]) + ' - ' + nome(p[1], p[2]):28s}" + "".join(f"{x:7.2f}" for x in c["frazioni"]) + (f"   {fin[p]}" if FIN else ""))
    if not nuovi:
        L.append("  nessuna")
    if FIN:
        L.append("Per residuo dell'epitopo (4.0 A, scarto 50 ns): contatti del cristallo persistenti / totali")
    else:
        L.append("Per residuo dell'epitopo (4.0 A, scarto 50 ns): contatti del cristallo con frazione >= 0.5, per replica / totali")
    for a in C.EPITOPO:
        cr = [p for p in tutte if p[0] == a and r50["coppie"].get(p, {}).get("cristallo")]
        if FIN:
            per = [p for p in cr if fin[p] == "persistente"]
            L.append(f"  {nome('A', a):10s} {len(per)}/{len(cr)}")
        else:
            L.append(f"  {nome('A', a):10s} " + ", ".join(f"{r} {sum(1 for p in cr if r50['coppie'][p]['frazioni'][i] >= 0.5)}/{len(cr)}"
                                                       for i, r in enumerate(r50["repliche"])))
L += ["", "Note dal documento: interfaccia stabile partendo dal cristallo = attesa, pesa poco; i contatti del cristallo",
      "possono essere sottostimati (atomi di un altro monomero o del tag rimossi dopo il raffinamento)."]
open(os.path.join(arg.out, "scopo1_contatti.txt"), "w").write("\n".join(L) + "\n")

with open(os.path.join(arg.out, "scopo1_contatti.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["distanza_A", "scarto_ns", "resA", "catena_fab", "resF", "coppia", "cristallo", "repliche", "frazioni", "n_ge_0.5", "soglia", "esito_scarto", "esito_finale"])
    for cut in C.CONTATTI:
        for sc in C.SCARTI:
            for p, c in ris[cut][sc]["coppie"].items():
                w.writerow([cut, int(sc), p[0], p[1], p[2], f"{nome('A', p[0])}-{nome(p[1], p[2])}", int(c["cristallo"]),
                            " ".join(ris[cut][sc]["repliche"]), " ".join(f"{x:.3f}" for x in c["frazioni"]), c["n_ge_0.5"],
                            ris[cut][sc]["k"], c["esito"], finale(cut, p)])
js = {"provenienza": P, "stato": stato, "T3_ns": T3,
      "risultati": {str(cut): {str(int(sc)): {"repliche": v["repliche"], "soglia": v["k"],
                                              "coppie": {f"{p[0]}|{p[1]}|{p[2]}": c for p, c in v["coppie"].items()}}
                               for sc, v in x.items()} for cut, x in ris.items()}}
json.dump(js, open(os.path.join(arg.out, "scopo1_contatti.json"), "w"), indent=1)
print("\n".join(L))
