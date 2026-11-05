"""Scopo 2: la conformazione legata dell'epitopo e' popolata nel libero? (criteri, sezione "Scopo 2").

Uso (dalla cartella complesso_fab, dopo 11_controlli.py):
    python3 analisi/13_scopo2_epitopo.py [--cache analisi/cache] [--out analisi/risultati]

Misura: RMSD del backbone (N, CA, C, O) dell'epitopo (wt 42, 86, 88, 90-94, 96, 167, 168) dal cristallo, dopo fit
sui 43 Ca di strand (calcolata da 10_serie.py, uguale nei due sistemi).
Soglia (involucro del legato): 95 percentile della stessa RMSD sulle repliche valide del complesso unite, finestra
dallo scarto a T3; sensibilita' 90 e 99 percentile e 95 percentile per singola replica; segnalata se > 2 A.
Per ogni replica valida del libero, finestra dallo scarto a T2: frazione dei fotogrammi entro soglia (RMSD <= soglia),
uscite e rientri (almeno 1 ns di fila fuori / dentro; i tratti piu' corti non cambiano lo stato), decorso a blocchi.
Esiti (k = ceil(0.8 x valide), con 10 valide 8): popolata = frazione >= 0.5 in >= k repliche e almeno un ciclo
uscita-rientro in >= 5 repliche; non distinguibile dalla memoria della partenza = frazione >= 0.5 in >= k ma cicli
in meno di 5; non popolata = frazione <= 0.05 in >= k; intermedia = il resto; con meno di 8 valide: non concludibile.
Vale solo se l'esito e' lo stesso con scarto 50 e 100 ns, altrimenti "non robusto".
NON NEL DOC: percentile con interpolazione lineare (numpy); stato iniziale = primo tratto di almeno 1 ns; un
"ingresso" senza uscita precedente nella finestra non conta come rientro; blocchi del decorso da 50 ns.
Scrive OUT/scopo2_epitopo.txt e .json.
"""
import os, json, argparse
import numpy as np
import comune_fab as C

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="analisi/cache")
ap.add_argument("--out", default="analisi/risultati")
arg = ap.parse_args()

ctrl = C.leggi_controlli(arg.out)
serie = C.leggi_cache(arg.cache)
T3, T2 = ctrl["T3_ns"], ctrl["T2_ns"]; stato = ctrl["stato"]["scopo2"]


def transizioni(dentro, nmin):
    """Stato con isteresi: cambia solo con almeno nmin fotogrammi di fila dall'altra parte."""
    stato = None; iniziale = None; uscite = rientri = ingressi = 0
    i = 0; n = len(dentro)
    while i < n:
        j = i
        while j < n and dentro[j] == dentro[i]:
            j += 1
        v = bool(dentro[i])
        if j - i >= nmin:
            if stato is None:
                stato = iniziale = v
            elif v != stato:
                if not v:
                    uscite += 1
                elif uscite > rientri:
                    rientri += 1
                else:
                    ingressi += 1
                stato = v
        i = j
    return {"stato_iniziale": None if iniziale is None else ("dentro" if iniziale else "fuori"),
            "uscite": uscite, "rientri": rientri, "ingressi_senza_uscita": ingressi, "cicli": rientri}


def esito(fr, cicli, n, k):
    if n < C.MIN_VALIDE["libero"]:
        return "non concludibile"
    ge = sum(f >= 0.5 for f in fr); le = sum(f <= 0.05 for f in fr); cy = sum(c >= 1 for c in cicli)
    if ge >= k:
        return "popolata" if cy >= C.CICLI_MIN else "non distinguibile dalla memoria della partenza"
    if le >= k:
        return "non popolata"
    return "intermedia"


ris = {}
for sc in C.SCARTI:
    vc = C.valide(ctrl, "complesso", "complesso", sc)["valide"]
    vl = C.valide(ctrl, "libero", "scopo2", sc)["valide"]
    rc = {r: serie[("complesso", r)]["rmsd_epi"][C.finestra(serie[("complesso", r)]["t_ns"], sc, T3)] for r in vc}
    if not rc:
        raise SystemExit(f"scarto {sc:.0f} ns: nessuna replica valida del complesso, soglia non calcolabile")
    tutto = np.concatenate(list(rc.values()))
    soglie = {"P95": float(np.percentile(tutto, 95)), "P90": float(np.percentile(tutto, 90)), "P99": float(np.percentile(tutto, 99))}
    per_rep = {r: float(np.percentile(x, 95)) for r, x in rc.items()}
    soglie["P95_replica_min"] = min(per_rep.values()); soglie["P95_replica_max"] = max(per_rep.values())
    n = len(vl); k = C.soglia_repliche(n) if n else 0
    lib = {}
    for r in vl:
        d = serie[("libero", r)]
        w = C.finestra(d["t_ns"], sc, T2)
        x = d["rmsd_epi"][w]; t = d["t_ns"][w]
        nmin = int(round(C.ISTERESI_NS * 1000.0 / float(d["dt_ps"])))
        lib[r] = {"rmsd_media": float(x.mean()), "rmsd_mediana": float(np.median(x)), "fotogrammi": int(len(x)), "n_isteresi": nmin}
        for nome, s in soglie.items():
            dentro = x <= s
            blocchi = []
            b0 = sc
            while b0 + C.BLOCCO_SCOPO2_NS <= T2 + 1e-6:
                wb = (t >= b0 - 1e-6) & (t < b0 + C.BLOCCO_SCOPO2_NS - 1e-6)
                blocchi.append([b0, b0 + C.BLOCCO_SCOPO2_NS, float(dentro[wb].mean()) if wb.any() else None])
                b0 += C.BLOCCO_SCOPO2_NS
            lib[r][nome] = {"frazione": float(dentro.mean()), "blocchi": blocchi, **transizioni(dentro, nmin)}
    es = {nome: esito([lib[r][nome]["frazione"] for r in vl], [lib[r][nome]["cicli"] for r in vl], n, k) for nome in soglie}
    ris[sc] = {"complesso_valide": vc, "libero_valide": vl, "k": k, "soglie": soglie, "soglia_per_replica": per_rep,
               "complesso_media": {r: float(x.mean()) for r, x in rc.items()}, "libero": lib, "esiti": es}

fin = {}
for nome in ris[50.0]["soglie"]:
    e = [ris[sc]["esiti"][nome] for sc in C.SCARTI]
    fin[nome] = e[0] if e[0] == e[1] else f"non robusto ({e[0]} a 50 ns, {e[1]} a 100 ns)"

P = C.prov(__file__, [v["_file"] for v in serie.values()] + [os.path.join(arg.out, "controlli.json")],
           {"cache": arg.cache, "out": arg.out}, {"CACHE": arg.cache, "OUT": arg.out})
L = C.testo_prov(P) + C.avviso_intermedia(stato)
L += ["", f"Scopo 2 ({stato}): soglia dal complesso fino a T3 = {T3:.1f} ns; libero fino a T2 = {T2:.1f} ns."]
for sc in C.SCARTI:
    r = ris[sc]; s = r["soglie"]
    L += ["", f"=== Scarto {sc:.0f} ns: complesso valide {r['complesso_valide']}, libero valide {r['libero_valide']} (soglia {r['k']})",
          f"Soglia P95 = {s['P95']:.3f} A (P90 {s['P90']:.3f}, P99 {s['P99']:.3f}; P95 per replica {s['P95_replica_min']:.3f}-{s['P95_replica_max']:.3f})"
          + ("  ATTENZIONE: soglia > 2 A" if s["P95"] > C.SOGLIA_EPI_FLAG else ""),
          "RMSD media nel complesso: " + ", ".join(f"{x} {v:.2f}" for x, v in r["complesso_media"].items()),
          f"{'libero':8s}{'RMSD med':>9s}{'fraz P95':>9s}{'usc':>5s}{'rientri':>8s}{'ingr':>5s}{'iniz':>7s}{'fraz P90':>9s}{'fraz P99':>9s}   blocchi P95 (da {C.BLOCCO_SCOPO2_NS:.0f} ns)"]
    for x, v in r["libero"].items():
        p = v["P95"]
        L.append(f"{x:8s}{v['rmsd_media']:9.2f}{p['frazione']:9.2f}{p['uscite']:5d}{p['rientri']:8d}{p['ingressi_senza_uscita']:5d}"
                 f"{(p['stato_iniziale'] or '-'):>7s}{v['P90']['frazione']:9.2f}{v['P99']['frazione']:9.2f}   "
                 + " ".join(f"{b[2]:.2f}" if b[2] is not None else "-" for b in p["blocchi"]))
    L.append("Esiti a questo scarto: " + "; ".join(f"{k} {v}" for k, v in r["esiti"].items()))
L += ["", "Esito (stesso a 50 e 100 ns): " + "; ".join(f"{k}: {v}" for k, v in fin.items()),
      "", "Da riportare con i risultati (documento dei criteri): l'acqua interna E206 sta nel cuore dell'epitopo",
      "(N di Asp86 e Cys87, O di Asp168); OD1 di Asp40 lega l'N di Ala90; l'epitopo circonda il ponte Cys43-Cys87 e",
      "le MD sono solo ossidate. Formulazioni: mai 'selezione conformazionale' ne' 'induced fit'."]
open(os.path.join(arg.out, "scopo2_epitopo.txt"), "w").write("\n".join(L) + "\n")
json.dump({"provenienza": P, "stato": stato, "T3_ns": T3, "T2_ns": T2, "per_scarto": {str(int(k)): v for k, v in ris.items()},
           "esito": fin}, open(os.path.join(arg.out, "scopo2_epitopo.json"), "w"), indent=1)
print("\n".join(L))
