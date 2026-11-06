"""Scopo 2: la conformazione legata dell'epitopo e' popolata nel libero? (criteri, sezione "Scopo 2").

Uso (dalla cartella complesso_fab, dopo 11_controlli.py):
    python3 analisi/13_scopo2_epitopo.py [--cache analisi/cache] [--out analisi/risultati]
                                         [--trj analisi/traiettorie] [--pdb mrka-fab_completo.pdb] [--controllo run1]

Misura PRIMARIA (doc): RMSD del backbone (N, CA, C, O) dell'epitopo (wt 42, 86, 88, 90-94, 96, 167, 168) dal
cristallo, dopo fit sui 43 Ca di strand (calcolata da 10_serie.py, uguale nei due sistemi).
Soglia (involucro del legato): 95 percentile della stessa RMSD sulle repliche valide del complesso unite, finestra
dallo scarto a T3; sensibilita' 90 e 99 percentile e 95 percentile per singola replica; segnalata se > 2 A.
Per ogni replica valida del libero, finestra dallo scarto a T2: frazione dei fotogrammi entro soglia (RMSD <= soglia),
uscite e rientri (almeno 1 ns di fila fuori / dentro; i tratti piu' corti non cambiano lo stato), decorso a blocchi.
Esiti (k = ceil(0.8 x valide), con 10 valide 8; c = ceil(0.5 x valide), con 10 valide 5): popolata = frazione >= 0.5
in >= k repliche e almeno un ciclo uscita-rientro in >= c repliche; non distinguibile dalla memoria della partenza =
frazione >= 0.5 in >= k ma cicli in meno di c; non popolata = frazione <= 0.05 in >= k; intermedia = il resto; con
meno di 8 valide: non concludibile. Vale solo se l'esito e' lo stesso con scarto 50 e 100 ns, altrimenti "non robusto".
NON NEL DOC (decise dal revisore il 5 nov): percentile con interpolazione lineare (numpy); stato iniziale = primo
tratto di almeno 1 ns; un "ingresso" senza uscita precedente nella finestra non conta come rientro; blocchi del
decorso da 50 ns; c = ceil(0.5 x valide) invece di 5 fisse.

Misura SECONDARIA (aggiunta il 6 nov 2026 dopo l'intermedia a 250 ns, per la rotazione Fv-core; proposta del
revisore, approvata da francesco; registro del doc): RMSD del backbone degli stessi 11 residui dopo fit sull'epitopo
stesso (comune_fab.rmsd_epitopo_proprio, unica definizione), con tutte le regole della primaria. Si riportano
entrambe; se gli esiti differiscono si scrive il disaccordo, senza scegliere. Condizioni del revisore (6 nov):
 - complesso: dalla cache di 10; controllo: ricalcolata qui dall'xtc di 09 di una replica del complesso con la stessa
   funzione deve coincidere con la cache (scarto < 1e-6 A, stampato nell'uscita);
 - libero: si calcola (dagli xtc di 09, senza cache su disco) SOLO con lo stato FINALE determinato dalle lunghezze
   (5 x 500 ns e 10 x 1 us); se lo stato FINALE viene da --finale di 11 prima di allora, 13 si ferma con un errore;
 - nell'intermedia si stampano solo la soglia secondaria e le medie del complesso.
Analisi INTERMEDIA (condizioni del revisore, 5 nov): soglia e numeri per replica, nessun esito.
Scrive OUT/scopo2_epitopo.txt e .json.
"""
import os, json, argparse
import numpy as np
import comune_fab as C

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="analisi/cache")
ap.add_argument("--out", default="analisi/risultati")
ap.add_argument("--trj", default="analisi/traiettorie", help="xtc di 09 e proteina_<sistema>.gro (secondaria)")
ap.add_argument("--pdb", default="mrka-fab_completo.pdb")
ap.add_argument("--controllo", default=None, help="replica del complesso per il controllo della secondaria (default: la prima)")
arg = ap.parse_args()
TOLL_CONTROLLO = 1e-6

ctrl = C.leggi_controlli(arg.out)
serie = C.leggi_cache(arg.cache)
T3, T2 = ctrl["T3_ns"], ctrl["T2_ns"]; stato = ctrl["stato"]["scopo2"]
FIN = stato == "FINALE"
_, _, lung = C.durate(serie)
FIN_LUNG = C.finale_per_lunghezze(lung)
if FIN and not FIN_LUNG:
    raise SystemExit("stato FINALE da --finale, ma non tutte le repliche del disegno sono alla durata massima: la misura "
                     "secondaria del libero si sblocca solo con FINALE per lunghezze (revisore, 6 nov 2026). 13 si ferma.")
mancano = [r for (s, r), d in serie.items() if s == "complesso" and "rmsd_epi_proprio" not in d]
if mancano:
    raise SystemExit(f"cache del complesso senza la misura secondaria ({', '.join(sorted(mancano))}): rilanciare 10_serie.py")


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
    if not FIN:
        return "-"
    if n < C.MIN_VALIDE["libero"]:
        return "non concludibile"
    ge = sum(f >= 0.5 for f in fr); le = sum(f <= 0.05 for f in fr); cy = sum(c >= 1 for c in cicli)
    if ge >= k:
        return "popolata" if cy >= C.soglia_cicli(n) else "non distinguibile dalla memoria della partenza"
    if le >= k:
        return "non popolata"
    return "intermedia"


def analizza(misura_c, misura_l):
    """misura_c(r) / misura_l(r): serie completa (allineata a t_ns della cache) di una replica del complesso / libero.
    misura_l = None: solo soglia e medie del complesso (secondaria nell'intermedia)."""
    ris = {}
    for sc in C.SCARTI:
        vc = C.valide(ctrl, "complesso", "complesso", sc)["valide"]
        vl = C.valide(ctrl, "libero", "scopo2", sc)["valide"] if misura_l else []
        rc = {r: misura_c(r)[C.finestra(serie[("complesso", r)]["t_ns"], sc, T3)] for r in vc}
        if not rc:
            raise SystemExit(f"scarto {sc:.0f} ns: nessuna replica valida del complesso, soglia non calcolabile")
        tutto = np.concatenate(list(rc.values()))
        soglie = {"P95": float(np.percentile(tutto, 95)), "P90": float(np.percentile(tutto, 90)), "P99": float(np.percentile(tutto, 99))}
        per_rep = {r: float(np.percentile(x, 95)) for r, x in rc.items()}
        soglie["P95_replica_min"] = min(per_rep.values()); soglie["P95_replica_max"] = max(per_rep.values())
        n = len(vl); k = (C.soglia_repliche(n) if n else 0) if FIN else None
        kc = (C.soglia_cicli(n) if n else 0) if FIN else None
        lib = {}
        for r in vl:
            d = serie[("libero", r)]
            w = C.finestra(d["t_ns"], sc, T2)
            x = misura_l(r)[w]; t = d["t_ns"][w]
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
        es = {nome: esito([lib[r][nome]["frazione"] for r in vl], [lib[r][nome]["cicli"] for r in vl], n, k) for nome in soglie} \
            if misura_l else None
        ris[sc] = {"complesso_valide": vc, "libero_valide": vl, "k": k, "k_cicli": kc, "soglie": soglie, "soglia_per_replica": per_rep,
                   "complesso_media": {r: float(x.mean()) for r, x in rc.items()}, "libero": lib, "esiti": es}
    fin = None
    if FIN and misura_l:
        fin = {}
        for nome in ris[C.SCARTI[0]]["soglie"]:
            e = [ris[sc]["esiti"][nome] for sc in C.SCARTI]
            fin[nome] = e[0] if e[0] == e[1] else f"non robusto ({e[0]} a 50 ns, {e[1]} a 100 ns)"
    return ris, fin


def xtc_della_cache(d):
    P = json.loads(str(d["provenienza"]))
    x = [k for k in P["input"] if k.endswith(".xtc")]
    if len(x) != 1 or not x[0].startswith("<TRJ>/"):
        raise SystemExit(f"{d['_file']}: xtc della serie non riconosciuto nella provenienza")
    f = os.path.join(arg.trj, x[0][len("<TRJ>/"):])
    if not os.path.exists(f) or C.sha256(f) != P["input"][x[0]]:
        raise SystemExit(f"{f}: assente o diverso da quello usato per la cache (rilanciare 10)")
    return f


def secondaria_da_xtc(s, r):
    d = serie[(s, r)]
    t, x = C.serie_epitopo_proprio(os.path.join(arg.trj, f"proteina_{s}.gro"), xtc_della_cache(d), s, arg.pdb)
    if len(t) != len(d["t_ns"]) or np.abs(t - d["t_ns"]).max() > 1e-3:
        raise SystemExit(f"{s} {r}: fotogrammi dell'xtc diversi da quelli della cache")
    return x


# --- primaria
prim, fin_p = analizza(lambda r: serie[("complesso", r)]["rmsd_epi"], lambda r: serie[("libero", r)]["rmsd_epi"])

# --- secondaria: controllo sul complesso (stessa funzione di 10, ricalcolata dall'xtc)
rep_c = sorted((r for s, r in serie if s == "complesso"), key=lambda r: int(r[3:]))
rc0 = arg.controllo or rep_c[0]
if ("complesso", rc0) not in serie:
    raise SystemExit(f"--controllo {rc0}: replica del complesso assente nella cache")
print(f"   secondaria, controllo: ricalcolo complesso {rc0} dall'xtc di 09 ...", flush=True)
diff_ctrl = float(np.abs(secondaria_da_xtc("complesso", rc0) - serie[("complesso", rc0)]["rmsd_epi_proprio"]).max())
if not diff_ctrl < TOLL_CONTROLLO:
    raise SystemExit(f"secondaria: ricalcolo dall'xtc diverso dalla cache di 10 per complesso {rc0} ({diff_ctrl:.2e} A >= {TOLL_CONTROLLO})")
ingressi_sec = [xtc_della_cache(serie[("complesso", rc0)])]
if FIN_LUNG:
    vl_tutte = sorted({r for sc in C.SCARTI for r in C.valide(ctrl, "libero", "scopo2", sc)["valide"]}, key=lambda r: int(r[3:]))
    sec_l = {}
    for r in vl_tutte:                                   # niente cache su disco della secondaria del libero
        print(f"   secondaria, libero {r} dall'xtc di 09 ...", flush=True)
        sec_l[r] = secondaria_da_xtc("libero", r); ingressi_sec.append(xtc_della_cache(serie[("libero", r)]))
    sec, fin_s = analizza(lambda r: serie[("complesso", r)]["rmsd_epi_proprio"], lambda r: sec_l[r])
else:
    sec, fin_s = analizza(lambda r: serie[("complesso", r)]["rmsd_epi_proprio"], None)

P = C.prov(__file__, [v["_file"] for v in serie.values()] + [os.path.join(arg.out, "controlli.json"), arg.pdb] + ingressi_sec,
           {"cache": arg.cache, "out": arg.out, "trj": arg.trj, "pdb": arg.pdb, "controllo": rc0},
           {"CACHE": arg.cache, "OUT": arg.out, "TRJ": arg.trj})
L = C.testo_prov(P) + C.avviso_intermedia(stato)
L += ["", f"Scopo 2 ({stato}): soglia dal complesso fino a T3 = {T3:.1f} ns; libero fino a T2 = {T2:.1f} ns."]


def righe(ris, titolo, con_libero):
    out = ["", f"##### {titolo}"]
    for sc in C.SCARTI:
        r = ris[sc]; s = r["soglie"]
        out += ["", f"=== Scarto {sc:.0f} ns: complesso valide {r['complesso_valide']}"
                + (f", libero valide {r['libero_valide']}" if con_libero else "")
                + (f" (soglia frazioni {r['k']}, cicli {r['k_cicli']})" if FIN and con_libero else ""),
                f"Soglia P95 = {s['P95']:.3f} A (P90 {s['P90']:.3f}, P99 {s['P99']:.3f}; P95 per replica {s['P95_replica_min']:.3f}-{s['P95_replica_max']:.3f})"
                + ("  ATTENZIONE: soglia > 2 A" if s["P95"] > C.SOGLIA_EPI_FLAG else ""),
                "RMSD media nel complesso: " + ", ".join(f"{x} {v:.2f}" for x, v in r["complesso_media"].items())]
        if not con_libero:
            continue
        out.append(f"{'libero':8s}{'RMSD med':>9s}{'fraz P95':>9s}{'usc':>5s}{'rientri':>8s}{'ingr':>5s}{'iniz':>7s}{'fraz P90':>9s}{'fraz P99':>9s}   blocchi P95 (da {C.BLOCCO_SCOPO2_NS:.0f} ns)")
        for x, v in r["libero"].items():
            p = v["P95"]
            out.append(f"{x:8s}{v['rmsd_media']:9.2f}{p['frazione']:9.2f}{p['uscite']:5d}{p['rientri']:8d}{p['ingressi_senza_uscita']:5d}"
                       f"{(p['stato_iniziale'] or '-'):>7s}{v['P90']['frazione']:9.2f}{v['P99']['frazione']:9.2f}   "
                       + " ".join(f"{b[2]:.2f}" if b[2] is not None else "-" for b in p["blocchi"]))
        if FIN:
            out.append("Esiti a questo scarto: " + "; ".join(f"{k} {v}" for k, v in r["esiti"].items()))
    return out


L += righe(prim, "PRIMARIA (doc): RMSD del backbone dell'epitopo dopo fit sui 43 Ca del core", True)
L += righe(sec, f"SECONDARIA ({C.SECONDARIA_ETICHETTA}): RMSD del backbone dell'epitopo dopo fit sull'epitopo stesso", FIN_LUNG)
L += ["", f"Controllo della secondaria: complesso {rc0} ricalcolato dall'xtc di 09 con comune_fab.rmsd_epitopo_proprio = cache di 10, "
      f"scarto massimo {diff_ctrl:.1e} A (tolleranza {TOLL_CONTROLLO:.0e} A)."]
if not FIN_LUNG:
    L.append("Secondaria nel libero: non calcolata (si sblocca solo con FINALE per lunghezze, condizione del revisore del 6 nov 2026).")
L.append("")
if FIN:
    L.append("Esito PRIMARIA (stesso a 50 e 100 ns): " + "; ".join(f"{k}: {v}" for k, v in fin_p.items()))
    L.append(f"Esito SECONDARIA ({C.SECONDARIA_ETICHETTA}): " + "; ".join(f"{k}: {v}" for k, v in fin_s.items()))
    conc = {k: fin_p[k] == fin_s[k] for k in fin_p}
    L.append("Primaria e secondaria: " + "; ".join(f"{k} esiti {'concordi' if v else 'discordi'}" for k, v in conc.items())
             + ". Se discordi si riportano entrambi, senza scegliere.")
else:
    conc = None
    L.append("Nessun esito: analisi intermedia (solo numeri per replica).")
L += ["", "Da riportare con i risultati (documento dei criteri): l'acqua interna E206 sta nel cuore dell'epitopo",
      "(N di Asp86 e Cys87, O di Asp168); OD1 di Asp40 lega l'N di Ala90; l'epitopo circonda il ponte Cys43-Cys87 e",
      "le MD sono solo ossidate. Formulazioni: mai 'selezione conformazionale' ne' 'induced fit'."]
open(os.path.join(arg.out, "scopo2_epitopo.txt"), "w").write("\n".join(L) + "\n")
json.dump({"provenienza": P, "stato": stato, "finale_per_lunghezze": FIN_LUNG, "T3_ns": T3, "T2_ns": T2,
           "per_scarto": {str(int(k)): v for k, v in prim.items()}, "esito": fin_p,
           "secondaria": {"etichetta": C.SECONDARIA_ETICHETTA, "per_scarto": {str(int(k)): v for k, v in sec.items()},
                          "esito": fin_s, "controllo": {"replica": rc0, "diff_max_A": diff_ctrl, "tolleranza_A": TOLL_CONTROLLO},
                          "libero_calcolato": FIN_LUNG},
           "esiti_concordi": conc}, open(os.path.join(arg.out, "scopo2_epitopo.json"), "w"), indent=1)
print("\n".join(L))
