"""Controlli tecnici e validita' delle repliche (criteri, sezione "Controlli tecnici e validita' delle repliche").

Uso (dalla cartella complesso_fab):
    python3 analisi/11_controlli.py [--cache analisi/cache] [--out analisi/risultati] [--pdb mrka-fab_completo.pdb] [--finale]

Finestre (tempo dalla produzione): dallo scarto (50 ns; 100 ns di sensibilita') alla durata confermativa.
T3 = minimo fra le repliche del complesso (max 500 ns): finestra del complesso e del libero per lo scopo 3.
T2 = minimo fra le repliche del libero (max 1 us): finestra del libero per lo scopo 2.
Replica non valida (criteri): RMSD dei 43 Ca del core > 2.0 A in > 50% dei fotogrammi; STOP tecnico (controllo_T.txt);
distanza dalle immagini periodiche < 12 A in > 5% dei fotogrammi (sotto 12 A in un fotogramma: solo segnalata);
solo complesso: persi > 50% dei contatti MrkA-Fab del cristallo (4.0 A); distanza dei centri cresciuta > 5 A.
Il documento non dice su che base si contano i contatti persi e la crescita dei centri (NON NEL DOC): si calcolano
due letture per ciascuno e, se danno esiti diversi, la replica e' marcata "da decidere" (decide il revisore).
  contatti persi, lettura 1: frazione dei contatti del cristallo con occupanza < 0.5 nella finestra (come "perso" nello scopo 1)
  contatti persi, lettura 2: frazione dei fotogrammi in cui sono presenti meno della meta' dei contatti del cristallo, > 50%
  centri, lettura 1: crescita > 5 A in > 50% dei fotogrammi (come la regola dell'RMSD)
  centri, lettura 2: crescita > 5 A in almeno un fotogramma
Centri = centri geometrici dei Ca della catena A e dei Ca del Fab (L + H); crescita rispetto al cristallo completato.
Si registrano anche (non decidono la validita'): RMSD dei Ca del Fab e dell'Fv, proline cis, temperature per tratto.
Scrive OUT/controlli.txt e OUT/controlli.json (letto dagli script degli scopi).
"""
import os, sys, json, argparse
import numpy as np
import comune_fab as C

ap = argparse.ArgumentParser()
ap.add_argument("--cache", default="analisi/cache")
ap.add_argument("--out", default="analisi/risultati")
ap.add_argument("--pdb", default="mrka-fab_completo.pdb")
ap.add_argument("--radice", default=".", help="cartella complesso_fab (per <sistema>/run<r>/controllo_T.txt)")
ap.add_argument("--finale", action="store_true", help="dichiara finale un disegno ridotto (decisione di francesco)")
arg = ap.parse_args()
os.makedirs(arg.out, exist_ok=True)

serie = C.leggi_cache(arg.cache)
T3, T2, lung = C.durate(serie)
stato = C.stato_analisi(lung, T3, T2, arg.finale)
cristallo = sorted(C.contatti_cristallo(arg.pdb, C.CONTATTO_PRINCIPALE))
icr = {p: i for i, p in enumerate(cristallo)}


def temperature(s, r):
    f = os.path.join(arg.radice, s, r, "controllo_T.txt")
    if not os.path.exists(f):
        return None, "controllo_T.txt non trovato"
    righe = [l.split() for l in open(f) if l.strip()]
    stop = [l for l in righe if len(l) >= 5 and l[4].startswith("STOP")]
    avv = [l for l in righe if len(l) >= 5 and l[4].startswith("AVVISO")]
    txt = "; ".join(f"{l[0]}-{l[1]} ns: prot {l[2]} K, acqua {l[3]} K, {l[4]}" for l in righe)
    return (len(stop) > 0, len(avv)), txt


def misure(d, s, w, temp):
    m, motivi, dubbi = {}, [], []
    n = int(w.sum()); m["fotogrammi"] = n
    rc = d["rmsd_core"][w]
    m["rmsd_core_media"] = float(rc.mean()); m["rmsd_core_max"] = float(rc.max())
    m["rmsd_core_frazione_sopra_2A"] = float((rc > C.RMSD_CORE_MAX).mean())
    if m["rmsd_core_frazione_sopra_2A"] > C.RMSD_CORE_FRAZ:
        motivi.append("RMSD core > 2 A in > 50% dei fotogrammi")
    pi = d["mindist_pi"][w]
    m["immagini_min_A"] = float(np.nanmin(pi)) if np.any(~np.isnan(pi)) else None     # None = oltre 25 A
    m["immagini_frazione_sotto_12A"] = float((np.nan_to_num(pi, nan=1e9) < C.PI_MIN).mean())
    m["immagini_segnalata"] = m["immagini_frazione_sotto_12A"] > 0
    if m["immagini_frazione_sotto_12A"] > C.PI_FRAZ:
        motivi.append("immagini periodiche < 12 A in > 5% dei fotogrammi")
    if temp[0] is None:
        m["temperatura"] = "non letta"
    else:
        m["temperatura"] = {"stop": temp[0][0], "avvisi": temp[0][1], "tratti": temp[1]}
        if temp[0][0]:
            motivi.append("STOP tecnico (temperatura)")
    if s == "complesso":
        fr = np.where(w)[0]
        sel = np.isin(d["c_f"], fr) & (d["c_d"] <= C.CONTATTO_PRINCIPALE)
        pres = np.zeros((len(fr), len(cristallo)), bool)
        pos = {f: i for i, f in enumerate(fr)}
        cc = np.where(d["c_c"] == 0, "L", "H")
        for f, a, c, k in zip(d["c_f"][sel], d["c_a"][sel], cc[sel], d["c_k"][sel]):
            j = icr.get((int(a), c, int(k)))
            if j is not None:
                pres[pos[f], j] = True
        occ = pres.mean(0); Q = pres.mean(1)
        m["contatti_cristallo"] = len(cristallo)
        m["Q_media"] = float(Q.mean())
        m["contatti_persi_lettura1"] = float((occ < 0.5).mean())
        m["fotogrammi_Q_sotto_0.5"] = float((Q < 0.5).mean())
        l1 = m["contatti_persi_lettura1"] > C.CONTATTI_PERSI_FRAZ
        l2 = m["fotogrammi_Q_sotto_0.5"] > 0.5
        if l1 and l2:
            motivi.append("persi > 50% dei contatti del cristallo (entrambe le letture)")
        elif l1 != l2:
            dubbi.append("contatti persi: le due letture discordano")
        dc = d["dcentri"][w] - float(d["dcentri_cristallo"])
        m["centri_crescita_media"] = float(dc.mean()); m["centri_crescita_max"] = float(dc.max())
        m["centri_frazione_sopra_5A"] = float((dc > C.CENTRI_MAX).mean())
        l1 = m["centri_frazione_sopra_5A"] > 0.5
        l2 = m["centri_crescita_max"] > C.CENTRI_MAX
        if l1 and l2:
            motivi.append("centri cresciuti > 5 A (entrambe le letture)")
        elif l1 != l2:
            dubbi.append("centri: le due letture discordano")
        m["rmsd_fab_media"] = float(d["rmsd_fab"][w].mean()); m["rmsd_fv_media"] = float(d["rmsd_fv"][w].mean())
        om = d["omega"][w]
        m["pro_cis_frazione"] = {str(p): float((np.abs(om[:, i]) < 90).mean()) for i, p in enumerate(d["pro_cis"])}
    m["motivi"] = motivi; m["dubbi"] = dubbi
    m["esito"] = "non valida" if motivi else ("da decidere" if dubbi else "valida")
    return m


ris, righe = {}, []
for (s, r), d in sorted(serie.items(), key=lambda x: (x[0][0], int(x[0][1][3:]))):
    temp = temperature(s, r)
    ris[f"{s}_{r}"] = {}
    for sc in C.SCARTI:
        fin = {"complesso": {"complesso": T3}, "libero": {"scopo3": T3, "scopo2": T2}}[s]
        ris[f"{s}_{r}"][str(int(sc))] = {k: misure(d, s, C.finestra(d["t_ns"], sc, T), temp) for k, T in fin.items()}

riep = {}
for s, finestre in (("complesso", ["complesso"]), ("libero", ["scopo3", "scopo2"])):
    for k in finestre:
        for sc in C.SCARTI:
            reps = [x for x in ris if x.startswith(s + "_")]
            val = [x.split("_")[1] for x in reps if ris[x][str(int(sc))][k]["esito"] != "non valida"]
            dub = [x.split("_")[1] for x in reps if ris[x][str(int(sc))][k]["esito"] == "da decidere"]
            n = len(val)
            riep[f"{s}|{k}|{int(sc)}"] = {"presenti": len(reps), "valide": C.ordina_repliche(val), "da_decidere": dub,
                                          "soglia_repliche": C.soglia_repliche(n) if n else None,
                                          "concludibile": n >= C.MIN_VALIDE[s]}

P = C.prov(__file__, [v["_file"] for v in serie.values()] + [arg.pdb],
           {"cache": arg.cache, "out": arg.out, "pdb": arg.pdb, "finale": arg.finale}, {"CACHE": arg.cache, "OUT": arg.out})
json.dump({"provenienza": P, "stato": stato, "T3_ns": T3, "T2_ns": T2, "lunghezze_ns": lung,
           "contatti_cristallo_4A": [list(p) for p in cristallo], "repliche": ris, "riepilogo": riep},
          open(os.path.join(arg.out, "controlli.json"), "w"), indent=1)

L = C.testo_prov(P) + C.avviso_intermedia("FINALE" if all(v == "FINALE" for v in stato.values()) else "INTERMEDIA")
L += ["", f"Stato: scopo 1 {stato['scopo1']}, scopo 2 {stato['scopo2']}, scopo 3 {stato['scopo3']}",
      f"Lunghezze (ns di produzione): " + "; ".join(f"{s} " + ", ".join(f"{r} {v:.1f}" for r, v in sorted(x.items(), key=lambda y: int(y[0][3:]))) for s, x in lung.items()),
      f"Durate confermative: T3 = {T3:.1f} ns (complesso; libero per lo scopo 3), T2 = {T2:.1f} ns (libero, scopo 2)",
      f"Contatti MrkA-Fab del cristallo completato a 4.0 A: {len(cristallo)} coppie di residui", ""]
for sc in C.SCARTI:
    L.append(f"=== Scarto {sc:.0f} ns")
    L.append(f"{'replica':16s}{'finestra':10s}{'fotog':>6s}{'RMSDcore':>9s}{'>2A':>6s}{'immag.min':>10s}{'<12A':>6s}"
             f"{'Q':>6s}{'persi1':>7s}{'Q<.5':>6s}{'dCentri':>8s}{'max':>6s}{'>5A':>6s}{'Fab':>6s}{'Fv':>6s}  esito")
    for x, v in ris.items():
        for k, m in v[str(int(sc))].items():
            im = f"{m['immagini_min_A']:.1f}" if m["immagini_min_A"] is not None else ">25"
            riga = f"{x:16s}{k:10s}{m['fotogrammi']:6d}{m['rmsd_core_media']:9.2f}{m['rmsd_core_frazione_sopra_2A']:6.2f}{im:>10s}{m['immagini_frazione_sotto_12A']:6.3f}"
            if "Q_media" in m:
                riga += (f"{m['Q_media']:6.2f}{m['contatti_persi_lettura1']:7.2f}{m['fotogrammi_Q_sotto_0.5']:6.2f}"
                         f"{m['centri_crescita_media']:8.2f}{m['centri_crescita_max']:6.2f}{m['centri_frazione_sopra_5A']:6.2f}"
                         f"{m['rmsd_fab_media']:6.2f}{m['rmsd_fv_media']:6.2f}")
            else:
                riga += " " * 64
            riga += f"  {m['esito']}" + (f" ({'; '.join(m['motivi'] + m['dubbi'])})" if m["motivi"] or m["dubbi"] else "")
            if m["immagini_segnalata"]:
                riga += " [immagini < 12 A in almeno un fotogramma]"
            L.append(riga)
    L.append("")
L.append("Proline cis del cristallo (frazione di fotogrammi con |omega| < 90, scarto 50 ns):")
for x, v in ris.items():
    m = v["50"].get("complesso")
    if m:
        L.append(f"  {x}: " + ", ".join(f"{p} {f:.3f}" for p, f in m["pro_cis_frazione"].items()))
L.append("Temperature per tratto (controllo_T.txt):")
for x, v in ris.items():
    t = list(v["50"].values())[0]["temperatura"]
    L.append(f"  {x}: " + (t if isinstance(t, str) else f"{'STOP' if t['stop'] else 'nessuno STOP'}, {t['avvisi']} avvisi; {t['tratti']}"))
L.append("")
L.append("Riepilogo (soglia = ceil(0.8 x valide); non concludibile con meno di 4 valide nel complesso o 8 nel libero):")
for k, v in riep.items():
    s, f, sc = k.split("|")
    L.append(f"  {s:10s} {f:10s} scarto {sc:>3s}: presenti {v['presenti']}, valide {len(v['valide'])} {v['valide']}"
             + (f", da decidere {v['da_decidere']}" if v["da_decidere"] else "")
             + f", soglia {v['soglia_repliche']}, {'concludibile' if v['concludibile'] else 'NON CONCLUDIBILE'}")
open(os.path.join(arg.out, "controlli.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
