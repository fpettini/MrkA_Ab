"""Copia di analisi_rest2/provenienza.py (cambia solo VERSIONE) per 05_rmsf_vs_bfactor.py.
Provenienza dei risultati di analisi_rest2: chi ha prodotto un file, da quali input e con quali pacchetti.

Ogni script chiama prov(__file__, input, argomenti) e mette il dizionario sotto "provenienza" nel json e
le righe di testo() in testa al txt. Lo sha256 dello script si calcola a runtime, quindi non dipende da un
numero aggiornato a mano; VERSIONE e' solo un'etichetta leggibile, da aggiornare a ogni modifica.
"""
import os, sys, json, hashlib, platform, datetime
import numpy as np
from importlib import metadata

VERSIONE = "bfactor 2026-10-29"
PACCHETTI = ["numpy", "scipy", "MDAnalysis", "pymbar"]

def sha256(path, blocco=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(blocco), b""):
            h.update(b)
    return h.hexdigest()

def etichetta(path, basi):
    """Percorso relativo a una cartella base, scritto <BASE>/resto; cosi' non compaiono percorsi della macchina."""
    a = os.path.abspath(path)
    for nome, d in (basi or {}).items():
        d = os.path.abspath(d)
        if a == d:
            return f"<{nome}>"
        if a.startswith(d + os.sep):
            return f"<{nome}>/" + os.path.relpath(a, d).replace(os.sep, "/")
    return os.path.basename(a)

def prov(script, inputs, argomenti, basi=None):
    """basi: {"ROOT": cartella del progetto, "CACHE": ..., "OUT": ...} per i percorsi relativi."""
    pk = {}
    for p in PACCHETTI:
        try:
            pk[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            pk[p] = None
    return {"versione": VERSIONE, "script": os.path.basename(script),
            "sha256_script": sha256(script),
            "sha256_provenienza_py": sha256(os.path.abspath(__file__)),
            "data": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
            "python": platform.python_version(), "pacchetti": pk,
            "argomenti": {k: (etichetta(v, basi) if isinstance(v, str) and os.path.exists(v) else v) for k, v in argomenti.items()},
            "input": {etichetta(i, basi): sha256(i) for i in inputs if os.path.exists(i)}}

def testo(p):
    r = [f"# {p['script']} ({p['versione']}), sha256 {p['sha256_script'][:16]}, provenienza.py {p['sha256_provenienza_py'][:16]}",
         f"# eseguito {p['data']}, Python {p['python']}, " + ", ".join(f"{k} {v}" for k, v in p["pacchetti"].items()),
         f"# argomenti: {p['argomenti']}"]
    r += [f"# input {k}: sha256 {v[:16]}" for k, v in p["input"].items()]
    for c, q in (p.get("cache") or {}).items():
        if q is None:
            r.append(f"# cache {c}: senza provenienza (fatta prima del 25 ott)")
        else:
            r += [f"# cache {c} <- {k}: sha256 {v[:16]} ({q['script']} {q['sha256_script'][:16]})" for k, v in q["input"].items()]
    return r

def provenienza_cache(cache, sistemi, repliche):
    """Provenienza scritta da 01_estrai.py dentro ogni cache (None per le cache vecchie, senza il campo)."""
    o = {}
    for s in sistemi:
        for r in repliche:
            f = os.path.join(cache, f"{s}_{r}.npz")
            if os.path.exists(f):
                d = np.load(f, allow_pickle=True)
                o[f"{s}_{r}"] = json.loads(str(d["provenienza"])) if "provenienza" in d.files else None
    return o
