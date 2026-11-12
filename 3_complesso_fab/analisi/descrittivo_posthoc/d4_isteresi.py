"""POST HOC, 12 nov 2026 (dopo i risultati finali): descrizione delle uscite dall'involucro e sensibilita' all'isteresi.
NON e' un criterio: l'isteresi registrata e' 1 ns e l'esito registrato resta quello di 13_scopo2_epitopo.py.
Legge le serie per fotogramma d1_libero_run<r>.npz (prim = fit sul core, sec = fit sull'epitopo) e le soglie dal json di 13."""
import json, sys, numpy as np
W = sys.argv[1]
J = json.load(open(f"{W}/risultati/scopo2_epitopo.json"))

def transizioni(dentro, nmin):
    stato = None; uscite = rientri = ingressi = 0; i = 0; n = len(dentro); ev = []
    while i < n:
        j = i
        while j < n and dentro[j] == dentro[i]:
            j += 1
        v = bool(dentro[i])
        if j - i >= nmin:
            if stato is None:
                stato = v
            elif v != stato:
                if not v: uscite += 1
                elif uscite > rientri: rientri += 1
                else: ingressi += 1
                stato = v
                ev.append((i, v))
        i = j
    return uscite, rientri, ev

out = []
for nome, chiave, ramo in [("PRIMARIA (fit sul core)", "prim", J), ("SECONDARIA (fit sull'epitopo)", "sec", J["secondaria"])]:
    for sc in ("50", "100"):
        P = ramo["per_scarto"][sc]; s95 = P["soglie"]["P95"]
        out.append(f"\n=== {nome}, scarto {sc} ns, soglia P95 {s95:.3f} A")
        for ist in (1.0, 2.0, 5.0):
            ge = 0; cy = 0; righe = []
            for r in range(1, 11):
                d = np.load(f"{W}/descrittivo/d1_libero_run{r}.npz")
                t = d["t_ns"]; w = (t >= float(sc) - 1e-6) & (t <= 1000.0 + 1e-6)
                x = d[chiave][w]; tt = t[w]; dentro = x <= s95
                nmin = int(round(ist * 1000.0 / 100.0))
                u, ri, ev = transizioni(dentro, nmin)
                fr = float(dentro.mean()); ge += fr >= 0.5; cy += ri >= 1
                if ist == 1.0:
                    rif = P["libero"][f"run{r}"]["P95"]
                    assert abs(rif["frazione"] - fr) < 1e-12 and rif["uscite"] == u and rif["rientri"] == ri, (r, rif["frazione"], fr, rif["uscite"], u, rif["rientri"], ri)
                    # durata dei tratti fuori (stato con isteresi): da ogni uscita al rientro successivo
                    dur = []; t_u = None
                    for i, v in ev:
                        if not v: t_u = tt[i]
                        elif t_u is not None: dur.append((t_u, tt[i] - t_u)); t_u = None
                    coda = f"; ultima uscita a {t_u:.1f} ns senza rientro" if t_u is not None else ""
                    righe.append(f"  run{r:<2d} frazione {fr:.3f}, uscite {u}, rientri {ri}; tratti fuori (inizio ns, durata ns): "
                                 + (", ".join(f"({a:.1f}, {b:.1f})" for a, b in dur) if dur else "nessuno") + coda)
                else:
                    righe.append(f"  run{r:<2d} uscite {u}, rientri {ri}")
            es = ("popolata" if cy >= 5 else "non distinguibile dalla memoria della partenza") if ge >= 8 else "intermedia (o non popolata)"
            out.append(f"-- isteresi {ist:.0f} ns{' (REGISTRATA; controllata uguale al json di 13)' if ist == 1.0 else ' (post hoc)'}: "
                       f"repliche con frazione >= 0.5: {ge}/10; con almeno un ciclo: {cy}/10; conteggio secondo la regola: {es}")
            out += righe
print("\n".join(out))
