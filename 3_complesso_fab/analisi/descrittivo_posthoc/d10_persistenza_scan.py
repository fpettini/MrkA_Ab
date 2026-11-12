"""POST HOC, 12 nov 2026 (dopo i risultati finali): numero di repliche del libero con almeno un ciclo uscita-rientro in funzione
della permanenza richiesta (la regola registrata chiede 1 ns = 10 fotogrammi di fila); tratti di fila fuori dall'involucro;
ultimo rientro; cicli nella sola finestra 500-1000 ns. NON e' un criterio: l'esito registrato resta quello di 13_scopo2_epitopo.py.
Nasce dalla lettura critica del 12 nov, che ha mostrato che i due punti (2 e 5 ns) di d4 non dicono dove cambia la categoria.
Legge le serie per fotogramma d1_libero_run<r>.npz (prim = fit sul core, sec = fit sull'epitopo) e le soglie dal json di 13.
Uso: python3 d10_persistenza_scan.py <cartella con descrittivo/ e risultati/>"""
import json, sys, numpy as np
W = sys.argv[1]
J = json.load(open(f"{W}/risultati/scopo2_epitopo.json"))
D = {r: np.load(f"{W}/descrittivo/d1_libero_run{r}.npz") for r in range(1, 11)}

def transizioni(dentro, nmin):
    """Stessa macchina a stati di 13_scopo2_epitopo.py e di d4: stato iniziale = primo tratto di almeno nmin fotogrammi;
    uscita = tratto fuori di almeno nmin dopo essere stato dentro; rientro = tratto dentro di almeno nmin dopo un'uscita."""
    stato = None; uscite = rientri = 0; i = 0; n = len(dentro); ev = []
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
                stato = v
                ev.append((i, v, j - i))
        i = j
    return uscite, rientri, ev

def tratti(dentro):
    """Tratti di fila (inizio, lunghezza, dentro?)."""
    out = []; i = 0; n = len(dentro)
    while i < n:
        j = i
        while j < n and dentro[j] == dentro[i]:
            j += 1
        out.append((i, j - i, bool(dentro[i]))); i = j
    return out

P = []; JS = {}
for nome, chiave, ramo in [("PRIMARIA (fit sul core)", "prim", J), ("SECONDARIA (fit sull'epitopo)", "sec", J["secondaria"])]:
    for sc in ("50", "100"):
        R = ramo["per_scarto"][sc]; s95 = R["soglie"]["P95"]
        P.append(f"\n=== {nome}, scarto {sc} ns, involucro P95 {s95:.3f} A")
        serie = {}
        for r in range(1, 11):
            t = D[r]["t_ns"]; w = (t >= float(sc) - 1e-6) & (t <= 1000.0 + 1e-6)
            serie[r] = (t[w], D[r][chiave][w] <= s95)
        # controllo con il json a 1 ns
        for r in range(1, 11):
            u, ri, _ = transizioni(serie[r][1], 10); rif = R["libero"][f"run{r}"]["P95"]
            assert rif["uscite"] == u and rif["rientri"] == ri and abs(rif["frazione"] - serie[r][1].mean()) < 1e-12
        dentro8 = sum(serie[r][1].mean() >= 0.5 for r in range(1, 11))
        P.append(f"repliche con frazione >= 0.5: {dentro8}/10 (non dipende dalla permanenza)")
        P.append("permanenza (ns): repliche con almeno un ciclo [quali]")
        prec = None
        for nmin in list(range(1, 31)) + [35, 40, 50, 60, 70, 80, 100, 150, 200]:
            chi = [r for r in range(1, 11) if transizioni(serie[r][1], nmin)[1] >= 1]
            riga = f"  {nmin/10:5.1f}: {len(chi)}  {chi}"
            if prec != chi or nmin in (10, 20, 50):
                P.append(riga + ("   <- regola registrata" if nmin == 10 else ""))
            prec = chi
        # intervalli di permanenza per ciascun conteggio
        conta = {}; perns = {}
        for nmin in range(1, 201):
            c = sum(transizioni(serie[r][1], nmin)[1] >= 1 for r in range(1, 11))
            conta.setdefault(c, []).append(nmin / 10); perns[f"{nmin / 10:.1f}"] = int(c)
        JS[f"{chiave}_{sc}"] = {"involucro": s95, "dentro": int(dentro8), "repliche_con_cicli": perns, "fuori_di_fila_max_ns": {}, "ultimo_rientro_ns": {}, "tratti_fuori_ge1ns": {}}
        P.append("conteggio -> intervallo di permanenza (ns): " + "; ".join(f"{c}: {min(v):.1f}-{max(v):.1f}" for c, v in sorted(conta.items(), reverse=True)))
        P.append("tratti di fila FUORI dall'involucro di almeno 1 ns (inizio ns, durata ns) e tratto piu' lungo, per replica")
        for r in range(1, 11):
            t, d = serie[r]
            tr = [(t[i], n / 10) for i, n, v in tratti(d) if not v and n >= 10]
            mx = max([n for i, n, v in tratti(d) if not v], default=0)
            u, ri, ev = transizioni(d, 10)
            rient = [t[i] for i, v, n in ev if v]
            JS[f"{chiave}_{sc}"]["fuori_di_fila_max_ns"][f"run{r}"] = mx / 10; JS[f"{chiave}_{sc}"]["ultimo_rientro_ns"][f"run{r}"] = float(rient[-1]) if rient else None
            JS[f"{chiave}_{sc}"]["tratti_fuori_ge1ns"][f"run{r}"] = [(float(a), float(b)) for a, b in tr]
            P.append(f"  run{r:<2d} fuori di fila piu' lungo {mx/10:.1f} ns; tratti >= 1 ns: {len(tr)}"
                     + (": " + ", ".join(f"({a:.1f}, {b:.1f})" for a, b in tr[:6]) + (" ..." if len(tr) > 6 else "") if tr else "")
                     + (f"; ultimo rientro a {rient[-1]:.1f} ns" if rient else ""))
        # sola finestra 500-1000 ns (stesso involucro)
        chi = []
        for r in range(1, 11):
            t, d = serie[r]; w = t >= 500.0 - 1e-6
            if transizioni(d[w], 10)[1] >= 1: chi.append(r)
        P.append(f"solo 500-1000 ns, stesso involucro, regola a 1 ns: cicli in {len(chi)}/10 {chi}")
        JS[f"{chiave}_{sc}"]["cicli_500_1000"] = chi
        for sc2 in (150, 200, 250, 300):
            chi = [r for r in range(1, 11) if transizioni(serie[r][1][serie[r][0] >= sc2 - 1e-6], 10)[1] >= 1]
            P.append(f"finestra {sc2}-1000 ns, stesso involucro: cicli in {len(chi)}/10 {chi}")
json.dump(JS, open(f"{W}/descrittivo/posthoc12/d10_persistenza_scan.json", "w"), indent=1)
print("\n".join(P))
