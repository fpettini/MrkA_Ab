"""Fase 3: figure di convergenza dai file prodotti da 02_analisi.py.

Uso (dalla cartella del progetto mrkA):
    python analisi_convergenza/03_figure.py analisi_convergenza/cache analisi_convergenza
"""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
C=sys.argv[1] if len(sys.argv)>1 else "analisi_convergenza/cache"
O=sys.argv[2] if len(sys.argv)>2 else "analisi_convergenza"
SIST=["xray","nmr_ox","nmr_red"]; TIT={"xray":"X-ray ossidato","nmr_ox":"NMR ossidato","nmr_red":"NMR ridotto"}
COL=["#2a78d6","#eb6834","#1baf7a"]          # replica 1,2,3 (palette validata, all-pairs)
SURF="#fcfcfb"; INK="#0b0b0b"; INK2="#52514e"; GRID="#e4e3df"; SHADE="#ecebe7"
plt.rcParams.update({"figure.facecolor":SURF,"axes.facecolor":SURF,"savefig.facecolor":SURF,
    "axes.edgecolor":INK2,"axes.labelcolor":INK,"xtick.color":INK2,"ytick.color":INK2,"text.color":INK,
    "axes.grid":True,"grid.color":GRID,"grid.linewidth":0.6,"axes.spines.top":False,"axes.spines.right":False,
    "font.size":9,"axes.titlesize":10,"axes.titleweight":"bold","legend.frameon":False})
def media_mobile(y,w=50): return np.convolve(y,np.ones(w)/w,mode="valid")

# ---------- Fig 1: serie temporali ----------
fig,ax=plt.subplots(3,3,figsize=(11,8),sharex=True,sharey="row")
fine={}; TMAX=0
righe=[("rmsd","RMSD Cα 44-199 (Å)"),("rb","RMSD core β (Å)"),("rg","Rg atomi pesanti (Å)")]
for j,s in enumerate(SIST):
    for r in (1,2,3):
        a=np.load(f"{C}/{s}_run{r}_serie.npz")
        t=a["t"]; TMAX=max(TMAX,t[-1]); ser={"rmsd":a["rmsd"],"rb":a["rb"],"rg":a["rg"]}
        for i,(k,lab) in enumerate(righe):
            y=ser[k]; ax[i,j].plot(t,y,color=COL[r-1],lw=0.6,alpha=0.22)
            ym=media_mobile(y); tm=media_mobile(t)
            ax[i,j].plot(tm,ym,color=COL[r-1],lw=2,label=f"replica {r}")
            fine.setdefault((i,j),[]).append((ym[-1],r))
    ax[0,j].set_title(TIT[s])
for (i,j),v in fine.items():
    v.sort(); lo,hi=ax[i,j].get_ylim(); gap=(hi-lo)*0.06; pos=[]
    for y,r in v:
        y=max(y,pos[-1]+gap) if pos else y; pos.append(y)
        ax[i,j].annotate(f"r{r}",(TMAX,y),xytext=(5,0),textcoords="offset points",va="center",fontsize=8,color=INK2,annotation_clip=False)
for i,(k,lab) in enumerate(righe): ax[i,0].set_ylabel(lab)
for j in range(3): ax[2,j].set_xlabel("tempo (ns)"); ax[2,j].set_xlim(0,TMAX*1.06)
h,l=ax[0,0].get_legend_handles_labels(); fig.legend(h,l,loc="upper center",ncol=3,bbox_to_anchor=(0.5,1.0))
fig.text(0.5,0.955,"Linee sottili: valore per frame (100 ps). Linee spesse: media mobile su 5 ns. Riferimento: struttura di partenza.",ha="center",fontsize=8,color=INK2)
fig.tight_layout(rect=(0,0,1,0.94)); fig.savefig(f"{O}/fig1_serie_temporali.png",dpi=200); plt.close(fig)

# ---------- Fig 2: RMSF per replica ----------
RES=np.arange(42,203)
fig,ax=plt.subplots(3,1,figsize=(11,7.5),sharex=True,sharey=True)
for j,s in enumerate(SIST):
    c=np.load(f"{C}/{s}_confronto.npz"); rm=c["rmsf"]
    for lo,hi,txt in ((41.5,43.5,"esclusi"),(199.5,202.5,"esclusi"),(115.5,117.5,"116-117\nricostruiti")):
        ax[j].axvspan(lo,hi,color=SHADE,zorder=0,lw=0)
    for r in range(3): ax[j].plot(RES,rm[r],color=COL[r],lw=1.6,label=f"replica {r+1}")
    sd=rm.std(0,ddof=1)
    ax[j].set_title(TIT[s],loc="left"); ax[j].set_ylabel("RMSF Cα (Å)")
ax[0].text(116.5,ax[0].get_ylim()[1]*0.97,"116-117 ricostruiti",ha="center",va="top",fontsize=7,color=INK2)
ax[0].text(42.5,ax[0].get_ylim()[1]*0.97,"esclusi",ha="left",va="top",fontsize=7,color=INK2)
ax[2].set_xlabel("residuo (numerazione UniProt P12267)"); ax[2].set_xlim(41,203)
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc="upper center",ncol=3,bbox_to_anchor=(0.5,1.0))
fig.text(0.5,0.955,f"Sovrapposizione sui Cα del core β, {TMAX:.0f} ns per replica. Bande grigie: residui esclusi dalle conclusioni per residuo.",ha="center",fontsize=8,color=INK2)
fig.tight_layout(rect=(0,0,1,0.94)); fig.savefig(f"{O}/fig2_rmsf_repliche.png",dpi=200); plt.close(fig)

# ---------- Fig 3: RMSD incrociato ----------
blu=LinearSegmentedColormap.from_list("blu",["#cde2fb","#86b6ef","#3987e5","#1c5cab","#0d366b"])
fig,ax=plt.subplots(1,3,figsize=(12,4.4))
vmax=5.0
for j,s in enumerate(SIST):
    D=np.load(f"{C}/{s}_confronto.npz")["D"]
    im=ax[j].imshow(D,cmap=blu,vmin=1.0,vmax=vmax,origin="lower",interpolation="nearest")
    n3=len(D)//3
    for b in (n3,2*n3): ax[j].axhline(b-0.5,color=SURF,lw=2); ax[j].axvline(b-0.5,color=SURF,lw=2)
    ax[j].set_xticks([n3/2,1.5*n3,2.5*n3]); ax[j].set_xticklabels(["r1","r2","r3"])
    ax[j].set_yticks([n3/2,1.5*n3,2.5*n3]); ax[j].set_yticklabels(["r1","r2","r3"])
    ax[j].grid(False); ax[j].set_title(TIT[s])
cb=fig.colorbar(im,ax=ax,fraction=0.025,pad=0.02); cb.set_label("RMSD Cα 44-199 (Å), scala 1-5 Å"); cb.outline.set_visible(False)
fig.suptitle(f"RMSD fra coppie di frame (1 frame/ns, {len(D)//3} per replica). Blocchi diagonali: dentro la replica; fuori diagonale: fra repliche.",fontsize=8,color=INK2,fontweight="normal",y=0.98)
fig.savefig(f"{O}/fig3_rmsd_incrociato.png",dpi=200,bbox_inches="tight"); plt.close(fig)

# ---------- Fig 4: distanza SG-SG nel ridotto ----------
fig,ax=plt.subplots(figsize=(7,3.6))
bins=np.arange(2.8,8.6,0.1)
for r in (1,2,3):
    ss=np.load(f"{C}/nmr_red_run{r}_serie.npz")["ss"]
    ax.hist(ss,bins=bins,histtype="step",lw=2,color=COL[r-1],density=True,label=f"replica {r}: {np.mean(ss<3.5)*100:.0f}% dei frame < 3.5 Å")
ax.axvline(3.5,color=INK2,lw=1,ls="--"); ax.text(3.45,ax.get_ylim()[1]*0.95,"3.5 Å",ha="right",va="top",fontsize=8,color=INK2)
ax.set_xlabel("distanza SG43-SG87 (Å)"); ax.set_ylabel("densità"); ax.set_title("NMR ridotto: distanza fra gli zolfi",loc="left")
ax.legend(loc="upper right",fontsize=8)
fig.tight_layout(); fig.savefig(f"{O}/fig4_ss_ridotto.png",dpi=200); plt.close(fig)

# ---------- Fig 5: convergenza in funzione della durata ----------
COLS={"xray":"#2a78d6","nmr_ox":"#eb6834","nmr_red":"#1baf7a"}   # qui il colore e' il sistema, non la replica
BREVE={"xray":"X-ray","nmr_ox":"NMR ox","nmr_red":"NMR red"}
fig,ax=plt.subplots(1,3,figsize=(12,4.2))
voci=[(1,"RMSIP medio fra repliche (10 PC)","più alto = stesse direzioni di moto"),
      (2,"RMSD fra repliche / dentro replica","più vicino a 1 = stesse conformazioni"),
      (3,"r medio dell'RMSF fra repliche","più alto = stesso profilo di flessibilità")]
for j,(c,tit,sot) in enumerate(voci):
    fin=[]
    for s in SIST:
        e=np.load(f"{C}/{s}_confronto.npz")["evo"]
        ax[j].plot(e[:,0],e[:,c],color=COLS[s],lw=2,marker="o",ms=4,label=TIT[s])
        fin.append((e[-1,c],s,e[-1,0]))
    fin.sort(); lo,hi=ax[j].get_ylim(); gap=(hi-lo)*0.07; prev=None
    for y,s,x in fin:
        yy=y if prev is None else max(y,prev+gap); prev=yy
        ax[j].annotate(BREVE[s],(x,yy),xytext=(6,0),textcoords="offset points",va="center",fontsize=7.5,color=INK2,annotation_clip=False)
    ax[j].set_title(f"{tit}\n",loc="left",fontsize=9)
    ax[j].text(0,1.02,sot,transform=ax[j].transAxes,fontsize=7.5,color=INK2,va="bottom")
    ax[j].set_xlabel("durata analizzata 0-T (ns)"); ax[j].set_xlim(80,e[-1,0]*1.2)
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc="upper center",ncol=3,bbox_to_anchor=(0.5,1.0))
fig.tight_layout(rect=(0,0,1,0.9)); fig.savefig(f"{O}/fig5_convergenza_vs_durata.png",dpi=200); plt.close(fig)
print("ok")
