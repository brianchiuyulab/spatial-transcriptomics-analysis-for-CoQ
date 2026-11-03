"""Render proposal panels at their printed sizes. Input is the analysis release directory."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--release", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
root, out = args.release, args.out
out.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family":"Arial","font.size":7,"axes.linewidth":.55,
    "xtick.major.size":2,"ytick.major.size":2,"xtick.major.width":.5,"ytick.major.width":.5,
    "xtick.major.pad":1.3,"ytick.major.pad":1.4,"axes.labelpad":1.7,
    "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none",
    "axes.unicode_minus":False,"legend.fontsize":7,"legend.frameon":False})
class Panel:
    def __init__(self, label, title, width, height):
        self.label,self.w,self.h=label,width,height
        self.fig=plt.figure(figsize=(width/72,height/72))
        self.text(0,9,label[-1],9,weight="bold")
        self.text(12,9,title,8,weight="bold")
    def text(self,x,y,s,size=7,**kw):
        return self.fig.text(x/self.w,1-y/self.h,s,fontsize=size,**kw)
    def ax(self,x,y,w,h):
        return self.fig.add_axes([x/self.w,1-(y+h)/self.h,w/self.w,h/self.h])
    def save(self):
        for ext in ("pdf","svg"):
            self.fig.savefig(out/f"{self.label}.{ext}",dpi=420,facecolor="white")
        plt.close(self.fig)
    def key(self,x,y,color,label):
        self.text(x,y,"●",7,color=color)
        self.text(x+7,y,label)
def tab(name):
    return pd.read_csv(root/"tables"/name,sep="\t")

import sys
sys.path.insert(0,str(root/"code"))
from spatial_plot_config import CELL_COLORS
from render_proposal_spatial_maps import positive_barcodes
from spatial_plot_config import SAMPLES
meta=pd.read_csv(root/"data/GSE266933_compact_metadata.csv.gz",usecols=["sample","barcode","x_um","y_um","max_pred_celltype","injury_classification"])
span=meta.groupby("sample")[["x_um","y_um"]].agg(lambda x:x.max()-x.min()).to_numpy().max()*1.08
ages=["Young","Geriatric"];styles=["#006fa6","#c4550b"]
for code,width,side in [("E",270,68),("F",202,68)]:
    p=Panel("5"+code,"Spatial cell labels" if code=="E" else "Coq8a detection in tissue",width,102)
    for i,age in enumerate(ages):
        d=meta[meta["sample"]==age];cx=(d.x_um.min()+d.x_um.max())/2;cy=(d.y_um.min()+d.y_um.max())/2
        assert len(d)==SAMPLES[age]["final_beads"]
        x=1+i*(side+4);p.text(x+side/2,20,age,ha="center")
        ax=p.ax(x,23,side,side);ax.set(xlim=(cx-span/2,cx+span/2),ylim=(cy+span/2,cy-span/2),xticks=[],yticks=[])
        for sp in ax.spines.values():sp.set_visible(True);sp.set_linewidth(.35)
        def dots(mask,color,sz,**kw):
            ax.scatter(d.loc[mask,"x_um"],d.loc[mask,"y_um"],s=sz,c=color,rasterized=True,linewidths=0,**kw)
        if code=="E":
            for cell,col in CELL_COLORS.items():dots(d.max_pred_celltype==cell,col,.12 if cell!="Myonuclei" else .06)
        else:
            injury=d.injury_classification=="injury";dots(~injury,"#e8e8e8",.05);dots(injury,"#c9c9c9",.08);dots(injury&(d.max_pred_celltype=="Fusing Myocytes"),"#f28e2b",.2)
            mu=injury&(d.max_pred_celltype=="MuSCs");ax.scatter(d.loc[mu,"x_um"],d.loc[mu,"y_um"],s=.7,facecolors="none",edgecolors="#2f6db0",lw=.12,rasterized=True)
            raw,pos=positive_barcodes(root/"data"/SAMPLES[age]["file"]);assert set(d.barcode)<=raw;dots(d.barcode.isin(pos),"#7b2cbf",.55)
        # Horizontal calibration uses the same square physical field in both sections.
        line=p.ax(x+side*.14,23+side+1,side*500/span,5);line.plot([0,1],[.5,.5],lw=1,c="k");line.set(xlim=(0,1),ylim=(0,1));line.axis("off")
        p.text(x+side*.14+side*500/span+2,23+side+5,"500 µm")
    if code=="E":
        labels=["B cells","Dendritic","Endothelial","Erythrocytes","FAPs","Fusing myo.","Mono./mac.","MuSCs","Myonuclei","NK cells","Neutrophils","Pericytes/SMC","Schwann/glial","T cells","Tenocytes"]
        for n,(col,lab) in enumerate(zip(CELL_COLORS.values(),labels)):p.key(151+60*(n//8),27+(n%8)*9.1,col,lab)
    else:
        for n,(col,lab) in enumerate([("#7b2cbf","Coq8a+ bead"),("#f28e2b","Fusing myocyte"),("#2f6db0","MuSC (ring)"),("#c9c9c9","Other injury"),("#e8e8e8","Outside injury")]):p.key(145,37+n*11,col,lab)
    p.save()
p=Panel("5G","Marker expression within injury",218,82)
t=tab("fusing_marker_audit.tsv");genes=t.gene.drop_duplicates().tolist();ax=p.ax(29,20,96,44)
for i,age in enumerate(ages):
    v=t[t["sample"]==age].set_index("gene").loc[genes]
    for j,suffix in enumerate(["fusing","other_injury"]):
        sc=ax.scatter(np.full(len(v),2*i+j),range(len(v)),s=2+24*v["detection_"+suffix],c=v["mean_log1p10k_"+suffix],cmap="magma",vmin=0,vmax=t[["mean_log1p10k_fusing","mean_log1p10k_other_injury"]].max().max(),linewidths=.15,edgecolors=".5")
ax.set(xlim=(-.5,3.5),ylim=(len(genes)-.5,-.5),xticks=range(4),yticks=range(len(genes)));ax.set_xticklabels(["Fusing","Other"]*2);ax.set_yticklabels(genes);ax.tick_params(length=0);p.text(52,19,"Young",ha="center");p.text(102,19,"Geriatric",ha="center")
p.text(132,22,"Mean ln(1 + CP10k)");p.fig.colorbar(sc,cax=p.ax(132,27,70,4),orientation="horizontal",ticks=[0,.6,1.2]);p.text(132,54,"Area: detection (%)")
for j,v in enumerate([10,30,50]):
    a=p.ax(133+25*j,57,10,10);a.scatter(0,0,s=2+24*v/100,c="#888");a.axis("off");p.text(142+25*j,67,str(v))
p.save()
p=Panel("5H","Coq8a detection by label",254,82);ax=p.ax(77,15,116,57)
labels=["Fusing Myocytes","Myonuclei","MuSCs","FAPs","Monocytes/Macrophages","Dendritic cells","Endothelial cells","Tenocytes"]
t=tab("coq8a_celltype_association.tsv")
for i,age in enumerate(ages):
    z=t[t["sample"]==age].set_index("cell_type").loc[labels];v=z.odds_ratio_plot
    ax.errorbar(v,np.arange(8)+(-.15 if i==0 else .15),xerr=np.vstack([v-z.odds_ratio_ci_low,z.odds_ratio_ci_high-v]),fmt="o",ms=2,c=styles[i],lw=.55)
ax.axvline(1,ls=":",c=".45",lw=.5);ax.set_xscale("log");ax.set(ylim=(7.5,-.5),yticks=range(8),xticks=[.1,1,10]);ax.set_xticklabels(["0.1","1","10"]);ax.set_yticklabels(["Fusing myocytes","Myonuclei","MuSC","FAP","Monocytes/mac.","Dendritic","Endothelial","Tenocytes"]);ax.tick_params(axis="y",length=0);p.text(202,28,"Odds ratio");p.text(202,39,"with 95% CI")
for i,age in enumerate(ages):p.key(202,55+11*i,styles[i],age)
p.save()
for code,w,name in [("I",218,"Coq8a detection by region"),("J",254,"Enrichment near fusing myocytes")]:
    p=Panel("5"+code,name,w,72);ax=p.ax(33,16,w-114,35)
    for i,age in enumerate(ages):
        if code=="I":
            z=tab("coq8a_region_detection.tsv");z=z[z["sample"]==age].set_index("region").loc[["Outside","Injury"]];ax.plot([0,1],z.coq8a_detection_pct,"o-",ms=3,c=styles[i],lw=.9);ax.set(xticks=[0,1],ylim=(.6,1.6),yticks=[.6,1,1.4]);ax.set_xticklabels(["Outside","Injury"]);ax.set_ylabel("Beads (%)")
        else:
            z=tab("coq8a_fusion_neighborhood.tsv");z=z[z["sample"]==age].sort_values("radius_um");ax.plot(z.radius_um,z.observed_to_random_ratio,"o-" if i==0 else "s--",ms=2.5,c=styles[i],lw=.9);sig=z[z.p_adj_bh_across_10_radius_tests<.05];ax.scatter(sig.radius_um,sig.observed_to_random_ratio+.015,s=12,marker="*",c=styles[i],lw=0)
            ax.set(xticks=[25,50,75,100,150],yticks=[1,1.1,1.2],ylim=(.98,1.28),xlabel="Radius (µm)",ylabel="Obs./random")
        p.key(w-71,28+10*i,styles[i],age)
    if code=="J":p.text(w-71,52,"* BH q < 0.05")
    p.save()
