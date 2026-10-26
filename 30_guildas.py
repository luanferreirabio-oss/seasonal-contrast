"""Reconstrucao das guildas de transporte de cations a partir da tabela curada COG-2024 (M1).

REGRAS EXPLICITAS DE INCLUSAO
  1. O COG deve constar da definicao curada COG-2024 (Galperin et al., 2025).
  2. O nome curado deve identificar simultaneamente (a) o cation e (b) uma funcao de
     transporte de membrana: transportador, permease, canal, subunidade de ATPase
     translocadora, receptor de membrana externa ou componente periplasmico de
     sistema de transporte.
  3. Sao EXCLUIDOS: subunidades reguladoras e sensoras citoplasmaticas; proteinas de
     ligacao sem funcao translocadora; enzimas de biossintese ou degradacao;
     simportadores que usam o gradiente do cation para mover outros solutos;
     proteinas de categoria R ou S (funcao mal caracterizada).
  4. Sistemas de efluxo/antiporte sao tabulados em grupo separado do de captacao.
"""
import pandas as pd, numpy as np, csv, json, itertools, math
from scipy import stats
pd.set_option('display.width',250)

DEF={r[0]:(r[1],r[2],r[3]) for r in csv.reader(open('/mnt/user-data/uploads/cog-24_def.tab',encoding='latin-1'),delimiter='\t') if len(r)>=4}

# ---------------- guildas ----------------
G={
 'K_captacao':      ['COG0168','COG3158','COG2060','COG2216','COG2156','COG1226'],
 'K_regulatorio':   ['COG0569','COG0490','COG1652','COG2205'],
 'K_efluxo':        ['COG0475','COG4651','COG0025','COG3263','COG0387'],
 'Mg_captacao':     ['COG0598','COG2239','COG0474','COG1285','COG3174'],
 'Mg_homeostase':   ['COG1253','COG4535','COG4536'],
 'Na_homeostase':   ['COG3004','COG3067','COG0025','COG3263','COG1757','COG1055','COG0387',
                     'COG0530','COG1006','COG1320','COG1863','COG2111','COG2212','COG0651',
                     'COG1668','COG4555'],
 'Fe_captacao':     ['COG0370','COG1918','COG0672','COG1914','COG3470',
                     'COG0609','COG0614','COG1120','COG4592','COG4594','COG4779',
                     'COG4604','COG4605','COG4606','COG4607',
                     'COG4771','COG4772','COG4773',
                     'COG1178','COG0390','COG4619','COG4558','COG4559',
                     'COG3182','COG2375','COG4114','COG4097'],
}
PRIM=['K_captacao','Mg_captacao','Na_homeostase','Fe_captacao']

print("="*104); print("COMPOSICAO DAS GUILDAS (tabela suplementar)"); print("="*104)
sup=[]
for g,cogs in G.items():
    print(f"\n### {g}  (n = {len(cogs)})")
    for c in cogs:
        cat,nome,gene=DEF.get(c,('?','NAO ENCONTRADO',''))
        print(f"  {c:<9}{cat:<5}{gene:<10}{nome[:70]}")
        sup.append(dict(guilda=g,COG=c,categoria=cat,gene=gene,descricao=nome))
pd.DataFrame(sup).to_csv('a3/TS5_guildas_composicao.csv',index=False)

# ---------------- abundancias ----------------
cog=pd.read_csv("a3/COG_matriz.csv",index_col=0)
ISL=[f'P{i}' for i in range(1,7)]
print(f"\n{len(cog)} COGs com abundancia (CPM); soma por amostra: "
      +', '.join(f'{i} {cog[i].sum():.0f}' for i in ISL))

ab={}; falt={}
for g,cogs in G.items():
    pres=[c for c in cogs if c in cog.index]; falt[g]=[c for c in cogs if c not in cog.index]
    ab[g]=cog.loc[pres,ISL].sum()
AB=pd.DataFrame(ab).T
print("\n"+"="*104); print("ABUNDANCIA DAS GUILDAS (CPM) - detectados / definidos"); print("="*104)
AB['media']=AB[ISL].mean(1); AB['CV_pct']=100*AB[ISL].std(1,ddof=1)/AB[ISL].mean(1)
AB['n_det']=[len(G[g])-len(falt[g]) for g in AB.index]; AB['n_def']=[len(G[g]) for g in AB.index]
print(AB.round(3).to_string())
for g in G:
    if falt[g]: print(f"  {g}: nao detectados -> {', '.join(falt[g])}")
AB.round(3).to_csv('a3/T6_guildas_abundancia.csv')

# ---------------- correlacao guilda x solo (p exato) ----------------
s=pd.read_csv('/mnt/user-data/outputs/solo_quimica_FINAL.csv')
si=s.groupby('ilha')[['K','Mg','Na','Fe']].mean().loc[ISL]
def exact_p(x,y):
    r=stats.spearmanr(x,y).statistic; n=len(x); cnt=0
    for pm in itertools.permutations(range(n)):
        if abs(stats.spearmanr(x,[y[i] for i in pm]).statistic)>=abs(r)-1e-12: cnt+=1
    return r,cnt/math.factorial(n)
print("\n"+"="*104); print("GUILDA vs TEOR NO SOLO - Spearman com p exato (720 permutacoes, n=6)"); print("="*104)
rows=[]
MAP={'K_captacao':'K','Mg_captacao':'Mg','Na_homeostase':'Na','Fe_captacao':'Fe'}
for g,el in MAP.items():
    x=AB.loc[g,ISL].values.astype(float); y=si[el].values.astype(float)
    r,p=exact_p(x,y)
    rows.append(dict(guilda=g,elemento=el,n_COG=AB.loc[g,'n_det'],
                     CPM_medio=round(AB.loc[g,'media'],2),CV_guilda=round(AB.loc[g,'CV_pct'],1),
                     CV_solo=round(100*y.std(ddof=1)/y.mean(),1),rho=round(r,3),p_exato=round(p,4)))
gt=pd.DataFrame(rows); print(gt.to_string(index=False)); gt.to_csv('a3/T7_guilda_solo.csv',index=False)
print(f"\n  significativas (p<0,05): {(gt.p_exato<0.05).sum()} de {len(gt)}")
print("  |rho| minimo para p<0,05 com n=6 por permutacao exata: 0,886")
json.dump(dict(guildas={g:dict(n=int(AB.loc[g,'n_det']),cpm=round(float(AB.loc[g,'media']),2),
                               cv=round(float(AB.loc[g,'CV_pct']),1)) for g in PRIM},
               corr=gt.to_dict('records')),open('a3/R_guildas.json','w'),indent=1)

# ---------------- normalizacao por equivalentes de genoma ----------------
MK=['COG0048','COG0049','COG0052','COG0080','COG0081','COG0087','COG0088','COG0090','COG0091',
    'COG0092','COG0093','COG0094','COG0096','COG0097','COG0098','COG0099','COG0100','COG0102',
    'COG0103','COG0184','COG0185','COG0186','COG0197','COG0200','COG0244','COG0256','COG0522']
ge=cog.loc[[m for m in MK if m in cog.index],ISL].median()
print("\n"+"="*104); print("GUILDAS EM COPIAS POR GENOMA (normalizacao por 27 marcadores ribossomais)"); print("="*104)
CG=AB.loc[PRIM,ISL].div(ge,axis=1)
CG['media']=CG.mean(axis=1); CG['CV_pct']=100*CG[ISL].std(axis=1,ddof=1)/CG[ISL].mean(axis=1)
print(CG.round(3).to_string())
CG.round(3).to_csv('a3/T6b_guildas_por_genoma.csv')
print("\nCorrelacao guilda-solo com a guilda em copias por genoma:")
for g,el in MAP.items():
    x=CG.loc[g,ISL].values.astype(float); y=si[el].values.astype(float)
    r,p=exact_p(x,y); print(f"  {g:<15} vs {el:<3} rho={r:+.3f}  p={p:.4f}")
