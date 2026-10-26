import pandas as pd, numpy as np, itertools, json
from scipy import stats
pd.set_option('display.width',260); pd.set_option('display.max_columns',70)
w=pd.read_csv('a3/agua_corrigida.csv')
RET=['Na','K','Ca','Mg','Mn','Cu','Ba','Zn','Sr','SO4','Cl','H4SiO4','Dureza','Alcalinidade','COT']
R={}

# ================= 1. censura (% < LQ) =================
rows=[]
for v in RET:
    r={'var':v}
    for c in ['C1','C2']:
        s=w.loc[w.camp==c,v]
        n=s.notna().sum(); z=int((s==0).sum())
        r[f'{c}_n']=n; r[f'{c}_det']=round(100*(n-z)/n,1)
    rows.append(r)
cens=pd.DataFrame(rows); cens.to_csv('a3/TS1_deteccao.csv',index=False)
print("="*80); print("1. DETECCAO (% de amostras >= LQ)"); print("="*80); print(cens.to_string(index=False))

# ================= 2. CV intra-ciclo =================
rows=[]
for v in RET:
    r={'var':v}
    for c in ['C1','C2']:
        for e in ['M','J']:
            s=w.loc[(w.camp==c)&(w.est==e),v].dropna().astype(float)
            r[f'{c}{e}']=round(100*s.std(ddof=1)/s.mean(),1) if s.mean()>0 else np.nan
    rows.append(r)
cv=pd.DataFrame(rows); cv.to_csv('a3/TS2_CV_intraciclo.csv',index=False)
print("\n"+"="*80); print("2. CV INTRA-CICLO DE MARE (%) - dados corrigidos"); print("="*80)
print(cv.to_string(index=False))
allcv=cv[['C1M','C1J','C2M','C2J']].values.ravel(); allcv=allcv[~np.isnan(allcv)]
print(f"\n  faixa global: {allcv.min():.1f}% a {allcv.max():.1f}%")
BAS=['Na','K','Ca','Mg']; sub=cv[cv['var'].isin(BAS)][['C1M','C1J','C2M','C2J']].values.ravel()
print(f"  cations basicos (Na,K,Ca,Mg): {sub.min():.1f}% a {sub.max():.1f}%")
for g,nm in [(['Dureza','Alcalinidade'],'dureza/alcalinidade'),(['Mn','Zn'],'Mn/Zn'),(['Cu'],'Cu'),
             (['SO4','Cl','H4SiO4'],'anions/silica'),(['COT'],'COT'),(['Ba','Sr'],'Ba/Sr')]:
    s=cv[cv['var'].isin(g)][['C1M','C1J','C2M','C2J']].values.ravel(); s=s[~np.isnan(s)]
    print(f"  {nm}: {s.min():.1f}% a {s.max():.1f}%")
R['cv_faixa']=[float(allcv.min()),float(allcv.max())]

# ================= 3. contraste estacoes + FDR =================
rows=[]
for v in RET:
    for c in ['C1','C2']:
        sub=w[w.camp==c]
        up=sub[sub.est=='M'].set_index('h')[v]; dn=sub[sub.est=='J'].set_index('h')[v]
        idx=up.index.intersection(dn.index)
        pair=pd.DataFrame({'u':up.loc[idx],'d':dn.loc[idx]}).dropna()
        u,d=pair.u.astype(float),pair.d.astype(float)
        try: _,p=stats.wilcoxon(u,d)
        except ValueError: p=np.nan
        rows.append(dict(var=v,camp=c,n=len(pair),med_M=round(u.median(),2),med_J=round(d.median(),2),p=p))
st=pd.DataFrame(rows)
pv=st.p.values; m=len(pv); o=np.argsort(pv)
q=np.minimum.accumulate((pv[o]*m/(np.arange(m)+1))[::-1])[::-1]
fdr=np.empty(m); fdr[o]=np.minimum(q,1.0)
st['q']=fdr.round(3); st['p']=st.p.round(4)
st.to_csv('a3/T1_estacoes.csv',index=False)
print("\n"+"="*80); print("3. CONTRASTE ENTRE ESTACOES + FDR (dados corrigidos)"); print("="*80)
print(st.to_string(index=False))
print(f"\n  {m} comparacoes | p<0,05: {(pv<0.05).sum()} | q_BH<0,05: {(fdr<0.05).sum()} | menor q = {fdr.min():.3f}")
print(f"  Bonferroni (limiar {0.05/m:.4f}): {(pv<0.05/m).sum()} significativas")
sig=st[st.p<0.05]
print("  p<0,05 brutos:", ', '.join(f"{r['var']} {r.camp} (p={r.p:.3f}, q={r.q:.2f})" for _,r in sig.iterrows()))
pC1=st.loc[st.camp=='C1','p']; print(f"  faixa de p no periodo menos chuvoso: {pC1.min():.3f} a {pC1.max():.3f}")
R['n_sig_bruto']=int((pv<0.05).sum()); R['q_min']=float(fdr.min())

# ================= 4. contraste entre campanhas (tamanho de efeito) =================
def cliff(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    gt=sum((y>x) for x in a for y in b); lt=sum((y<x) for x in a for y in b)
    return (gt-lt)/(len(a)*len(b))
rows=[]
for v in RET:
    a=w.loc[w.camp=='C1',v].dropna().astype(float); b=w.loc[w.camp=='C2',v].dropna().astype(float)
    q1a,q3a=np.percentile(a,[25,75]); q1b,q3b=np.percentile(b,[25,75])
    rows.append(dict(var=v,n_C1=len(a),med_C1=round(a.median(),2),AIQ_C1=f'{q1a:.2f}-{q3a:.2f}',
                     n_C2=len(b),med_C2=round(b.median(),2),AIQ_C2=f'{q1b:.2f}-{q3b:.2f}',
                     var_pct=round(100*(b.median()-a.median())/a.median(),1) if a.median()>0 else np.nan,
                     cliff=round(cliff(a,b),3)))
se=pd.DataFrame(rows).sort_values('cliff',key=abs,ascending=False)
se.to_csv('a3/T2_campanhas.csv',index=False)
print("\n"+"="*80); print("4. CONTRASTE ENTRE CAMPANHAS (delta de Cliff: chuvoso vs menos chuvoso)"); print("="*80)
print(se.to_string(index=False))
gr=se[se.cliff.abs()>=.474]; un=se[se.cliff.abs()==1]
print(f"\n  |delta|>=0,474 (grande): {len(gr)} de 15 -> {', '.join(gr['var'])}")
print(f"  |delta|=1,000 (separacao completa): {len(un)} -> {', '.join(un['var'])}")
R['n_grande']=len(gr); R['uni']=list(un['var'])
json.dump(R,open('a3/R_agua.json','w'),indent=1)
