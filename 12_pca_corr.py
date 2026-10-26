import pandas as pd, numpy as np, itertools, json
from scipy import stats
pd.set_option('display.width',260)
w=pd.read_csv('a3/agua_corrigida.csv')
RET=['Na','K','Ca','Mg','Mn','Cu','Ba','Zn','Sr','SO4','Cl','H4SiO4','Dureza','Alcalinidade','COT']
# imputa as 3 celulas ausentes (SO4 e Cl de C1) pela mediana da campanha, so para a PCA
wp=w.copy()
for v in ['SO4','Cl']:
    for c in ['C1','C2']:
        m=wp.loc[wp.camp==c,v].median()
        wp.loc[(wp.camp==c)&(wp[v].isna()),v]=m

def pca(df,cols):
    X=df[cols].astype(float).values; Z=(X-X.mean(0))/X.std(0,ddof=1)
    U,S,Vt=np.linalg.svd(Z,full_matrices=False); ev=S**2/(len(Z)-1)
    return 100*ev/ev.sum(), Vt, U*S

print("="*80); print("PCA - sensibilidade a dureza e ao cobre"); print("="*80)
VAR={'15 var (com dureza e Cu)':RET,
     '14 var (sem dureza)':[c for c in RET if c!='Dureza'],
     '13 var (sem dureza e sem Cu)':[c for c in RET if c not in ('Dureza','Cu')]}
for nm,cols in VAR.items():
    expl,Vt,sc=pca(wp,cols)
    s1=sc[(wp.camp=='C1').values,0]; s2=sc[(wp.camp=='C2').values,0]
    sM=sc[(wp.est=='M').values,0];  sJ=sc[(wp.est=='J').values,0]
    print(f"\n{nm}: PC1={expl[0]:.1f}% PC2={expl[1]:.1f}% (soma {expl[0]+expl[1]:.1f}%)")
    l=pd.Series(Vt[0],index=cols).sort_values(key=abs,ascending=False)
    print("  PC1:", ', '.join(f'{k} {v:+.3f}' for k,v in list(l.items())[:6]))
    print(f"  PC1 mediana: C1={np.median(s1):+.2f} C2={np.median(s2):+.2f} | M={np.median(sM):+.2f} J={np.median(sJ):+.2f}")
    if nm.startswith('14'):
        out=pd.DataFrame(sc[:,:2],columns=['PC1','PC2']); out['camp']=wp.camp.values; out['est']=wp.est.values
        out.to_csv('a3/T_PCA_scores.csv',index=False)
        pd.DataFrame(Vt[:2].T,index=cols,columns=['PC1','PC2']).round(3).to_csv('a3/TS3_PCA_cargas.csv')
        json.dump(dict(pc1=round(expl[0],1),pc2=round(expl[1],1),
                       med_C1=round(float(np.median(s1)),2),med_C2=round(float(np.median(s2)),2),
                       med_M=round(float(np.median(sM)),2),med_J=round(float(np.median(sJ)),2),
                       cargas={k:round(float(v),3) for k,v in list(l.items())[:6]}),
                  open('a3/R_pca.json','w'),indent=1)

print("\n"+"="*80); print("CORRELACOES: agrupada (n=51) vs dentro de cada campanha"); print("="*80)
rows=[]
for a,b in itertools.combinations(RET,2):
    d=w[[a,b,'camp']].dropna()
    rp,pp=stats.spearmanr(d[a],d[b])
    d1=d[d.camp=='C1']; d2=d[d.camp=='C2']
    r1,p1=stats.spearmanr(d1[a],d1[b]); r2,p2=stats.spearmanr(d2[a],d2[b])
    rows.append(dict(par=f'{a}-{b}',rho=round(rp,3),p=round(pp,4),
                     rho_C1=round(r1,3),p_C1=round(p1,4),rho_C2=round(r2,3),p_C2=round(p2,4)))
co=pd.DataFrame(rows)
# FDR sobre os 105 pares agrupados
for col,q in [('p','q'),('p_C1','q_C1'),('p_C2','q_C2')]:
    pv=co[col].values; m=len(pv); o=np.argsort(pv)
    qq=np.minimum.accumulate((pv[o]*m/(np.arange(m)+1))[::-1])[::-1]
    f=np.empty(m); f[o]=np.minimum(qq,1.0); co[q]=f.round(3)
co.to_csv('a3/TS4_correlacoes.csv',index=False)
print("\n-- pares usados na Discussao --")
for pr in ['Cu-COT','K-COT','K-Cu','Mn-COT','Ca-Dureza','Ca-Mg','K-Mg','Na-Mn','Mn-Cu','Alcalinidade-COT']:
    r=co[co.par==pr]
    if len(r): print(r[['par','rho','q','rho_C1','p_C1','rho_C2','p_C2']].to_string(index=False,header=(pr=='Cu-COT')))
BL=['Na','K','Ca','Mg','Sr','Ba','Mn','Dureza']
sub=co[co.par.apply(lambda s: all(p in BL for p in s.split('-')))]
print(f"\n-- bloco cations basicos + Mn + dureza ({len(sub)} pares) --")
print(f"  agrupado: rho {sub.rho.min():.3f} a {sub.rho.max():.3f} | todos positivos: {(sub.rho>0).all()}")
print(f"  dentro de C1: {sub.rho_C1.min():.3f} a {sub.rho_C1.max():.3f}")
print(f"  dentro de C2: {sub.rho_C2.min():.3f} a {sub.rho_C2.max():.3f}")
print(f"  significativos apos FDR (agrupado): {(sub.q<0.05).sum()}/{len(sub)}")
inv=co[(np.sign(co.rho)!=np.sign(co.rho_C1))|(np.sign(co.rho)!=np.sign(co.rho_C2))]
print(f"\n-- pares em que o sinal se inverte dentro de alguma campanha: {len(inv)} de {len(co)} --")
json.dump(dict(bloco_min=float(sub.rho.min()),bloco_max=float(sub.rho.max()),
               bloco_C1=[float(sub.rho_C1.min()),float(sub.rho_C1.max())],
               bloco_C2=[float(sub.rho_C2.min()),float(sub.rho_C2.max())],
               n_inv=len(inv),n_pares=len(co)),open('a3/R_corr.json','w'),indent=1)
