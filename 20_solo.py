"""Tabela de solo, razoes molares e coeficientes de seletividade (S2, S4, S5, S6)."""
import pandas as pd, numpy as np, json
pd.set_option('display.width',280); pd.set_option('display.max_columns',40)
s=pd.read_csv('/mnt/user-data/outputs/solo_quimica_FINAL.csv')
w=pd.read_csv('a3/agua_corrigida.csv')
NOME={'P1':'Santana','P2':'Santaninha','P3':'Angapijó','P4':'Conceição','P5':'São Joaquim','P6':'Tauaré'}
s['Ilha']=s.ilha.map({f'P{i}':i for i in range(1,7)})

print("="*100); print("TABELA DE SOLO POR ILHA E PROFUNDIDADE (S2)"); print("="*100)
COLS=['pH','MO','N','P','Ca','Mg','K','Na','Al','Acidez','S','CTC','V','Fe','Mn','Cu','Zn']
t=s[['Ilha','prof']+COLS].sort_values(['Ilha','prof'])
print(t.to_string(index=False))
t.to_csv('a3/T3_solo_completa.csv',index=False)
print("\n-- medias do perfil 0-20 cm e CV entre ilhas --")
si=s.groupby('Ilha')[COLS].mean()
print(si.round(3).to_string())
cvs=(100*si.std(ddof=1)/si.mean()).round(1)
print("\nCV entre ilhas (%):"); print(cvs.to_string())
si.round(3).to_csv('a3/T3b_solo_media_perfil.csv')

# ---------------- razoes molares ----------------
mol=pd.DataFrame({'Ca':si.Ca/2,'Mg':si.Mg/2,'K':si.K,'Na':si.Na})
mol['Ca/Mg']=mol.Ca/mol.Mg; mol['Ca/K']=mol.Ca/mol.K; mol['Na/K']=mol.Na/mol.K
wr=pd.DataFrame({'Ca/Mg':w.Ca/w.Mg,'Ca/K':w.Ca/w.K,'Na/K':w.Na/w.K,'camp':w.camp})
print("\n"+"="*100); print("RAZOES MOLARES"); print("="*100)
rows=[]
for c,lab in [('C1','Água C1 (menos chuvoso)'),('C2','Água C2 (chuvoso)')]:
    sub=wr[wr.camp==c]
    rows.append(dict(Compartimento=lab,
        **{k:f"{sub[k].median():.3f} ({sub[k].min():.3f}-{sub[k].max():.3f})" for k in ['Ca/Mg','Ca/K','Na/K']}))
for i,r in mol.iterrows():
    rows.append(dict(Compartimento=f"Solo - Ilha {i} ({NOME['P%d'%i]})",
        **{k:f"{r[k]:.3f}" for k in ['Ca/Mg','Ca/K','Na/K']}))
rt=pd.DataFrame(rows); print(rt.to_string(index=False)); rt.to_csv('a3/T4_razoes.csv',index=False)

print("\nSolo (6 ilhas):  Ca/Mg %.3f-%.3f | Ca/K %.3f-%.3f | Na/K %.3f-%.3f"%(
    mol['Ca/Mg'].min(),mol['Ca/Mg'].max(),mol['Ca/K'].min(),mol['Ca/K'].max(),mol['Na/K'].min(),mol['Na/K'].max()))
for c in ['C1','C2']:
    sub=wr[wr.camp==c]
    print("Agua %s (n=%d): Ca/Mg %.3f-%.3f | Ca/K %.3f-%.3f | Na/K %.3f-%.3f  [medianas %.3f / %.3f / %.3f]"%(
        c,len(sub),sub['Ca/Mg'].min(),sub['Ca/Mg'].max(),sub['Ca/K'].min(),sub['Ca/K'].max(),
        sub['Na/K'].min(),sub['Na/K'].max(),sub['Ca/Mg'].median(),sub['Ca/K'].median(),sub['Na/K'].median()))

# ---------------- seletividade ----------------
print("\n"+"="*100); print("SELETIVIDADE APARENTE E COEFICIENTE DE GAPON"); print("="*100)
out={}
for c in ['C1','C2']:
    sub=w[w.camp==c]
    Ca_w=sub.Ca.median()*1e-6; K_w=sub.K.median()*1e-6
    Mg_w=sub.Mg.median()*1e-6; Na_w=sub.Na.median()*1e-6
    print(f"\n--- solucao equilibrante = mediana da campanha {c} "
          f"(Ca {Ca_w*1e6:.2f}, Mg {Mg_w*1e6:.2f}, K {K_w*1e6:.2f}, Na {Na_w*1e6:.2f} umol/L) ---")
    rows=[]
    for i in range(1,7):
        r=si.loc[i]
        # homovalentes: razao de razoes molares (adimensional)
        selCaMg=(r.Ca/2/(r.Mg/2))/(Ca_w/Mg_w)
        selKNa=(r.K/r.Na)/(K_w/Na_w)
        # heterovalente K-Ca: coeficiente de Gapon, (L/mol)^0.5, base de carga
        KG=(r.K/r.Ca)*(np.sqrt(Ca_w)/K_w)
        # fator de enriquecimento Ca/K (razao de razoes molares, nao e coef. de seletividade)
        enr=((r.Ca/2)/r.K)/(Ca_w/K_w)
        rows.append(dict(Ilha=i,sel_CaMg=round(selCaMg,3),sel_KNa=round(selKNa,3),
                         K_Gapon=round(KG,2),enriq_CaK=round(enr,2)))
    d=pd.DataFrame(rows); print(d.to_string(index=False))
    print(f"  seletividade aparente Ca/Mg: {d.sel_CaMg.min():.3f} a {d.sel_CaMg.max():.3f}  (<1 = solo favorece Mg)")
    print(f"  seletividade aparente K/Na : {d.sel_KNa.min():.2f} a {d.sel_KNa.max():.2f}   (>1 = solo favorece K)")
    print(f"  coeficiente de Gapon K-Ca  : {d.K_Gapon.min():.2f} a {d.K_Gapon.max():.2f} (L/mol)^0,5")
    print(f"  enriquecimento da razao Ca/K: {d.enriq_CaK.min():.1f} a {d.enriq_CaK.max():.1f} vezes")
    out[c]=d; d.to_csv(f'a3/T5_seletividade_{c}.csv',index=False)

# faixa global do enriquecimento Ca/K contra as medianas das duas campanhas
alle=pd.concat([out['C1'],out['C2']])
print(f"\n>> ENRIQUECIMENTO Ca/K sobre as duas campanhas: {alle.enriq_CaK.min():.1f} a {alle.enriq_CaK.max():.1f} vezes")
print(f">> seletividade Ca/Mg global: {alle.sel_CaMg.min():.3f} a {alle.sel_CaMg.max():.3f}")
print(f">> seletividade K/Na global : {alle.sel_KNa.min():.1f} a {alle.sel_KNa.max():.1f}")
print(f">> Gapon K-Ca global        : {alle.K_Gapon.min():.2f} a {alle.K_Gapon.max():.2f}")
# reducao Na/K agua->solo
red=[wr.loc[wr.camp==c,'Na/K'].median()/mol['Na/K'] for c in ['C1','C2']]
allred=pd.concat(red)
print(f">> reducao da razao Na/K agua->solo: {allred.min():.1f} a {allred.max():.1f} vezes")
json.dump(dict(solo_CaMg=[round(mol['Ca/Mg'].min(),3),round(mol['Ca/Mg'].max(),3)],
  solo_CaK=[round(mol['Ca/K'].min(),3),round(mol['Ca/K'].max(),3)],
  solo_NaK=[round(mol['Na/K'].min(),3),round(mol['Na/K'].max(),3)],
  agua_CaMg={c:[round(wr.loc[wr.camp==c,'Ca/Mg'].min(),3),round(wr.loc[wr.camp==c,'Ca/Mg'].max(),3),round(wr.loc[wr.camp==c,'Ca/Mg'].median(),3)] for c in ['C1','C2']},
  agua_CaK={c:[round(wr.loc[wr.camp==c,'Ca/K'].min(),3),round(wr.loc[wr.camp==c,'Ca/K'].max(),3),round(wr.loc[wr.camp==c,'Ca/K'].median(),3)] for c in ['C1','C2']},
  agua_NaK={c:[round(wr.loc[wr.camp==c,'Na/K'].min(),3),round(wr.loc[wr.camp==c,'Na/K'].max(),3),round(wr.loc[wr.camp==c,'Na/K'].median(),3)] for c in ['C1','C2']},
  enriq_CaK=[round(alle.enriq_CaK.min(),1),round(alle.enriq_CaK.max(),1)],
  sel_CaMg=[round(alle.sel_CaMg.min(),3),round(alle.sel_CaMg.max(),3)],
  sel_KNa=[round(alle.sel_KNa.min(),1),round(alle.sel_KNa.max(),1)],
  gapon=[round(alle.K_Gapon.min(),2),round(alle.K_Gapon.max(),2)],
  red_NaK=[round(allred.min(),1),round(allred.max(),1)],
  cv_solo={k:float(v) for k,v in cvs.items()}),open('a3/R_solo.json','w'),indent=1)
