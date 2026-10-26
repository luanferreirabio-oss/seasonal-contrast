"""Constroi o dataset de agua corrigido contra a fonte primaria (Cordeiro, 2024, Tabela 5)."""
import pandas as pd, numpy as np
w=pd.read_csv('a3/agua.csv')
w['est']=w.Ponto.map({'Ponto 2':'M','Ponto 1':'J'})
log=[]
# 1) SO4 C1 = 687,07 -> invalido (a dissertacao reporta max 77,04; n=24)
i=w.index[(w.camp=='C1')&(w.SO4>600)][0]; log.append(('SO4','C1',w.at[i,'SO4'],'excluido')); w.at[i,'SO4']=np.nan
# 2) Cl C1 = 817,98 -> invalido (a dissertacao reporta max 70,52; n=24)
i=w.index[(w.camp=='C1')&(w.Cl>600)][0]; log.append(('Cl','C1',w.at[i,'Cl'],'excluido')); w.at[i,'Cl']=np.nan
# 3) COD C2 = 21,00 -> 0 (abaixo do LQ; reproduz min=0 max=4,40 media=3,78 dp=0,80)
i=w.index[(w.camp=='C2')&(w.COD>10)][0]; log.append(('COD','C2',w.at[i,'COD'],'-> 0 (<LQ)')); w.at[i,'COD']=0.0
# 4) Fe: 100% abaixo do LQ (0,100 mg/L) em C1 -> fora do conjunto retido
w=w.drop(columns=['Fe'])
w=w.rename(columns={'COD':'COT'})   # metodo SMWW 5310 B em amostra nao filtrada = C organico total
w.to_csv('a3/agua_corrigida.csv',index=False)
print("Correcoes aplicadas:")
for v,c,old,acao in log: print(f"  {v:<5} {c}  {old:>8.2f}  {acao}")
print("  Fe    removido da agua (100% < LQ em C1)")
print("  COD   renomeado COT (SMWW 5310 B, amostra nao filtrada)")
print(f"\nn = {len(w)}  |  variaveis retidas = 15")
RET=['Na','K','Ca','Mg','Mn','Cu','Ba','Zn','Sr','SO4','Cl','H4SiO4','Dureza','Alcalinidade','COT']
print("  ", ', '.join(RET))
