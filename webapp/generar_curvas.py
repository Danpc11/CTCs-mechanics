import sys,json,yaml,pandas as pd; import pathlib,os; R=pathlib.Path(__file__).resolve().parents[1]; os.chdir(R); sys.path.insert(0,'src')
from ctcfem import geometria as G
from ctcfem.modelo_celula import cerrar_apertura
from ctcfem import poblaciones as P; c=P.muestrear(P.cargar('data/poblaciones.yaml'),'MCF-7',1,20261015).iloc[0].to_dict()
mats=yaml.safe_load(open('data/materiales.yaml'))
out={}
for nd in ['membrana_PC','canal_PDMS']:
    disp=mats['dispositivos'][nd]; mat=mats['materiales'][disp['material']]
    perfil=G.perfil_de_dispositivo(disp)
    radios=[c['a']/x for x in [1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4]]
    l_w=disp['espesor_um']['nominal'] if disp['tipo']=='membrana' else None
    r=cerrar_apertura(c,perfil,radios,conf_max_intento=2.4,E_pared_kPa=mat['E_GPa']['nominal']*1e6,nu_pared=mat['nu'],l_pared=l_w)
    a=c['a']; Gs=c['Gs']; cur=sorted(r['curva'],key=lambda p:p['rho'])
    out[nd]=dict(u=[p['rho']/a for p in cur],g=[p['W_total']/(Gs*a*a) for p in cur],perfil=dict(zt=perfil.zt,L=perfil.L,R_ext=perfil.R_ext),
                 conf_max=a/r['rho_min_alcanzado'],motivo=r['motivo_parada'],seg=r['segundos'])
    print(nd,'conf max alcanzado',round(a/r['rho_min_alcanzado'],3),r['motivo_parada'],r['segundos'],'s',flush=True)
json.dump(out,open('webapp/curvas_fem.json','w'),indent=0)
