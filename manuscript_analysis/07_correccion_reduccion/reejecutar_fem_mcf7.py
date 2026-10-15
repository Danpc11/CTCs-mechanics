import sys,time,json; sys.path.insert(0,'../01_solver'); sys.path.insert(0,'../03_modalidades'); sys.path.insert(0,'.')
import paso1_fem as P
for c in [1.70,1.42,1.20,1.30,1.55,1.90,2.10,2.43,1.12]:
    t=time.time()
    r=P.tarea(("MCF-7","referencia",4.0,5.0,c))
    keep={k:v for k,v in r.items() if isinstance(v,(int,float,str,bool)) or v is None}
    keep['seg_total']=time.time()-t
    open('../datos/fem_mcf7_referencia_reejecucion.jsonl','a').write(json.dumps(keep)+'\n')
