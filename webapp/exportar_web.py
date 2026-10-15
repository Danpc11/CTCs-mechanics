"""Genera webapp/index.html desde data/poblaciones.yaml y webapp/curvas_fem.json.

    python webapp/exportar_web.py

Para regenerar curvas_fem.json con otro dispositivo o malla, ver webapp/README.md.
"""
import json, pathlib, yaml
H = pathlib.Path(__file__).resolve().parent
R = H.parent
pob = yaml.safe_load((R / 'data' / 'poblaciones.yaml').read_text(encoding='utf-8'))
NOMBRES = {'neutrofilo': 'Neutrófilo', 'linfocito': 'Linfocito', 'monocito': 'Monocito', 'eosinofilo': 'Eosinófilo',
           'basofilo': 'Basófilo'}
def m(s): return s['valor'] if s['dist'] == 'fijo' else s['media']
fen = []
for k, p in pob['poblaciones'].items():
    fen.append(dict(id=k, nombre=NOMBRES.get(k, k), grupo=p['grupo'], D=round(m(p['D']), 2), E_app=round(m(p['E_app']), 3),
                    gamma_pN_um=round(1000 * m(p['gamma']), 1), tau=round(m(p['tau']), 4), an_a=round(m(p['an_a']), 3),
                    fuente=f"E_app: {p['E_app']['estado']}; diámetro: {p['D']['estado']}; γ: {p['gamma']['estado']}."))
cur = json.loads((H / 'curvas_fem.json').read_text())
NOMDEV = {'membrana_PC': ('Membrana track-etched de policarbonato (10 µm)', 'membrana'),
          'canal_PDMS': ('Constricción en canal de PDMS', 'canal')}
dev = {k: dict(nombre=NOMDEV[k][0], tipo=NOMDEV[k][1], perfil=v['perfil'],
               u=[round(x, 5) for x in v['u']], g=[round(x, 7) for x in v['g']]) for k, v in cur.items()}
datos = dict(K_E=pob['modelo']['K_E'], fenotipos=fen, dispositivos=dev)
nucleo = (H / 'nucleo.js').read_text().replace("if (typeof module !== 'undefined') module.exports", "// export: ")
html = (H / 'plantilla.html').read_text(encoding='utf-8').replace('__DATOS__', json.dumps(datos, ensure_ascii=False)).replace('__NUCLEO__', nucleo)
(H / 'index.html').write_text(html, encoding='utf-8')
print('index.html:', len(html) // 1024, 'kB;', len(fen), 'fenotipos;', len(dev), 'dispositivos')
