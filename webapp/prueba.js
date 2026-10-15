const N = require('./nucleo.js'), C = require('./referencia_python.json');
let ok = true;
for (const [nd, est] of [['canal_PDMS', 'local'], ['membrana_PC', 'trabajo']]) {
  const c = C[nd], cv = N.curva(c), pf = N.perfil(c.perfil);
  for (const [D, ref] of Object.entries(c.ref)) {
    const Rn = D / 2; if (Rn / c.a < cv.xmin) { console.log(nd, D, 'um fuera de la curva de prueba: excluido (Python da cota inferior)'); continue; }
    const js = est === 'local' ? N.dPlocal(pf, Rn, c.a, c.Gs, cv) : N.dPtrabajo(Rn, c.a, c.Gs, cv);
    const py = est === 'local' ? ref[0] : ref[1], e = Math.abs(js - py) / py;
    if (e > 0.03) ok = false;
    console.log(nd, est, D, 'um  JS', js.toFixed(4), ' Python', py, ' error', (100 * e).toFixed(2) + '%');
  }
}
console.log(ok ? 'PRUEBA SUPERADA (<3 %)' : 'PRUEBA FALLIDA');
