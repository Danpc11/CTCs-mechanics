// The JavaScript core reproduces the Python reduction (node webapp/test_core.js).
const N = require('./core.js'), C = require('./python_reference.json');
let ok = true;
for (const [dev, est] of [['channel_PDMS', 'local'], ['membrane_PC', 'work']]) {
  const c = C[dev], cv = N.energyCurve(c), pf = N.wallProfile(c.perfil);
  for (const [D, ref] of Object.entries(c.ref)) {
    const Rn = D / 2; if (Rn / c.a < cv.xmin) continue;           // outside the reference curve
    const js = est === 'local' ? N.dPlocal(pf, Rn, c.a, c.Gs, cv) : N.dPwork(Rn, c.a, c.Gs, cv);
    const py = est === 'local' ? ref[0] : ref[1], e = Math.abs(js - py) / py;
    if (e > 0.03) ok = false;
    console.log(`${dev} ${est} ${D} um  JS ${js.toFixed(4)}  Python ${py}  error ${(100 * e).toFixed(2)}%`);
  }
}
console.log(ok ? 'PASSED (< 3 %)' : 'FAILED'); process.exitCode = ok ? 0 : 1;
