// Calculation core of the web app (no dependencies). Same algorithm as src/ctcfem/reduction.py, using the
// verified collapse W = Gs a^2 g(rho/a) with g tabulated from the finite-element model.
function pchip(x, y) {                       // monotone cubic interpolation (Fritsch-Carlson)
  const n = x.length, h = [], d = [], m = new Array(n);
  for (let i = 0; i < n - 1; i++) { h.push(x[i + 1] - x[i]); d.push((y[i + 1] - y[i]) / h[i]); }
  m[0] = d[0]; m[n - 1] = d[n - 2];
  for (let i = 1; i < n - 1; i++) {
    if (d[i - 1] * d[i] <= 0) m[i] = 0;
    else { const w1 = 2 * h[i] + h[i - 1], w2 = h[i] + 2 * h[i - 1]; m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i]); }
  }
  const seg = t => { let i = 0; while (i < n - 2 && t > x[i + 1]) i++; return i; };
  return {
    f(t) { const i = seg(t), s = (t - x[i]) / h[i], s2 = s * s, s3 = s2 * s;
      return (2*s3-3*s2+1)*y[i] + (s3-2*s2+s)*h[i]*m[i] + (-2*s3+3*s2)*y[i+1] + (s3-s2)*h[i]*m[i+1]; },
    df(t) { const i = seg(t), s = (t - x[i]) / h[i], s2 = s * s;
      return ((6*s2-6*s)*y[i] + (3*s2-4*s+1)*h[i]*m[i] + (-6*s2+6*s)*y[i+1] + (3*s2-2*s)*h[i]*m[i+1]) / h[i]; },
    xmin: x[0], xmax: x[n - 1]
  };
}
function energyCurve(c) {                    // c.u ascending, c.g; adds g(1) = 0 and enforces monotonicity
  const u = [], g = [];
  c.u.forEach((v, i) => { if (v < 1) { u.push(v); g.push(Math.max(c.g[i], 0)); } });
  u.push(1); g.push(0);
  for (let i = g.length - 2; i >= 0; i--) g[i] = Math.max(g[i], g[i + 1]);
  return pchip(u, g);
}
function wallProfile(p) {
  const G = z => 1 + 0.5 * (Math.tanh((z - p.zt) / p.L) - Math.tanh((z + p.zt) / p.L)), g0 = G(0);
  return {
    Rw: (z, Rn) => Rn + (p.R_ext - Rn) * (G(Math.abs(z)) - g0) / (1 - g0),
    dRw: (z, Rn) => { const s1 = 1 / Math.cosh((z - p.zt) / p.L) ** 2, s2 = 1 / Math.cosh((z + p.zt) / p.L) ** 2;
      return (p.R_ext - Rn) * 0.5 * (s1 - s2) / p.L / (1 - g0); },
    zmax: p.zt + 8 * p.L, ...p
  };
}
function dPlocal(pf, Rn, a, Gs, cv, N = 6000) {
  if (Rn >= a) return 0;
  let best = 0;
  for (let k = 1; k <= N; k++) {
    const z = pf.zmax * k / N, R = pf.Rw(z, Rn), u = R / a;
    if (R > a || u < cv.xmin) continue;
    const q = Gs * a * Math.abs(cv.df(u)) * Math.abs(pf.dRw(z, Rn)) / (Math.PI * R * R);
    if (q > best) best = q;
  }
  return best;
}
function dPwork(Rn, a, Gs, cv) { return Rn >= a ? 0 : Gs * a * cv.f(Rn / a) / (Math.PI * Rn * Rn); }
function dPdrop(Rn, a, gamma) { return (Rn >= a || !gamma) ? 0 : 2 * gamma * (1 / Rn - 1 / a); }
if (typeof module !== 'undefined') module.exports = { pchip, energyCurve, wallProfile, dPlocal, dPwork, dPdrop };
