/* SAXS-FFT core: scattering intensity, radial averaging and analytic curves.
 *
 * I(q) = |a^2 sum_x drho(x) w(x) exp(-i q.x)|^2 / (A <w^2>),  A = (N a)^2,
 * i.e. intensity per unit sample area in units of [SLD^2 nm^2].
 * With q = 2 pi k / (M a) on the (possibly zero-padded, M = pad N) grid.
 * Parseval: sum_q I(q) (dq / 2 pi)^2 = <drho^2>  (the invariant Q).
 */
(function (root) {
  root.SAXS_PARTS = root.SAXS_PARTS || [];
  root.SAXS_PARTS.push(function partAnalysis(S) {
    var PI = Math.PI;

    /* Bessel J1 (Numerical Recipes rational + asymptotic approximations). */
    function besselJ1(x) {
      var ax = Math.abs(x), y, a1, a2;
      if (ax < 8) {
        y = x * x;
        a1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
        a2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))));
        return a1 / a2;
      }
      var z = 8 / ax, xx = ax - 2.356194491;
      y = z * z;
      a1 = 1 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
      a2 = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
      var ans = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * a1 - z * Math.sin(xx) * a2);
      return x < 0 ? -ans : ans;
    }

    /* Normalised disk amplitude 2 J1(x) / x. */
    function diskAmp(x) {
      if (Math.abs(x) < 1e-4) return 1 - x * x / 8;
      return 2 * besselJ1(x) / x;
    }

    /* Orientation-averaged |F(q)|^2 for one particle shape at scale s. */
    function makeFormFn(form, rhoM, rhoP) {
      if (form.kind === 'disk') return function (q, s) {
        var r = form.R * s, F = (rhoP - rhoM) * PI * r * r * diskAmp(q * r);
        return F * F;
      };
      if (form.kind === 'cs') return function (q, s) {
        var rc = form.rc * s, ro = form.ro * s;
        var F = (form.rhoC - form.rhoS) * PI * rc * rc * diskAmp(q * rc) + (form.rhoS - rhoM) * PI * ro * ro * diskAmp(q * ro);
        return F * F;
      };
      if (form.kind === 'ell') {
        var NPSI = 48, cs2 = [], sn2 = [];
        for (var k = 0; k < NPSI; k++) {
          var psi = (k + 0.5) * (PI / 2) / NPSI;
          cs2.push(Math.cos(psi) * Math.cos(psi)); sn2.push(Math.sin(psi) * Math.sin(psi));
        }
        return function (q, s) {
          var ea = form.a * s, eb = form.b * s, V = (rhoP - rhoM) * PI * ea * eb, sum = 0;
          for (var k = 0; k < NPSI; k++) {
            var F = V * diskAmp(q * Math.sqrt(ea * ea * cs2[k] + eb * eb * sn2[k]));
            sum += F * F;
          }
          return sum / NPSI;
        };
      }
      return null;
    }

    /* Compress the particle size scales into (scale, weight) groups. */
    function scaleGroups(parts) {
      var ss = parts.map(function (p) { return p.s; }), uniq = {}, nu = 0;
      for (var i = 0; i < ss.length && nu <= 200; i++) if (!uniq[ss[i]]) { uniq[ss[i]] = 1; nu++; }
      var groups = [];
      if (nu <= 200) {
        var cnt = {};
        ss.forEach(function (s) { cnt[s] = (cnt[s] || 0) + 1; });
        for (var k in cnt) groups.push({ s: +k, w: cnt[k] });
        return groups;
      }
      var mn = Math.min.apply(null, ss), mx = Math.max.apply(null, ss), nb = 200, h = new Float64Array(nb), hs = new Float64Array(nb);
      ss.forEach(function (s) {
        var b = Math.min(nb - 1, Math.floor((s - mn) / (mx - mn) * nb));
        h[b]++; hs[b] += s;
      });
      for (var b = 0; b < nb; b++) if (h[b]) groups.push({ s: hs[b] / h[b], w: h[b] });
      return groups;
    }

    /* Dilute-limit expectation I(q) = (1/A) sum_i <|F_i(q)|^2>; returns a function of q. */
    function dilute(gen, P, A) {
      if (!gen.form || !gen.particles || !gen.particles.length) return null;
      var fn = makeFormFn(gen.form, P.rhoM, P.rhoP);
      if (!fn) return null;
      var groups = scaleGroups(gen.particles);
      return function (q) {
        var s = 0;
        for (var g = 0; g < groups.length; g++) s += groups[g].w * fn(q, groups[g].s);
        return s / A;
      };
    }

    /* Azimuthal average of the pixel transfer function [sinc(qx a/2) sinc(qy a/2)]^2:
       rasterising with area coverage convolves the structure with one pixel. */
    function pixelMTF(q, a) {
      var NA = 32, s = 0;
      function sinc(x) { return Math.abs(x) < 1e-8 ? 1 : Math.sin(x) / x; }
      for (var k = 0; k < NA; k++) {
        var ph = (k + 0.5) * (PI / 2) / NA, v = sinc(q * Math.cos(ph) * a / 2) * sinc(q * Math.sin(ph) * a / 2);
        s += v * v;
      }
      return s / NA;
    }

    /* Exponentially scaled modified Bessel function e^-x I0(x), x >= 0 (Numerical Recipes). */
    function besselI0e(x) {
      if (x < 3.75) {
        var y = (x / 3.75) * (x / 3.75);
        return Math.exp(-x) * (1 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2))))));
      }
      var z = 3.75 / x;
      return (0.39894228 + z * (0.1328592e-1 + z * (0.225319e-2 + z * (-0.157565e-2 + z * (0.916281e-2 + z * (-0.2057706e-1 +
             z * (0.2635537e-1 + z * (-0.1647633e-1 + z * 0.392377e-2)))))))) / Math.sqrt(x);
    }

    /* Resolution smearing: periodic separable Gaussian (sigma in grid pixels) on the
       unshifted M x M intensity. The kernel sums to 1, so the invariant is kept. */
    function smear2D(I, M, sig) {
      var R = Math.ceil(4 * sig), K = new Float64Array(2 * R + 1), ks = 0, i, x, y, j;
      for (i = -R; i <= R; i++) { K[i + R] = Math.exp(-0.5 * i * i / (sig * sig)); ks += K[i + R]; }
      for (i = 0; i <= 2 * R; i++) K[i] /= ks;
      var row = new Float64Array(M), tmp = new Float64Array(M * M);
      for (y = 0; y < M; y++) {
        for (x = 0; x < M; x++) row[x] = I[y * M + x];
        for (x = 0; x < M; x++) {
          var v = 0;
          for (j = -R; j <= R; j++) v += K[j + R] * row[(x + j + M * R) % M];
          tmp[y * M + x] = v;
        }
      }
      for (x = 0; x < M; x++) {
        for (y = 0; y < M; y++) row[y] = tmp[y * M + x];
        for (y = 0; y < M; y++) {
          var w = 0;
          for (j = -R; j <= R; j++) w += K[j + R] * row[(y + j + M * R) % M];
          I[y * M + x] = w;
        }
      }
    }

    /* The same 2D Gaussian smearing for an isotropic curve f(q): the radial kernel is
       (q'/s^2) exp(-(q-q')^2/2s^2) e^(-qq'/s^2) I0(qq'/s^2). f is tabulated once. */
    function smearRadial(f, qs, sq) {
      var h = sq / 6, qend = qs[qs.length - 1] + 7 * sq, nt = Math.ceil(qend / h) + 1, tab = new Float64Array(nt), t;
      for (t = 0; t < nt; t++) tab[t] = f(t * h);
      return qs.map(function (q) {
        var t0 = Math.max(0, Math.floor((q - 7 * sq) / h)), t1 = Math.min(nt - 1, Math.ceil((q + 7 * sq) / h)), num = 0, den = 0;
        for (var t = t0; t <= t1; t++) {
          var qp = t * h, w = qp * Math.exp(-0.5 * (q - qp) * (q - qp) / (sq * sq)) * besselI0e(q * qp / (sq * sq));
          num += w * tab[t]; den += w;
        }
        return den > 0 ? num / den : f(q);
      });
    }

    function hann(n) {
      var w = new Float64Array(n);
      for (var i = 0; i < n; i++) w[i] = 0.5 * (1 - Math.cos(2 * PI * (i + 0.5) / n));
      return w;
    }

    /* Main entry: rho = SLD map (N x N), O = scattering options. */
    function scatter(rho, gen, P, O) {
      var N = P.N, a = P.a, pad = O.pad || 1, M = N * pad, A = (N * a) * (N * a);
      var n2 = N * N, mean = 0, i, x, y;
      for (i = 0; i < n2; i++) mean += rho[i];
      mean /= n2;
      var variance = 0;
      for (i = 0; i < n2; i++) variance += (rho[i] - mean) * (rho[i] - mean);
      variance /= n2;

      var w = O.window === 'hann' ? hann(N) : null, w2 = 1;
      if (w) { w2 = 0; for (i = 0; i < N; i++) w2 += w[i] * w[i]; w2 = (w2 / N) * (w2 / N); }
      var re = new Float64Array(M * M), im = new Float64Array(M * M), sub = O.meanSub ? mean : 0;
      for (y = 0; y < N; y++) for (x = 0; x < N; x++)
        re[y * M + x] = (rho[y * N + x] - sub) * (w ? w[x] * w[y] : 1);
      S.fft2(re, im, M, false);

      var scale = a * a * a * a / (A * w2), I = new Float32Array(M * M);
      for (i = 0; i < M * M; i++) I[i] = (re[i] * re[i] + im[i] * im[i]) * scale;
      re = im = null;

      // instrumental resolution: Gaussian in q with FWHM = O.smear * (2 pi / L)
      var sigQ = O.smear > 0 ? O.smear * (2 * PI / (N * a)) / (2 * Math.sqrt(2 * Math.LN2)) : 0;
      if (sigQ > 0) smear2D(I, M, sigQ / (2 * PI / (M * a)));

      // radial average in bins of width dq centred on k = 0, 1, 2, ...
      var dq = 2 * PI / (M * a), kmax = O.range === 'corner' ? Math.ceil(M / Math.SQRT2) : M / 2;
      var sI = new Float64Array(kmax + 1), sQ = new Float64Array(kmax + 1), cnt = new Float64Array(kmax + 1), Q2D = 0;
      for (y = 0; y < M; y++) {
        var ky = y <= M / 2 ? y : y - M;
        for (x = 0; x < M; x++) {
          var kx = x <= M / 2 ? x : x - M, kr = Math.sqrt(kx * kx + ky * ky), b = Math.round(kr), v = I[y * M + x];
          Q2D += v;
          if (b > kmax) continue;
          sI[b] += v; sQ[b] += kr * dq; cnt[b]++;
        }
      }
      Q2D *= (dq / (2 * PI)) * (dq / (2 * PI));
      var qs = [], Is = [], ns = [];
      for (var b2 = 1; b2 <= kmax; b2++) if (cnt[b2]) { qs.push(sQ[b2] / cnt[b2]); Is.push(sI[b2] / cnt[b2]); ns.push(cnt[b2]); }
      var Q1D = 0;                                       // (1/2pi) int q I dq, trapezoid
      for (i = 1; i < qs.length; i++) Q1D += 0.5 * (qs[i] * Is[i] + qs[i - 1] * Is[i - 1]) * (qs[i] - qs[i - 1]);
      Q1D /= 2 * PI;

      // Porod: total variation sum |grad rho| a = |drho| * perimeter
      var tv = 0;
      for (y = 0; y < N; y++) for (x = 0; x < N; x++) {
        var gx = (rho[y * N + (x + 1) % N] - rho[y * N + (x + N - 1) % N]) / 2;
        var gy = (rho[((y + 1) % N) * N + x] - rho[((y + N - 1) % N) * N + x]) / 2;
        tv += Math.sqrt(gx * gx + gy * gy);
      }
      var contrast = Math.abs(P.rhoP - P.rhoM), perimeter = contrast > 0 ? tv * a / contrast : 0;
      var porodK = 2 * contrast * contrast * perimeter / A;    // I -> porodK / q^3

      var res = {
        I2D: S.fftshift(I, M, Float32Array), M: M, dq: dq,
        profile: { q: qs, I: Is, n: ns },
        stats: { mean: mean, variance: variance, Q2D: Q2D, Q1D: Q1D, perimeter: perimeter, porodK: porodK,
                 qmin: 2 * PI / (N * a), qmax: PI / a }
      };

      var dil = dilute(gen, P, A);
      if (dil) {
        var fq = function (q) { return dil(q) * (O.pixelMTF ? pixelMTF(q, a) : 1); };
        res.dilute = sigQ > 0 ? smearRadial(fq, qs, sigQ) : qs.map(fq);
        var I0 = dil(0), rb = gen.form.kind === 'disk' ? gen.form.R : gen.form.kind === 'cs' ? gen.form.ro : Math.max(gen.form.a, gen.form.b);
        var q1 = 0.02 / rb, rg2 = 2 * (1 - dil(q1) / I0) / (q1 * q1);
        res.stats.I0 = I0; res.stats.Rg = Math.sqrt(Math.max(rg2, 0)); res.stats.guinierOf = 'particles';
        res.stats.np = gen.particles.length;
      }
      if (gen.guinier !== undefined) {                   // generator knows better (clusters / none)
        res.stats.I0 = gen.guinier ? gen.guinier.I0 : 0;
        res.stats.Rg = gen.guinier ? gen.guinier.Rg : 0;
        res.stats.guinierOf = gen.guinier ? gen.guinier.label : '';
      }
      if (w) res.stats.windowNote = true;
      if (sigQ > 0) res.stats.smearFWHM = sigQ * 2 * Math.sqrt(2 * Math.LN2);
      return res;
    }

    S.scatter = scatter;
    S.besselJ1 = besselJ1;
    S.diskAmp = diskAmp;
  });
})(typeof self !== 'undefined' ? self : globalThis);
