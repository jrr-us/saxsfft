/* SAXS-FFT core: antialiased rasterisation of particles into an SLD map.
 *
 * Pixel (i, j) covers [i a, (i+1) a) x [j a, (j+1) a). Pixels fully inside or
 * outside a shape get coverage 1 / 0; pixels near an edge are 4x4
 * supersampled. Shapes wrap around the periodic box. Shapes are composited
 * in order: rho = rho * (1 - cov) + rho_shape * cov.
 */
(function (root) {
  root.SAXS_PARTS = root.SAXS_PARTS || [];
  root.SAXS_PARTS.push(function partRaster(S) {
    var SUB = 4;

    /* Composite one ellipse (semi-axes ea, eb, angle th) centred at (cx, cy). */
    function drawEllipse(rho, N, a, cx, cy, ea, eb, th, val) {
      var c = Math.cos(th), s = Math.sin(th);
      var hx = Math.sqrt(ea * ea * c * c + eb * eb * s * s), hy = Math.sqrt(ea * ea * s * s + eb * eb * c * c);
      var x0 = Math.floor((cx - hx) / a), x1 = Math.floor((cx + hx) / a);
      var y0 = Math.floor((cy - hy) / a), y1 = Math.floor((cy + hy) / a);
      if (x1 - x0 >= N) { x0 = Math.floor(cx / a) - (N >> 1); x1 = x0 + N - 1; }
      if (y1 - y0 >= N) { y0 = Math.floor(cy / a) - (N >> 1); y1 = y0 + N - 1; }
      // Pixels within the shrunk copy (1 - m/bmin) E are fully inside and pixels
      // outside the grown copy (1 + m/bmin) E fully outside (E convex, contains
      // a disk of radius bmin), with margin m >= half the pixel diagonal.
      var m = 0.75 * a, bmin = Math.min(ea, eb), kin = 1 - m / bmin, kout = 1 + m / bmin;
      var ai2 = 1 / (ea * ea), bi2 = 1 / (eb * eb), hasIn = kin > 0;
      var kin2 = hasIn ? kin * kin : 0, kout2 = kout * kout;
      for (var j = y0; j <= y1; j++) {
        var jj = ((j % N) + N) % N, py = (j + 0.5) * a - cy;
        for (var i = x0; i <= x1; i++) {
          var px = (i + 0.5) * a - cx;
          var u = c * px + s * py, v = -s * px + c * py, cov;
          var e = u * u * ai2 + v * v * bi2;
          if (e > kout2) continue;
          if (hasIn && e <= kin2) cov = 1;
          else {
            var hit = 0;
            for (var sy = 0; sy < SUB; sy++) {
              var qy = (j + (sy + 0.5) / SUB) * a - cy;
              for (var sx = 0; sx < SUB; sx++) {
                var qx = (i + (sx + 0.5) / SUB) * a - cx;
                var uu = c * qx + s * qy, vv = -s * qx + c * qy;
                if (uu * uu * ai2 + vv * vv * bi2 <= 1) hit++;
              }
            }
            if (!hit) continue;
            cov = hit / (SUB * SUB);
          }
          var k = jj * N + (((i % N) + N) % N);
          rho[k] += (val - rho[k]) * cov;
        }
      }
    }

    function renderParticles(parts, N, a, rhoM) {
      var rho = new Float64Array(N * N).fill(rhoM);
      for (var n = 0; n < parts.length; n++) {
        var p = parts[n];
        if (p.kind === 'cs') {
          drawEllipse(rho, N, a, p.x, p.y, p.a, p.a, 0, p.rhoS);
          drawEllipse(rho, N, a, p.x, p.y, p.rc, p.rc, 0, p.rhoC);
        } else drawEllipse(rho, N, a, p.x, p.y, p.a, p.b, p.th, p.rho);
      }
      return rho;
    }

    /* Coverage field (0..1) -> SLD map. */
    function renderField(cov, rhoM, rhoP) {
      var rho = new Float64Array(cov.length);
      for (var i = 0; i < cov.length; i++) rho[i] = rhoM + (rhoP - rhoM) * cov[i];
      return rho;
    }

    /* Build the SLD map from a generator result. invert = Babinet swap. */
    function rasterize(gen, P) {
      var rho = gen.particles ? renderParticles(gen.particles, P.N, P.a, P.rhoM)
                              : renderField(gen.field, P.rhoM, P.rhoP);
      if (P.invert) {
        var s = P.rhoM + P.rhoP;
        for (var i = 0; i < rho.length; i++) rho[i] = s - rho[i];
      }
      return rho;
    }

    S.rasterize = rasterize;
    S.renderParticles = renderParticles;
  });
})(typeof self !== 'undefined' ? self : globalThis);
