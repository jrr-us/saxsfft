/* SAXS-FFT core: seeded random numbers, so every figure is reproducible. */
(function (root) {
  root.SAXS_PARTS = root.SAXS_PARTS || [];
  root.SAXS_PARTS.push(function partRNG(S) {
    /* mulberry32 generator; returns an object with uniform/normal helpers. */
    function makeRng(seed) {
      var a = (seed >>> 0) || 1, spare = null;
      function uniform() {
        a = (a + 0x6D2B79F5) | 0;
        var t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
      }
      function normal() {
        if (spare !== null) { var s = spare; spare = null; return s; }
        var u, v, r;
        do { u = 2 * uniform() - 1; v = 2 * uniform() - 1; r = u * u + v * v; }
        while (r >= 1 || r === 0);
        var f = Math.sqrt(-2 * Math.log(r) / r);
        spare = v * f;
        return u * f;
      }
      return { uniform: uniform, normal: normal };
    }

    /* Draw a size scale factor with mean 1 and relative std sigma.
       dist: 'gaussian' (truncated at 0.05) or 'lognormal'. */
    function sizeScale(rng, sigma, dist) {
      if (!(sigma > 0)) return 1;
      if (dist === 'lognormal') {
        var s2 = Math.log(1 + sigma * sigma);
        return Math.exp(-s2 / 2 + Math.sqrt(s2) * rng.normal());
      }
      var x;
      do { x = 1 + sigma * rng.normal(); } while (x < 0.05);
      return x;
    }

    S.makeRng = makeRng;
    S.sizeScale = sizeScale;
  });
})(typeof self !== 'undefined' ? self : globalThis);
