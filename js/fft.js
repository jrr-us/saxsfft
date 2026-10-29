/* SAXS-FFT core: radix-2 complex FFT (1D and 2D) and fftshift.
 *
 * Every core file registers a self-contained "part" function into
 * SAXS_PARTS. The parts are executed once on the main thread and their
 * source is also stitched into a Blob Web Worker (works from file://).
 * Parts must therefore not reference anything outside their argument S.
 */
(function (root) {
  root.SAXS_PARTS = root.SAXS_PARTS || [];
  root.SAXS_PARTS.push(function partFFT(S) {
    var twiddleCache = {};

    function twiddles(n) {
      if (twiddleCache[n]) return twiddleCache[n];
      var c = new Float64Array(n / 2), s = new Float64Array(n / 2);
      for (var k = 0; k < n / 2; k++) {
        c[k] = Math.cos(-2 * Math.PI * k / n);
        s[k] = Math.sin(-2 * Math.PI * k / n);
      }
      return (twiddleCache[n] = { c: c, s: s });
    }

    function isPow2(n) { return n > 0 && (n & (n - 1)) === 0; }

    /* In-place iterative FFT of (re, im), length n = power of two.
       inverse=true computes the unnormalised inverse transform. */
    function fft1(re, im, inverse) {
      var n = re.length, i, j, k, len, half, step;
      if (!isPow2(n)) throw new Error('FFT length must be a power of two');
      for (i = 1, j = 0; i < n; i++) {           // bit reversal
        var bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
          var t = re[i]; re[i] = re[j]; re[j] = t;
          t = im[i]; im[i] = im[j]; im[j] = t;
        }
      }
      var tw = twiddles(n), sgn = inverse ? -1 : 1;
      for (len = 2; len <= n; len <<= 1) {
        half = len >> 1; step = n / len;
        for (i = 0; i < n; i += len) {
          for (k = 0; k < half; k++) {
            var wr = tw.c[k * step], wi = sgn * tw.s[k * step];
            var a = i + k, b = a + half;
            var xr = re[b] * wr - im[b] * wi;
            var xi = re[b] * wi + im[b] * wr;
            re[b] = re[a] - xr; im[b] = im[a] - xi;
            re[a] += xr; im[a] += xi;
          }
        }
      }
    }

    /* In-place 2D FFT of an n x n row-major complex array. */
    function fft2(re, im, n, inverse) {
      var rr = new Float64Array(n), ri = new Float64Array(n), x, y, o;
      for (y = 0; y < n; y++) {
        o = y * n;
        for (x = 0; x < n; x++) { rr[x] = re[o + x]; ri[x] = im[o + x]; }
        fft1(rr, ri, inverse);
        for (x = 0; x < n; x++) { re[o + x] = rr[x]; im[o + x] = ri[x]; }
      }
      for (x = 0; x < n; x++) {
        for (y = 0; y < n; y++) { rr[y] = re[y * n + x]; ri[y] = im[y * n + x]; }
        fft1(rr, ri, inverse);
        for (y = 0; y < n; y++) { re[y * n + x] = rr[y]; im[y * n + x] = ri[y]; }
      }
    }

    /* Return a copy of an n x n array with q = 0 moved to (n/2, n/2). */
    function fftshift(a, n, Ctor) {
      var out = new (Ctor || a.constructor)(n * n), h = n >> 1;
      for (var y = 0; y < n; y++) {
        var ys = ((y + h) % n) * n, o = y * n;
        for (var x = 0; x < n; x++) out[ys + ((x + h) % n)] = a[o + x];
      }
      return out;
    }

    S.fft1 = fft1;
    S.fft2 = fft2;
    S.fftshift = fftshift;
    S.isPow2 = isPow2;
  });
})(typeof self !== 'undefined' ? self : globalThis);
