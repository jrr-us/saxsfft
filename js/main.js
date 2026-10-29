/* SAXS-FFT UI: parameter panels, pipeline orchestration, drawing, export. */
(function () {
  'use strict';
  var PI = Math.PI, fmt = SAXS_plot.fmt;
  function $(id) { return document.getElementById(id); }

  /* ---------- parameter schema ---------- */

  var FIELDS = {
    R: { label: 'Radius R', unit: 'nm', min: 0.5, step: 0.5 },
    sigma: { label: 'Polydispersity σ/R', min: 0, max: 0.6, step: 0.01 },
    dist: { label: 'Size distribution', type: 'select', options: [['gaussian', 'Gaussian'], ['lognormal', 'log-normal']], show: function (v) { return v.sigma > 0; } },
    amount: { label: 'Amount given by', type: 'select', options: [['phi', 'area fraction φ'], ['count', 'number of particles']] },
    phi: { label: 'Area fraction φ', min: 0.0001, max: 0.85, step: 0.01, show: function (v) { return v.amount !== 'count'; } },
    count: { label: 'Number', min: 1, step: 1, show: function (v) { return v.amount === 'count'; } },
    packing: { label: 'Positions', type: 'select', options: [['rsa', 'random, no overlap (RSA)'], ['dense', 'hard-particle fluid (MC)'], ['overlap', 'random, overlaps allowed']] },
    sweeps: { label: 'MC sweeps', min: 0, step: 50, show: function (v) { return v.packing === 'dense'; } },
    ea: { label: 'Semi-axis a', unit: 'nm', min: 0.5, step: 0.5 },
    eb: { label: 'Semi-axis b', unit: 'nm', min: 0.5, step: 0.5 },
    orient: { label: 'Orientation', type: 'select', options: [['random', 'random'], ['aligned', 'aligned']] },
    theta: { label: 'Angle θ', unit: '°', step: 5, show: function (v) { return v.orient === 'aligned'; } },
    sigmaTheta: { label: 'Angular spread σθ', unit: '°', min: 0, step: 1, show: function (v) { return v.orient === 'aligned'; } },
    Rc: { label: 'Core radius', unit: 'nm', min: 0.5, step: 0.5 },
    t: { label: 'Shell thickness', unit: 'nm', min: 0, step: 0.5 },
    rhoC: { label: 'ρ core', step: 0.1 },
    rhoS: { label: 'ρ shell', step: 0.1 },
    filter: { label: 'Spectrum', type: 'select', options: [['ring', 'ring at q₀ = 2π/d₀ (spinodal)'], ['gauss', 'Gaussian (random blobs)']] },
    d0: { label: 'Length scale d₀', unit: 'nm', min: 1, step: 1 },
    width: { label: 'Ring width Δq/q₀', min: 0.02, max: 2, step: 0.02, show: function (v) { return v.filter === 'ring'; } },
    r0: { label: 'Particle radius', unit: 'nm', min: 0.5, step: 0.5 },
    M: { label: 'Particles / cluster', min: 2, step: 100 },
    nClusters: { label: 'Clusters', min: 1, step: 1 }
  };

  var STRUCTURES = {
    single: { label: 'Single disk', gen: 'single', fields: { R: 20 },
      hint: 'I(q) = (Δρ πR²)² [2J₁(qR)/(qR)]² / A. Minima at qR = 3.83, 7.02, 10.17; Guinier Rg² = R²/2; Porod I ∝ q⁻³ in 2D.' },
    dilute: { label: 'Dilute disks', gen: 'disks', fields: { R: 10, sigma: 0, dist: 'gaussian', amount: 'count', phi: 0.02, count: 30, packing: 'rsa', sweeps: 300 },
      hint: 'Few, uncorrelated disks: I(q) ≈ n·P(q)/A, S(q) ≈ 1. Add polydispersity to wash out the minima.' },
    dense: { label: 'Hard-disk fluid (dense)', gen: 'disks', fields: { R: 8, sigma: 0, dist: 'gaussian', amount: 'phi', phi: 0.5, count: 500, packing: 'dense', sweeps: 300 },
      hint: 'Excluded volume creates a correlation peak near q ≈ 2π/(2R) and suppresses I(q→0): S(0) < 1.' },
    ellipses: { label: 'Ellipses', gen: 'ellipses', fields: { ea: 20, eb: 7, sigma: 0, dist: 'gaussian', amount: 'phi', phi: 0.1, count: 50, packing: 'rsa', sweeps: 200, orient: 'aligned', theta: 30, sigmaTheta: 5 },
      hint: 'Aligned ellipses give an anisotropic pattern, elongated perpendicular to the long axis. Random orientation gives an isotropic one.' },
    coreshell: { label: 'Core–shell disks', gen: 'coreshell', fields: { Rc: 12, t: 5, rhoC: 1, rhoS: -0.5, sigma: 0.03, dist: 'gaussian', amount: 'phi', phi: 0.05, count: 50, packing: 'rsa', sweeps: 300 },
      hint: 'Two interfaces interfere: F = (ρc−ρs)πRc²f(qRc) + (ρs−ρm)πRo²f(qRo). Contrast matching (ρs = ρm or ρc = ρs) isolates one of them.' },
    spinodal: { label: 'Two-phase random (spinodal)', gen: 'spinodal', fields: { filter: 'ring', d0: 60, width: 0.3, phi: 0.5 },
      hint: 'Bicontinuous structure with a characteristic length d₀: broad peak at q* ≈ 2π/d₀, Porod q⁻³ at high q.' },
    dla: { label: 'Fractal aggregates (DLA)', gen: 'dla', minN: 1024, fields: { r0: 2, M: 1500, nClusters: 2 },
      hint: 'Mass fractal: I ∝ q^−Df between 2π/Rg and 2π/r₀ (Df ≈ 1.71 for 2D DLA). Small clusters show a smaller apparent slope.' }
  };

  var state = {
    struct: 'dilute',
    values: {},               // per-structure field values
    data: null,               // last structure: {rho, N, a, L, info}
    scat: null,               // last scattering result
    reqId: 0, busy: false, needGen: true, tab: 'I', lastP: null
  };
  Object.keys(STRUCTURES).forEach(function (k) { state.values[k] = Object.assign({}, STRUCTURES[k].fields); });

  /* ---------- building the controls ---------- */

  var selStruct = $('structure');
  Object.keys(STRUCTURES).forEach(function (k) {
    var o = document.createElement('option'); o.value = k; o.textContent = STRUCTURES[k].label; selStruct.appendChild(o);
  });
  selStruct.value = state.struct;
  SAXS_colormaps.names.forEach(function (n) {
    var o = document.createElement('option'); o.value = n; o.textContent = n; $('cmap').appendChild(o);
  });

  function buildFields() {
    var box = $('structFields'), def = STRUCTURES[state.struct], vals = state.values[state.struct];
    box.innerHTML = '';
    $('hint').textContent = def.hint;
    Object.keys(def.fields).forEach(function (key) {
      var f = FIELDS[key], el, wrap = document.createElement('label');
      if (f.type === 'check') {
        wrap.className = 'check';
        el = document.createElement('input'); el.type = 'checkbox'; el.checked = !!vals[key];
        wrap.appendChild(el); wrap.appendChild(document.createTextNode(' ' + f.label));
      } else {
        wrap.className = 'row';
        var sp = document.createElement('span'); sp.textContent = f.label; wrap.appendChild(sp);
        if (f.type === 'select') {
          el = document.createElement('select');
          f.options.forEach(function (op) { var o = document.createElement('option'); o.value = op[0]; o.textContent = op[1]; el.appendChild(o); });
          el.value = vals[key];
        } else {
          el = document.createElement('input'); el.type = 'number'; el.value = vals[key];
          if (f.min !== undefined) el.min = f.min;
          if (f.max !== undefined) el.max = f.max;
          el.step = f.step || 'any';
        }
        wrap.appendChild(el);
        if (f.unit) { var em = document.createElement('em'); em.textContent = f.unit; wrap.appendChild(em); }
      }
      wrap.dataset.key = key;
      el.addEventListener('change', function () {
        vals[key] = f.type === 'check' ? el.checked : f.type === 'select' ? el.value : parseFloat(el.value);
        updateVisibility();
        structureChanged();
      });
      box.appendChild(wrap);
    });
    $('rhoPRow').style.display = state.struct === 'coreshell' ? 'none' : '';
    updateVisibility();
  }

  function updateVisibility() {
    var vals = state.values[state.struct];
    Array.prototype.forEach.call($('structFields').children, function (w) {
      var f = FIELDS[w.dataset.key];
      w.style.display = !f.show || f.show(vals) ? '' : 'none';
    });
  }

  /* ---------- collecting parameters ---------- */

  function num(id) { return parseFloat($(id).value); }

  function genParams() {
    var vals = state.values[state.struct], P = Object.assign({}, vals);
    P.gen = STRUCTURES[state.struct].gen;
    P.N = parseInt($('N').value, 10);
    P.a = num('a');
    P.seed = parseInt($('seed').value, 10) || 1;
    P.rhoM = num('rhoM');
    P.rhoP = state.struct === 'coreshell' ? vals.rhoS : num('rhoP');
    P.invert = $('invert').checked;
    if (P.amount === 'count') P.count = Math.max(1, Math.round(P.count));
    if (P.M) P.M = Math.max(2, Math.round(P.M));
    if (P.nClusters) P.nClusters = Math.max(1, Math.round(P.nClusters));
    return P;
  }

  function scatParams() {
    return { meanSub: $('meanSub').checked, window: $('window').value, pad: parseInt($('pad').value, 10),
             range: $('range').value, pixelMTF: $('pixelMTF').checked,
             smear: $('smear').checked ? Math.max(0, parseFloat($('smearW').value) || 0) : 0 };
  }

  function validate(P) {
    var bad = [];
    Object.keys(STRUCTURES[state.struct].fields).forEach(function (k) {
      var v = P[k], f = FIELDS[k];
      if (f.type || (f.show && !f.show(P))) return;
      if (!isFinite(v)) bad.push(f.label + ' is not a number');
      else if (f.min !== undefined && v < f.min) bad.push(f.label + ' must be ≥ ' + f.min);
      else if (f.max !== undefined && v > f.max) bad.push(f.label + ' must be ≤ ' + f.max);
    });
    if (!(P.a > 0)) bad.push('Pixel size must be > 0');
    if (P.packing === 'dense' && P.amount !== 'count' && P.phi > 0.72) bad.push('A hard-particle fluid needs φ ≤ 0.72');
    return bad;
  }

  function updateGridInfo() {
    var N = parseInt($('N').value, 10), a = num('a'), L = N * a;
    $('gridInfo').textContent = 'Box L = ' + fmt(L) + ' nm · q from ' + fmt(2 * PI / L) + ' to ' + fmt(PI / a) + ' nm⁻¹ (Δq = 2π/L)';
  }

  /* ---------- engine ---------- */

  var engine = SAXS_createEngine(onReply);
  var genTimer = null, scatTimer = null;

  function setStatus(text, cls) { var s = $('status'); s.textContent = text; s.className = 'status ' + (cls || ''); }
  function setBusy(b) { state.busy = b; $('stop').disabled = !b || !engine.usesWorker(); }

  function request(withGen) {
    var P = genParams(), bad = validate(P);
    if (bad.length) { setStatus(bad.join('; '), 'error'); return; }
    if (state.busy && engine.restart()) state.needGen = true;      // abort the running job
    if (state.needGen) withGen = true;
    var msg = { id: ++state.reqId, gen: withGen ? P : null, scat: scatParams() };
    if (withGen) state.pendingP = P;
    setBusy(true);
    setStatus(withGen ? 'Generating…' : 'Computing FFT…', 'busy');
    engine.send(msg);
  }

  function onReply(m) {
    if (m.id !== state.reqId && m.id !== -1) return;              // stale
    if (m.progress) { setStatus(m.progress, 'busy'); return; }
    setBusy(false);
    if (!m.ok) { setStatus('Error: ' + m.error, 'error'); state.needGen = true; return; }
    if (m.rho) {
      var P = state.pendingP;
      state.data = { rho: m.rho, N: m.N, a: P.a, L: P.N * P.a, info: m.info, P: P };
      state.needGen = false;
      drawStructure();
    }
    state.scat = m.scat;
    state.scat.timings = m.timings;
    draw2D(true);
    draw1D(false);
    showInfo();
    setStatus('Done' + (m.timings.generate !== undefined ? ' · structure ' + m.timings.generate + ' ms' : '') + ' · scattering ' + m.timings.scatter + ' ms');
  }

  function structureChanged() {
    state.needGen = true;
    if (!$('auto').checked) { setStatus('Parameters changed: press Generate'); return; }
    clearTimeout(genTimer); genTimer = setTimeout(function () { request(true); }, 350);
  }
  function scatteringChanged() {
    if (!$('auto').checked) { setStatus('Parameters changed: press Generate'); return; }
    clearTimeout(scatTimer); scatTimer = setTimeout(function () { request(false); }, 150);
  }

  /* ---------- drawing ---------- */

  var structView = new SAXS_plot.ImageView($('cvStruct'), {
    xlabel: 'x (nm)', ylabel: 'y (nm)',
    onHover: function (x, y) {
      var d = state.data;
      if (x === null || !d) { $('roStruct').innerHTML = '&nbsp;'; return; }
      var i = Math.floor(x / d.a), j = Math.floor(y / d.a);
      var v = i >= 0 && j >= 0 && i < d.N && j < d.N ? d.rho[j * d.N + i] : NaN;
      $('roStruct').textContent = 'x = ' + x.toFixed(1) + ' nm, y = ' + y.toFixed(1) + ' nm, ρ = ' + (isFinite(v) ? v.toFixed(3) : '–');
    },
    overlay: function (ctx, F, V, toPx) {                         // scale bar
      var span = V.x1 - V.x0, nice = [1, 2, 5], len = 1;
      for (var e = -2; e < 7; e++) for (var k = 0; k < 3; k++) { var v = nice[k] * Math.pow(10, e); if (v <= span / 4) len = v; }
      var px = len / span * F.w, x = F.l + F.w - px - 12, y = F.t + F.h - 14;
      ctx.fillStyle = 'rgba(0,0,0,.75)'; ctx.fillRect(x - 6, y - 16, px + 12, 24);
      ctx.fillStyle = '#fff'; ctx.fillRect(x, y, px, 3);
      ctx.font = '11px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
      ctx.fillText(fmt(len) + ' nm', x + px / 2, y - 1);
    }
  });

  var view2D = new SAXS_plot.ImageView($('cv2D'), {
    xlabel: 'qx (nm⁻¹)', ylabel: 'qy (nm⁻¹)',
    onHover: function (qx, qy) {
      var s = state.scat;
      if (qx === null || !s) { $('ro2D').innerHTML = '&nbsp;'; return; }
      var M = s.M, i = Math.round(qx / s.dq) + M / 2, j = Math.round(qy / s.dq) + M / 2;
      var v = i >= 0 && j >= 0 && i < M && j < M ? s.I2D[j * M + i] : NaN, q = Math.hypot(qx, qy);
      $('ro2D').textContent = 'qx = ' + qx.toFixed(4) + ', qy = ' + qy.toFixed(4) + ', |q| = ' + q.toFixed(4) + ' nm⁻¹ (d = 2π/q = ' + (2 * PI / q).toFixed(1) + ' nm), I = ' + fmt(v);
    }
  });

  var plot1D = new SAXS_plot.LinePlot($('cv1D'), {
    xlog: true, ylog: true, maxDecades: 10, xlabel: 'q (nm⁻¹)', ylabel: 'I(q)  [Δρ² nm²]',
    onHover: function (q, series) {
      if (q === null || !series || !series.length) { $('ro1D').innerHTML = '&nbsp;'; return; }
      var s = series[0], best = 0;
      for (var i = 1; i < s.x.length; i++) if (Math.abs(s.x[i] - q) < Math.abs(s.x[best] - q)) best = i;
      var txt = 'q = ' + s.x[best].toFixed(4) + ' nm⁻¹, d = 2π/q = ' + (2 * PI / s.x[best]).toFixed(1) + ' nm';
      series.forEach(function (ss) { if (ss.label && ss.x === s.x) txt += ', ' + ss.label.split(' (')[0] + ' = ' + fmt(ss.y[best]); });
      $('ro1D').textContent = txt;
    }
  });

  function drawStructure() {
    var d = state.data;
    if (!d) return;
    var N = d.N, lo = Infinity, hi = -Infinity, i;
    for (i = 0; i < d.rho.length; i++) { if (d.rho[i] < lo) lo = d.rho[i]; if (d.rho[i] > hi) hi = d.rho[i]; }
    if (hi - lo < 1e-12) { lo -= 1; hi += 1; }
    var flipped = new Float32Array(N * N);                        // row 0 at the top = largest y
    for (var y = 0; y < N; y++) flipped.set(d.rho.subarray((N - 1 - y) * N, (N - y) * N), y * N);
    var off = document.createElement('canvas'); off.width = N; off.height = N;
    off.getContext('2d').putImageData(SAXS_colormaps.toImageData(flipped, N, N, lo, hi, 'gray'), 0, 0);
    structView.setImage(off, { x0: 0, x1: d.L, y0: 0, y1: d.L }, true);
  }

  function draw2D(newData) {
    var s = state.scat;
    if (!s) return;
    var M = s.M, n = M * M, c = (M / 2) * M + M / 2, lin = $('scale2D').value === 'lin';
    if (newData) s.img = {};
    var key = lin ? 'lin' : 'log';
    if (!s.img[key]) {
      var img = new Float32Array(n), mx = -Infinity;
      for (var y = 0; y < M; y++) {                    // flip so +qy points up
        var src = (M - 1 - y) * M, dst = y * M;
        for (var x = 0; x < M; x++) {
          var v = s.I2D[src + x], iv = lin ? v : (v > 0 ? Math.log10(v) : -300);
          img[dst + x] = iv;
          if (src + x !== c && iv > mx) mx = iv;
        }
      }
      // Max-pooled pyramid: isolated Bragg peaks stay visible when downscaled.
      var pyr = [{ v: img, n: M }];
      for (var lv = img, n2 = M; n2 > 128; ) {
        var h = n2 >> 1, pooled = new Float32Array(h * h);
        for (var py = 0; py < h; py++) for (var px = 0; px < h; px++) {
          var o = 2 * py * n2 + 2 * px;
          pooled[py * h + px] = Math.max(lv[o], lv[o + 1], lv[o + n2], lv[o + n2 + 1]);
        }
        pyr.push({ v: pooled, n: h });
        lv = pooled; n2 = h;
      }
      s.img[key] = { pyr: pyr, max: mx };
    }
    // log: 'top' lowers the ceiling by decades, 'decades' sets the span;
    // linear: 0 to max x 10^top (lowering it saturates the strong centre)
    var im = s.img[key], cm = $('cmap').value, hi, lo;
    if (lin) { hi = im.max * Math.pow(10, num('top')); lo = 0; }
    else { hi = im.max + num('top'); lo = hi - num('decades'); }
    var levels = im.pyr.map(function (L) {
      var off = document.createElement('canvas'); off.width = L.n; off.height = L.n;
      off.getContext('2d').putImageData(SAXS_colormaps.toImageData(L.v, L.n, L.n, lo, hi, cm), 0, 0);
      return off;
    });
    var e0 = (-M / 2 - 0.5) * s.dq, e1 = (M / 2 - 0.5) * s.dq;
    view2D.setImage(levels, { x0: e0, x1: e1, y0: e0, y1: e1 }, true);   // zoom kept unless the q grid changed
  }

  function sOfQ(s) {
    if (!s.dilute) return null;
    var D = s.dilute, I = s.profile.I, out = [];
    for (var i = 0; i < D.length; i++) {               // mask near the zeros of P(q)
      var m = 0;
      for (var k = Math.max(0, i - 8); k <= Math.min(D.length - 1, i + 8); k++) m = Math.max(m, D[k]);
      out.push(D[i] > 0.1 * m ? I[i] / D[i] : NaN);
    }
    return out;
  }

  function draw1D(keepView) {
    var s = state.scat;
    if (!s) return;
    var q = s.profile.q, series = [], st = s.stats;
    var c1 = cssVar('--series1'), c2 = cssVar('--series2'), c3 = cssVar('--series3'), c4 = cssVar('--series4');
    if (state.tab === 'I') {
      plot1D.o.ylabel = 'I(q)  [Δρ² nm²]'; plot1D.o.maxDecades = 10;
      series.push({ x: q, y: s.profile.I, color: c1, style: 'both', r: 1.6, width: 1, label: 'I(q) (FFT)' });
      if ($('ovDilute').checked && s.dilute) series.push({ x: q, y: s.dilute, color: c2, width: 1.6, dash: [6, 4], label: 'n⟨|F|²⟩/A (S = 1)' });
      if ($('ovPorod').checked && st.porodK > 0)
        series.push({ x: q, y: q.map(function (v) { return st.porodK / (v * v * v); }), color: c4, width: 1.4, dash: [2, 3], noAuto: true, label: 'Porod 2Δρ²(P/A)q⁻³' });
      if ($('ovGuinier').checked && st.Rg > 0) {
        var gq = [];                                     // own grid: big clusters leave few data bins
        for (var gv = q[0]; gv * st.Rg < 1.6; gv *= 1.05) gq.push(gv);
        series.push({ x: gq, y: gq.map(function (v) { return st.I0 * Math.exp(-v * v * st.Rg * st.Rg / 2); }), color: c3, width: 2, noAuto: true, label: 'Guinier, ' + st.guinierOf + ' (Rg = ' + st.Rg.toFixed(2) + ' nm)' });
      }
      plot1D.xlog = $('xlog').checked; plot1D.ylog = $('ylog').checked;
    } else {
      plot1D.o.ylabel = 'S(q) = I / (n⟨|F|²⟩/A)'; plot1D.o.maxDecades = 0;
      var S = sOfQ(s);
      if (S) {
        series.push({ x: q, y: S, color: c1, style: 'both', r: 1.6, width: 1, label: 'S(q)' });
        series.push({ x: [q[0], q[q.length - 1]], y: [1, 1], color: c4, width: 1, dash: [4, 4], noAuto: true });
      }
      plot1D.xlog = $('xlog').checked; plot1D.ylog = false;
    }
    plot1D.setData(series, keepView);
  }

  function cssVar(n) { return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }

  function showInfo() {
    var ul = $('info'), d = state.data, s = state.scat, lines = [];
    if (d) d.info.forEach(function (t) { lines.push(t); });
    if (s) {
      var st = s.stats, P = d ? d.P : null;
      lines.push('Box ' + fmt(d.L) + ' nm, ' + d.N + '² px · FFT grid ' + s.M + '² · q = ' + fmt(st.qmin) + ' … ' + fmt(st.qmax) + ' nm⁻¹');
      lines.push('Mean ρ = ' + st.mean.toFixed(4) + ' · ⟨Δρ²⟩ = ' + fmt(st.variance) + (P && !P.invert && P.gen !== 'coreshell' && P.rhoP !== P.rhoM ?
        '  (two-phase: φ(1−φ)Δρ² with φ = ' + phiFromMean(st.mean, P).toFixed(4) + ')' : ''));
      lines.push('Invariant: Σ I (Δq/2π)² = ' + fmt(st.Q2D) + ' (= ⟨Δρ²⟩, Parseval) · (1/2π)∫ q I dq to q_max = ' + fmt(st.Q1D) +
                 ' (' + (100 * st.Q1D / st.variance).toFixed(1) + '%)');
      if (st.perimeter > 0) lines.push('Interface length from ∑|∇ρ|: P = ' + fmt(st.perimeter) + ' nm, P/A = ' + fmt(st.perimeter / (d.L * d.L)) + ' nm⁻¹');
      if (st.Rg > 0) lines.push('Guinier (' + st.guinierOf + '): Rg = ' + st.Rg.toFixed(2) + ' nm, I(0) = ' + fmt(st.I0));
      if (st.smearFWHM) lines.push('Resolution smearing: Gaussian in q, FWHM = ' + fmt(st.smearFWHM) + ' nm⁻¹ (applied to the 2D data and to the analytic curve)');
      if (st.windowNote) lines.push('⚠ Hann window: I divided by ⟨w²⟩, which is right for homogeneous structures but not for a single centred particle');
      if (!engine.usesWorker()) lines.push('⚠ Web Worker unavailable: computing on the main thread');
    }
    ul.innerHTML = '';
    lines.forEach(function (t) {
      var li = document.createElement('li'); li.textContent = t;
      if (t.indexOf('⚠') === 0) li.className = 'warn';
      ul.appendChild(li);
    });
  }

  function phiFromMean(mean, P) { return (mean - P.rhoM) / (P.rhoP - P.rhoM); }

  /* ---------- export ---------- */

  function download(blob, name) {
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
  }
  function baseName() { return 'saxs_' + state.struct + '_N' + $('N').value + '_seed' + $('seed').value; }
  function exportPNG(canvas, suffix) { canvas.toBlob(function (b) { download(b, baseName() + '_' + suffix + '.png'); }); }

  function allParams() {
    var scal = {};
    ['N', 'a', 'seed', 'rhoM', 'rhoP', 'window', 'pad', 'range', 'smearW', 'cmap', 'scale2D', 'decades', 'top'].forEach(function (id) { scal[id] = $(id).value; });
    ['invert', 'meanSub', 'pixelMTF', 'smear', 'ovDilute', 'ovPorod', 'ovGuinier', 'xlog', 'ylog'].forEach(function (id) { scal[id] = $(id).checked; });
    return { app: 'SAXS-FFT 2D', version: 1, structure: state.struct, fields: state.values[state.struct], controls: scal };
  }

  function applyParams(p) {
    if (!p || !STRUCTURES[p.structure]) throw new Error('Not a SAXS-FFT parameter file');
    state.struct = p.structure; selStruct.value = p.structure;
    state.values[p.structure] = Object.assign({}, STRUCTURES[p.structure].fields, p.fields || {});
    var c = p.controls || {};
    Object.keys(c).forEach(function (id) {
      var el = $(id);
      if (!el) return;
      if (el.type === 'checkbox') el.checked = !!c[id]; else el.value = c[id];
    });
    $('smearW').disabled = !$('smear').checked;
    buildFields(); updateGridInfo(); updateDisplayLabels();
    state.needGen = true; request(true);
  }

  function exportCSV() {
    var s = state.scat;
    if (!s) return;
    var S = sOfQ(s), lines = [];
    lines.push('# SAXS-FFT 2D  structure=' + state.struct + '  params=' + JSON.stringify(allParams()));
    lines.push('# I in units of [rho^2 nm^2] (intensity per unit area); q in nm^-1');
    lines.push('q,I,I_particles_S1,S,n_pixels');
    for (var i = 0; i < s.profile.q.length; i++)
      lines.push([s.profile.q[i], s.profile.I[i], s.dilute ? s.dilute[i] : '', S && isFinite(S[i]) ? S[i] : '', s.profile.n[i]].join(','));
    download(new Blob([lines.join('\n')], { type: 'text/csv' }), baseName() + '.csv');
  }

  /* ---------- wiring ---------- */

  function updateDisplayLabels() {
    $('decadesOut').textContent = $('decades').value + ' dec';
    var t = parseFloat($('top').value);
    var lin = $('scale2D').value === 'lin';
    $('topOut').textContent = t === 0 ? 'max' : lin ? 'max × ' + fmt(Math.pow(10, t)) : t + ' dec';
    $('decades').disabled = lin;
    $('title2D').textContent = lin ? 'Scattering pattern I(qx, qy)' : 'Scattering pattern log₁₀ I(qx, qy)';
  }

  selStruct.addEventListener('change', function () {
    state.struct = selStruct.value;
    var minN = STRUCTURES[state.struct].minN;
    if (minN && parseInt($('N').value, 10) < minN) { $('N').value = String(minN); updateGridInfo(); }
    buildFields(); structureChanged();
  });
  ['N', 'a', 'seed', 'rhoM', 'rhoP', 'invert'].forEach(function (id) {
    $(id).addEventListener('change', function () { updateGridInfo(); structureChanged(); });
  });
  $('reseed').addEventListener('click', function () { $('seed').value = 1 + Math.floor(Math.random() * 99999); structureChanged(); });
  ['meanSub', 'window', 'pad', 'range', 'pixelMTF', 'smear', 'smearW'].forEach(function (id) { $(id).addEventListener('change', scatteringChanged); });
  $('smear').addEventListener('change', function () { $('smearW').disabled = !$('smear').checked; });
  ['cmap', 'scale2D', 'decades', 'top'].forEach(function (id) {
    $(id).addEventListener('input', function () { updateDisplayLabels(); draw2D(false); });
  });
  ['ovDilute', 'ovPorod', 'ovGuinier'].forEach(function (id) { $(id).addEventListener('change', function () { draw1D(true); }); });
  ['xlog', 'ylog'].forEach(function (id) { $(id).addEventListener('change', function () { draw1D(false); }); });
  $('run').addEventListener('click', function () { state.needGen = true; request(true); });
  $('stop').addEventListener('click', function () {
    if (engine.restart()) { state.needGen = true; state.reqId++; setBusy(false); setStatus('Stopped'); }
  });
  document.querySelectorAll('.tab').forEach(function (b) {
    b.addEventListener('click', function () {
      document.querySelectorAll('.tab').forEach(function (x) { x.classList.toggle('active', x === b); });
      state.tab = b.dataset.tab; draw1D(false);
    });
  });
  $('pngStruct').addEventListener('click', function () { exportPNG($('cvStruct'), 'structure'); });
  $('png2D').addEventListener('click', function () { exportPNG($('cv2D'), 'pattern2D'); });
  $('png1D').addEventListener('click', function () { exportPNG($('cv1D'), state.tab === 'I' ? 'Iq' : 'Sq'); });
  $('csv').addEventListener('click', exportCSV);
  $('saveJSON').addEventListener('click', function () {
    download(new Blob([JSON.stringify(allParams(), null, 2)], { type: 'application/json' }), baseName() + '.json');
  });
  $('loadJSON').addEventListener('click', function () { $('fileJSON').click(); });
  $('fileJSON').addEventListener('change', function () {
    var f = this.files[0]; if (!f) return;
    var r = new FileReader();
    r.onload = function () { try { applyParams(JSON.parse(r.result)); } catch (e) { setStatus('Load failed: ' + e.message, 'error'); } };
    r.readAsText(f); this.value = '';
  });
  if (window.matchMedia) window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', function () {
    structView.draw(); view2D.draw(); draw1D(true);
  });

  buildFields();
  updateGridInfo();
  updateDisplayLabels();
  request(true);

  window.SAXS_app = { state: state, request: request, applyParams: applyParams, STRUCTURES: STRUCTURES };
})();
