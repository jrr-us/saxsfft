/* SAXS-FFT core: microstructure generators.
 *
 * All lengths are in nm. The simulation box is L x L with L = N * a and
 * periodic boundary conditions. A generator returns
 *   { particles: [...] | null, field: Float64Array | null,
 *     form: description used for analytic curves | null,
 *     info: [text lines], rhoP: particle SLD used for Porod/invert }
 * Particles: {kind:'disk'|'ell'|'cs', x, y, a, b, th, s, rb, rho | rhoC,rhoS,rc}
 *   a, b = semi-axes (a = b = radius for disks / outer radius for core-shell),
 *   s = size scale factor (polydispersity), rb = bounding radius.
 */
(function (root) {
  root.SAXS_PARTS = root.SAXS_PARTS || [];
  root.SAXS_PARTS.push(function partStructures(S) {
    var PI = Math.PI;
    function now() { return (typeof performance !== 'undefined' ? performance : Date).now(); }
    function wrap(v, L) { v %= L; return v < 0 ? v + L : v; }

    /* ---------- particle geometry ---------- */

    function setEllipseMatrix(p) {
      var c = Math.cos(p.th), s = Math.sin(p.th), a2 = p.a * p.a, b2 = p.b * p.b;
      p.m11 = a2 * c * c + b2 * s * s;
      p.m22 = a2 * s * s + b2 * c * c;
      p.m12 = (a2 - b2) * c * s;
    }

    /* Largest common scale factor g such that p and q (scaled by g) do not
       overlap. Disks: d / (ra + rb). Ellipses: sqrt of the Perram-Wertheim
       contact function max_l l(1-l) r^T [(1-l)A + l B]^-1 r. */
    function pairRatio(p, q, dx, dy) {
      if (p.kind !== 'ell' && q.kind !== 'ell')
        return Math.sqrt(dx * dx + dy * dy) / (p.a + q.a);
      var lo = 0, hi = 1, gr = 0.6180339887498949;
      function f(l) {
        var m = 1 - l;
        var c11 = m * p.m11 + l * q.m11, c22 = m * p.m22 + l * q.m22, c12 = m * p.m12 + l * q.m12;
        var det = c11 * c22 - c12 * c12;
        return l * m * (c22 * dx * dx - 2 * c12 * dx * dy + c11 * dy * dy) / det;
      }
      var x1 = hi - gr * (hi - lo), x2 = lo + gr * (hi - lo), f1 = f(x1), f2 = f(x2);
      for (var it = 0; it < 28; it++) {
        if (f1 < f2) { lo = x1; x1 = x2; f1 = f2; x2 = lo + gr * (hi - lo); f2 = f(x2); }
        else { hi = x2; x2 = x1; f2 = f1; x1 = hi - gr * (hi - lo); f1 = f(x1); }
      }
      return Math.sqrt(Math.max(f1, f2));
    }

    /* ---------- periodic cell list ---------- */

    function CellGrid(L, minCell) {
      var nc = Math.max(1, Math.floor(L / minCell));
      if (nc < 3) nc = 1;                    // too few cells: one cell holds everything
      this.nc = nc; this.cs = L / nc; this.L = L;
      this.cells = [];
      for (var i = 0; i < nc * nc; i++) this.cells.push([]);
    }
    CellGrid.prototype.index = function (x, y) {
      var nc = this.nc;
      var ix = Math.floor(x / this.cs) % nc, iy = Math.floor(y / this.cs) % nc;
      if (ix < 0) ix += nc; if (iy < 0) iy += nc;
      return iy * nc + ix;
    };
    CellGrid.prototype.add = function (p) { p.cell = this.index(p.x, p.y); this.cells[p.cell].push(p); };
    CellGrid.prototype.remove = function (p) {
      var c = this.cells[p.cell], k = c.indexOf(p);
      c[k] = c[c.length - 1]; c.pop();
    };
    /* Call fn(q, dx, dy) for every particle q in the 3x3 neighbourhood of (x, y),
       with (dx, dy) the minimum-image vector from (x, y) to q. Stops if fn returns false. */
    CellGrid.prototype.forNeighbours = function (x, y, fn) {
      var nc = this.nc, L = this.L, cells = this.cells, c, k, q, dx, dy;
      if (nc === 1) {
        c = cells[0];
        for (k = 0; k < c.length; k++) {
          q = c[k]; dx = q.x - x; dy = q.y - y;
          dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L);
          if (fn(q, dx, dy) === false) return false;
        }
        return true;
      }
      var ix = Math.floor(x / this.cs), iy = Math.floor(y / this.cs);
      for (var oy = -1; oy <= 1; oy++) {
        var cy = ((iy + oy) % nc + nc) % nc;
        for (var ox = -1; ox <= 1; ox++) {
          c = cells[cy * nc + ((ix + ox) % nc + nc) % nc];
          for (k = 0; k < c.length; k++) {
            q = c[k]; dx = q.x - x; dy = q.y - y;
            dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L);
            if (fn(q, dx, dy) === false) return false;
          }
        }
      }
      return true;
    };

    function maxBound(parts) {
      var m = 0;
      for (var i = 0; i < parts.length; i++) if (parts[i].rb > m) m = parts[i].rb;
      return m;
    }

    /* Would particle p at (x, y), scaled by g, overlap anything in the grid? */
    function overlaps(grid, p, x, y, g) {
      return !grid.forNeighbours(x, y, function (q, dx, dy) {
        if (q === p) return true;
        var d2 = dx * dx + dy * dy, sb = g * (p.rb + q.rb);
        if (d2 >= sb * sb) return true;
        return pairRatio(p, q, dx, dy) >= g;
      });
    }

    /* ---------- placement strategies ---------- */

    /* Random sequential addition (largest first). Returns the placed particles. */
    function placeRSA(parts, L, rng, g, maxTries) {
      g = g || 1; maxTries = maxTries || 3000;
      parts.sort(function (p, q) { return q.rb - p.rb; });
      var grid = new CellGrid(L, 2 * g * maxBound(parts) + 1e-9), placed = [];
      for (var i = 0; i < parts.length; i++) {
        var p = parts[i];
        for (var t = 0; t < maxTries; t++) {
          var x = rng.uniform() * L, y = rng.uniform() * L;
          if (!overlaps(grid, p, x, y, g)) { p.x = x; p.y = y; grid.add(p); placed.push(p); break; }
        }
      }
      return placed;
    }

    function placeOverlapping(parts, L, rng) {
      for (var i = 0; i < parts.length; i++) { parts[i].x = rng.uniform() * L; parts[i].y = rng.uniform() * L; }
      return parts;
    }

    /* Remove overlaps by pushing overlapping pairs apart along their centre
       line (soft-particle relaxation). The contact ratio scales linearly with
       distance, so a pair with ratio r needs its distance multiplied by 1/r.
       Returns the number of iterations used, or -1 if overlaps remain. */
    function relaxOverlaps(parts, L, rng, maxIter, progress) {
      var n = parts.length, bmax = maxBound(parts), dxs = new Float64Array(n), dys = new Float64Array(n), last = now();
      parts.forEach(function (p, i) { p.id = i; });
      for (var it = 0; it < maxIter; it++) {
        var grid = new CellGrid(L, 2 * bmax + 1e-9), nover = 0;
        parts.forEach(function (p) { grid.add(p); });
        dxs.fill(0); dys.fill(0);
        for (var i = 0; i < n; i++) {
          var p = parts[i];
          grid.forNeighbours(p.x, p.y, function (q, dx, dy) {
            if (q.id <= p.id) return true;
            var d2 = dx * dx + dy * dy, sb = p.rb + q.rb;
            if (d2 >= sb * sb) return true;
            var d = Math.sqrt(d2);
            if (d < 1e-9) { var an = 2 * PI * rng.uniform(); dx = Math.cos(an) * 1e-3; dy = Math.sin(an) * 1e-3; d = 1e-3; }
            var r = pairRatio(p, q, dx, dy);
            if (r >= 1) return true;
            nover++;
            var push = Math.min(0.5 * (1.002 / Math.max(r, 1e-3) - 1) * d, 0.5 * sb) / d;
            dxs[p.id] -= push * dx; dys[p.id] -= push * dy;
            dxs[q.id] += push * dx; dys[q.id] += push * dy;
            return true;
          });
        }
        if (!nover) return it;
        for (i = 0; i < n; i++) { parts[i].x = wrap(parts[i].x + dxs[i], L); parts[i].y = wrap(parts[i].y + dys[i], L); }
        if (progress && now() - last > 300) { last = now(); progress('Removing overlaps: ' + nover + ' left (iteration ' + it + ')'); }
      }
      return -1;
    }

    /* Dense hard-particle fluid: random positions at full size, overlaps
       removed by relaxation, then eqSweeps Metropolis hard-particle MC sweeps
       to reach equilibrium. If relaxation fails, particles are shrunk by the
       smallest contact ratio g. Returns {g, sweeps, relax}. */
    function placeDense(parts, L, rng, phi, opt) {
      var n = parts.length, rotate = !!opt.rotate, eqSweeps = opt.sweeps || 300;
      var budget = opt.timeBudget || 20000, t0 = now(), lastProg = t0;
      placeOverlapping(parts, L, rng);
      var relax = relaxOverlaps(parts, L, rng, 3000, opt.progress);
      var grid = new CellGrid(L, 2 * maxBound(parts) + 1e-9);
      parts.forEach(function (p) { grid.add(p); });
      var g = relax < 0 ? Math.min(1, 0.9999 * minRatio(parts, grid)) : 1;

      var mean = 0; parts.forEach(function (p) { mean += p.rb; }); mean /= n;
      var delta = 0.1 * mean, dth = 0.1, eq = 0, sweeps = 0;
      while (eqSweeps > 0) {
        var acc = 0;
        for (var m = 0; m < n; m++) {
          var p = parts[(rng.uniform() * n) | 0];
          var ox = p.x, oy = p.y, oth = p.th;
          var nx = wrap(ox + delta * (2 * rng.uniform() - 1), L);
          var ny = wrap(oy + delta * (2 * rng.uniform() - 1), L);
          if (rotate) { p.th = oth + dth * (2 * rng.uniform() - 1); setEllipseMatrix(p); }
          if (!overlaps(grid, p, nx, ny, g)) {
            grid.remove(p); p.x = nx; p.y = ny; grid.add(p); acc++;
          } else if (rotate) { p.th = oth; setEllipseMatrix(p); }
        }
        sweeps++;
        var rate = acc / n;
        delta *= rate > 0.35 ? 1.08 : 0.92;
        delta = Math.min(Math.max(delta, 1e-4 * mean), L / 4);
        dth = Math.min(Math.max(dth * (rate > 0.35 ? 1.08 : 0.92), 1e-3), PI / 2);
        if (++eq >= eqSweeps) break;
        var t = now();
        if (t - t0 > budget) break;
        if (opt.progress && t - lastProg > 300) {
          lastProg = t;
          opt.progress('Equilibrating (hard-particle MC): sweep ' + eq + ' / ' + eqSweeps);
        }
      }
      return { g: g, sweeps: sweeps, relax: relax };
    }

    function minRatio(parts, grid) {
      var m = Infinity;
      for (var i = 0; i < parts.length; i++) {
        var p = parts[i];
        grid.forNeighbours(p.x, p.y, function (q, dx, dy) {
          if (q === p) return true;
          var d2 = dx * dx + dy * dy, sb = p.rb + q.rb;
          if (d2 >= sb * sb) return true;          // ratio >= 1 anyway
          var r = pairRatio(p, q, dx, dy);
          if (r < m) m = r;
          return true;
        });
      }
      return m;
    }

    function scaleParticles(parts, g) {
      if (g === 1) return;
      parts.forEach(function (p) {
        p.a *= g; p.b *= g; p.rb *= g; p.s *= g; if (p.rc) p.rc *= g;
        if (p.kind === 'ell') setEllipseMatrix(p);
      });
    }

    /* Draw particles until their total area reaches phi * L^2 (or count n). */
    function drawPopulation(P, L, rng, make, area) {
      var parts = [], A = 0, target = P.phi * L * L, maxN = 200000;
      if (P.amount === 'count') {
        for (var i = 0; i < P.count; i++) { var p = make(S.sizeScale(rng, P.sigma, P.dist)); parts.push(p); A += area(p); }
      } else {
        while (A < target && parts.length < maxN) { var q = make(S.sizeScale(rng, P.sigma, P.dist)); parts.push(q); A += area(q); }
      }
      return parts;
    }

    function arrange(parts, P, L, rng, progress, rotate) {
      var info = [], g = 1, totA = 0, placed;
      parts.forEach(function (p) { totA += PI * p.a * p.b; });
      var phiReq = totA / (L * L);
      if (P.packing === 'overlap') {
        placed = placeOverlapping(parts, L, rng);
        info.push('Positions: Poisson (overlaps allowed)');
      } else if (P.packing === 'dense') {
        if (phiReq > 0.85) throw new Error('Area fraction too high for a hard-particle fluid (max ≈ 0.72)');
        var r = placeDense(parts, L, rng, phiReq, { sweeps: P.sweeps, rotate: rotate, progress: progress });
        g = r.g; scaleParticles(parts, g); placed = parts;
        info.push('Overlap relaxation: ' + (r.relax < 0 ? 'incomplete' : r.relax + ' iterations') + '; hard-particle MC: ' + r.sweeps + ' sweeps');
        if (g < 1) info.push('⚠ could not remove all overlaps: sizes shrunk by ' + g.toFixed(3));
      } else {
        placed = placeRSA(parts, L, rng, 1, 3000);
        info.push('Positions: random sequential addition (no overlap)');
        if (placed.length < parts.length)
          info.push('⚠ RSA placed only ' + placed.length + ' of ' + parts.length + ' particles (jamming ≈ 0.547); use "dense"');
      }
      var A = 0; placed.forEach(function (p) { A += PI * p.a * p.b; });
      info.unshift('Particles: ' + placed.length + ',  area fraction φ = ' + (A / (L * L)).toFixed(4));
      return { parts: placed, info: info };
    }

    /* ---------- generators ---------- */

    function disk(x, y, r, rho, s) {
      return { kind: 'disk', x: x, y: y, a: r, b: r, th: 0, s: s, rb: r, rho: rho };
    }

    function genSingle(P, L) {
      return {
        particles: [disk(L / 2, L / 2, P.R, P.rhoP, 1)],
        form: { kind: 'disk', R: P.R },
        info: ['Single disk, R = ' + P.R + ' nm,  φ = ' + (PI * P.R * P.R / (L * L)).toFixed(5)]
      };
    }

    function genDisks(P, L, rng, progress) {
      var parts = drawPopulation(P, L, rng, function (s) { return disk(0, 0, P.R * s, P.rhoP, s); },
                                 function (p) { return PI * p.a * p.a; });
      var r = arrange(parts, P, L, rng, progress, false);
      return { particles: r.parts, form: { kind: 'disk', R: P.R }, info: r.info };
    }

    function genEllipses(P, L, rng, progress) {
      var parts = drawPopulation(P, L, rng, function (s) {
        var th = P.orient === 'aligned' ? (P.theta + P.sigmaTheta * rng.normal()) * PI / 180 : rng.uniform() * PI;
        var p = { kind: 'ell', x: 0, y: 0, a: P.ea * s, b: P.eb * s, th: th, s: s, rb: Math.max(P.ea, P.eb) * s, rho: P.rhoP };
        setEllipseMatrix(p); return p;
      }, function (p) { return PI * p.a * p.b; });
      var r = arrange(parts, P, L, rng, progress, P.orient !== 'aligned');
      return { particles: r.parts, form: { kind: 'ell', a: P.ea, b: P.eb }, info: r.info };
    }

    function genCoreShell(P, L, rng, progress) {
      var ro = P.Rc + P.t;
      var parts = drawPopulation(P, L, rng, function (s) {
        return { kind: 'cs', x: 0, y: 0, a: ro * s, b: ro * s, th: 0, s: s, rb: ro * s, rc: P.Rc * s, rhoC: P.rhoC, rhoS: P.rhoS };
      }, function (p) { return PI * p.a * p.a; });
      var r = arrange(parts, P, L, rng, progress, false);
      return { particles: r.parts, form: { kind: 'cs', rc: P.Rc, ro: ro, rhoC: P.rhoC, rhoS: P.rhoS }, info: r.info };
    }

    /* Off-lattice diffusion-limited aggregation of M disks of radius r0. */
    function dlaCluster(M, r0, rng) {
      // Fine grid (cell 2 r0, linked lists) for contact search and a coarse
      // occupancy grid (cell 8 r0) that lets walkers take long steps through
      // empty regions. Both are square, centred on the seed, and regrown.
      var xs = new Float64Array(M), ys = new Float64Array(M), next = new Int32Array(M), n = 0, Rmax = 0;
      var h = 2 * r0, C = 8 * r0, G, head, half, CG, occ, chalf;
      function build(R) {
        half = Math.ceil(R / h) + 3; G = 2 * half; head = new Int32Array(G * G).fill(-1);
        chalf = Math.ceil(R / C) + 6; CG = 2 * chalf; occ = new Uint8Array(CG * CG);
        for (var i = 0; i < n; i++) insert(i);
      }
      function insert(i) {
        var cx = Math.floor(xs[i] / h) + half, cy = Math.floor(ys[i] / h) + half, c = cy * G + cx;
        next[i] = head[c]; head[c] = i;
        occ[(Math.floor(ys[i] / C) + chalf) * CG + Math.floor(xs[i] / C) + chalf] = 1;
      }
      function addPt(x, y) {
        xs[n] = x; ys[n] = y; n++;
        var r = Math.sqrt(x * x + y * y);
        if (r > Rmax) Rmax = r;
        if (Rmax + 8 * C > half * h || Rmax + 8 * C > chalf * C) build(2 * Rmax + 10 * C);
        else insert(n - 1);
      }
      build(20 * r0);
      addPt(0, 0);
      while (n < M) {
        var Rl = Rmax + 5 * r0, ang = 2 * PI * rng.uniform();
        var x = Rl * Math.cos(ang), y = Rl * Math.sin(ang), Rkill = 2 * Rmax + 30 * r0;
        while (true) {
          var rr = Math.sqrt(x * x + y * y), step;
          if (rr > Rkill) { ang = 2 * PI * rng.uniform(); x = Rl * Math.cos(ang); y = Rl * Math.sin(ang); continue; }
          if (rr > Rmax + 4 * r0) step = rr - Rmax - 2 * r0;
          else {
            // nearest occupied coarse ring k: particles there are >= (k-1) C away
            var ccx = Math.floor(x / C) + chalf, ccy = Math.floor(y / C) + chalf, kr = 0, found = false;
            for (kr = 0; kr <= 4 && !found; kr++) {
              for (var oy = -kr; oy <= kr && !found; oy++) {
                var row = (ccy + oy) * CG, stride = (oy === -kr || oy === kr) ? 1 : 2 * kr;
                for (var ox = -kr; ox <= kr; ox += stride) if (occ[row + ccx + ox]) { found = true; break; }
              }
            }
            kr--;                                        // ring that contained a particle (or 4)
            var best = -1, dmin;
            if (found && kr <= 1) {
              var ix = Math.floor(x / h) + half, iy = Math.floor(y / h) + half;
              dmin = 4 * r0;
              for (oy = -2; oy <= 2; oy++) for (ox = -2; ox <= 2; ox++) {
                for (var k = head[(iy + oy) * G + ix + ox]; k >= 0; k = next[k]) {
                  var dx = x - xs[k], dy = y - ys[k], d = Math.sqrt(dx * dx + dy * dy);
                  if (d < dmin) { dmin = d; best = k; }
                }
              }
            } else dmin = (found ? kr - 1 : 4) * C;
            step = dmin - 2 * r0;
            if (best >= 0 && step < 0.1 * r0) {          // stick: snap to contact
              var f = 2 * r0 / dmin;
              addPt(xs[best] + (x - xs[best]) * f, ys[best] + (y - ys[best]) * f);
              break;
            }
          }
          var a2 = 2 * PI * rng.uniform();
          x += step * Math.cos(a2); y += step * Math.sin(a2);
        }
      }
      // mass-radius dimension from the seed: N(r) ~ r^Df between 5 r0 and Rmax / 2
      var rs = [];
      for (var i = 0; i < n; i++) rs.push(Math.sqrt(xs[i] * xs[i] + ys[i] * ys[i]));
      rs.sort(function (p, q) { return p - q; });
      var lx = [], ly = [];
      for (var r = 5 * r0; r < Rmax / 2; r *= 1.25) {
        var lo = 0, hi = n;
        while (lo < hi) { var mid = (lo + hi) >> 1; if (rs[mid] < r) lo = mid + 1; else hi = mid; }
        lx.push(Math.log(r)); ly.push(Math.log(lo));
      }
      var df = NaN;
      if (lx.length > 2) {
        var mx = 0, my = 0, sxy = 0, sxx = 0;
        for (i = 0; i < lx.length; i++) { mx += lx[i]; my += ly[i]; }
        mx /= lx.length; my /= lx.length;
        for (i = 0; i < lx.length; i++) { sxy += (lx[i] - mx) * (ly[i] - my); sxx += (lx[i] - mx) * (lx[i] - mx); }
        df = sxy / sxx;
      }
      var cx = 0, cy = 0, rg2 = 0, rb = 0;
      for (i = 0; i < n; i++) { cx += xs[i]; cy += ys[i]; }
      cx /= n; cy /= n;
      for (i = 0; i < n; i++) {
        xs[i] -= cx; ys[i] -= cy;
        var d2 = xs[i] * xs[i] + ys[i] * ys[i]; rg2 += d2; rb = Math.max(rb, Math.sqrt(d2));
      }
      return { xs: xs, ys: ys, rb: rb + r0, rg: Math.sqrt(rg2 / n + r0 * r0 / 2), df: df };
    }

    function genDLA(P, L, rng, progress) {
      var clusters = [], parts = [], info = [];
      for (var c = 0; c < P.nClusters; c++) {
        if (progress) progress('Growing cluster ' + (c + 1) + ' / ' + P.nClusters);
        clusters.push(dlaCluster(P.M, P.r0, rng));
      }
      // place clusters by RSA on their bounding circles, falling back to overlap
      var holders = clusters.map(function (cl, i) { return { kind: 'disk', a: cl.rb, b: cl.rb, rb: cl.rb, idx: i, th: 0 }; });
      var placed = placeRSA(holders.slice(), L, rng, 1, 2000), overl = placed.length < holders.length;
      if (overl) placeOverlapping(holders.filter(function (h) { return placed.indexOf(h) < 0; }), L, rng);
      var rgm = 0, dfm = 0;
      holders.forEach(function (h) {
        var cl = clusters[h.idx], rot = 2 * PI * rng.uniform(), cs = Math.cos(rot), sn = Math.sin(rot);
        rgm += cl.rg / clusters.length; dfm += cl.df / clusters.length;
        for (var i = 0; i < cl.xs.length; i++)
          parts.push(disk(wrap(h.x + cs * cl.xs[i] - sn * cl.ys[i], L), wrap(h.y + sn * cl.xs[i] + cs * cl.ys[i], L), P.r0, P.rhoP, 1));
      });
      info.push('Particles: ' + parts.length + ' in ' + clusters.length + ' DLA clusters');
      info.push('Mass–radius dimension N(r) ∝ r^Df: Df ≈ ' + dfm.toFixed(2) + ' (2D DLA: 1.71)');
      info.push('Mean cluster Rg ≈ ' + rgm.toFixed(1) + ' nm,  2π/Rg ≈ ' + (2 * PI / rgm).toFixed(4) + ' nm⁻¹');
      if (overl) info.push('⚠ clusters did not fit without overlapping; use a larger box');
      if (2 * clusters[0].rb > L) info.push('⚠ cluster larger than the box: it overlaps its periodic image');
      // Guinier of whole clusters: I(0) = sum_clusters (M pi r0^2 drho)^2 / A
      var m = P.M * PI * P.r0 * P.r0 * (P.rhoP - P.rhoM);
      return { particles: parts, form: { kind: 'disk', R: P.r0 }, info: info,
               guinier: { I0: clusters.length * m * m / (L * L), Rg: rgm, label: 'clusters' } };
    }

    /* Gaussian random field thresholded at area fraction phi ("spinodal"). */
    function genSpinodal(P, L, rng, progress, N, a) {
      var n2 = N * N, re = new Float64Array(n2), im = new Float64Array(n2), i, x, y;
      for (i = 0; i < n2; i++) re[i] = rng.normal();
      S.fft2(re, im, N, false);
      var dq = 2 * PI / L, q0 = 2 * PI / P.d0, w = P.width * q0;
      for (y = 0; y < N; y++) {
        var ky = y <= N / 2 ? y : y - N;
        for (x = 0; x < N; x++) {
          var kx = x <= N / 2 ? x : x - N, q = dq * Math.sqrt(kx * kx + ky * ky), f;
          if (P.filter === 'ring') f = Math.exp(-(q - q0) * (q - q0) / (2 * w * w));
          else f = Math.exp(-q * q * P.d0 * P.d0 / (16 * PI * PI));   // Gaussian, correlation ~ d0
          if (x === 0 && y === 0) f = 0;
          re[y * N + x] *= f; im[y * N + x] *= f;
        }
      }
      S.fft2(re, im, N, true);
      var sorted = Float64Array.from(re).sort();
      var t = sorted[Math.min(n2 - 1, Math.max(0, Math.round((1 - P.phi) * n2)))];
      var cov = new Float64Array(n2), sum = 0;
      for (y = 0; y < N; y++) for (x = 0; x < N; x++) {        // antialias with the local gradient
        var v = re[y * N + x] - t;
        var gx = (re[y * N + (x + 1) % N] - re[y * N + (x + N - 1) % N]) / 2;
        var gy = (re[((y + 1) % N) * N + x] - re[((y + N - 1) % N) * N + x]) / 2;
        var gm = Math.sqrt(gx * gx + gy * gy) + 1e-30;
        var c = 0.5 + v / gm;
        c = c < 0 ? 0 : c > 1 ? 1 : c;
        cov[y * N + x] = c; sum += c;
      }
      return {
        field: cov, form: null,
        info: ['Gaussian random field, ' + (P.filter === 'ring' ? 'ring filter at q₀ = 2π/' + P.d0 + ' nm' : 'Gaussian filter, length ' + P.d0 + ' nm'),
               'Thresholded: φ = ' + (sum / n2).toFixed(4)]
      };
    }

    var GENERATORS = {
      single: genSingle, disks: genDisks, ellipses: genEllipses, coreshell: genCoreShell,
      spinodal: genSpinodal, dla: genDLA
    };

    /* P: full parameter object (P.gen selects the generator). */
    function generate(P, progress) {
      var N = P.N, a = P.a, L = N * a, rng = S.makeRng(P.seed);
      var fn = GENERATORS[P.gen];
      if (!fn) throw new Error('Unknown structure ' + P.gen);
      var out = fn(P, L, rng, progress, N, a);
      out.particles = out.particles || null;
      out.field = out.field || null;
      out.L = L;
      return out;
    }

    S.generate = generate;
    S.pairRatio = pairRatio;
    S.setEllipseMatrix = setEllipseMatrix;
  });
})(typeof self !== 'undefined' ? self : globalThis);
