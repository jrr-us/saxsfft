/* Minimal canvas plotting: LinePlot (1D, log/lin axes, zoom, hover) and
 * ImageView (2D image with physical axes, square zoom, hover). Colours
 * come from CSS custom properties so plots follow the light/dark theme. */
(function (root) {
  function css(name, fallback) {
    var v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
    return v || fallback;
  }

  function niceStep(span, n) {
    var raw = span / Math.max(n, 1), p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
    return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p;
  }

  function linTicks(lo, hi, n) {
    var st = niceStep(hi - lo, n), t = [], v = Math.ceil(lo / st - 1e-9) * st;
    for (; v <= hi + st * 1e-9; v += st) t.push({ v: Math.abs(v) < st * 1e-9 ? 0 : v, major: true });
    return t;
  }

  function logTicks(lo, hi) {           // lo, hi in data units (> 0)
    var t = [], d0 = Math.floor(Math.log10(lo)), d1 = Math.ceil(Math.log10(hi)), dec = d1 - d0;
    var every = Math.max(1, Math.ceil(dec / 8));
    for (var d = d0; d <= d1; d++) for (var m = 1; m < 10; m++) {
      var v = m * Math.pow(10, d);
      if (v < lo * (1 - 1e-9) || v > hi * (1 + 1e-9)) continue;
      if (m === 1 && (d % every === 0)) t.push({ v: v, major: true });
      else if (dec <= 6) t.push({ v: v, major: false, label: dec <= 1.5 && (m === 2 || m === 5) });
    }
    return t;
  }

  function fmt(v) {
    if (v === 0) return '0';
    var a = Math.abs(v);
    if (a >= 1e5 || a < 1e-3) {
      var e = Math.floor(Math.log10(a)), m = v / Math.pow(10, e);
      return (Math.abs(m - 1) < 1e-9 ? '' : +m.toFixed(1) + '×') + '10' + sup(e);
    }
    return a >= 100 ? String(Math.round(v)) : String(+v.toPrecision(3));
  }
  function sup(n) {
    var map = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
    return String(n).split('').map(function (c) { return map[c]; }).join('');
  }

  function setupCanvas(canvas) {
    var dpr = window.devicePixelRatio || 1, w = canvas.clientWidth, h = canvas.clientHeight;
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
      canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
    }
    var ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx: ctx, w: w, h: h };
  }

  /* Attach drag-to-zoom (rubber band) + dblclick reset + hover. */
  function attachMouse(canvas, self) {
    var drag = null;
    function pos(e) { var r = canvas.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; }
    canvas.addEventListener('mousedown', function (e) { if (e.button === 0) { drag = { a: pos(e), b: pos(e) }; } });
    window.addEventListener('mousemove', function (e) {
      if (drag) { drag.b = pos(e); self.band = drag; self.draw(); }
    });
    canvas.addEventListener('mousemove', function (e) { if (!drag) self.hover(pos(e)); });
    canvas.addEventListener('mouseleave', function () { self.hover(null); });
    window.addEventListener('mouseup', function () {
      if (!drag) return;
      var d = drag; drag = null; self.band = null;
      if (Math.abs(d.a.x - d.b.x) > 6 && Math.abs(d.a.y - d.b.y) > 6) self.zoomTo(d.a, d.b);
      else self.draw();
    });
    canvas.addEventListener('dblclick', function () { self.resetView(); });
    if (window.ResizeObserver) new ResizeObserver(function () { self.draw(); }).observe(canvas);
  }

  /* ---------------- LinePlot ---------------- */

  function LinePlot(canvas, opts) {
    this.c = canvas; this.o = opts || {}; this.series = []; this.view = null; this.mouse = null;
    this.xlog = !!this.o.xlog; this.ylog = !!this.o.ylog;
    attachMouse(canvas, this);
  }
  LinePlot.prototype.setData = function (series, keepView) {
    this.series = series || [];
    if (!keepView) this.view = null;
    this.draw();
  };
  LinePlot.prototype.setScales = function (xlog, ylog) { this.xlog = xlog; this.ylog = ylog; this.view = null; this.draw(); };
  LinePlot.prototype.resetView = function () { this.view = null; this.draw(); };
  LinePlot.prototype.autoRange = function () {
    var xl = Infinity, xh = -Infinity, yl = Infinity, yh = -Infinity, self = this;
    this.series.forEach(function (s) {
      if (s.noAuto) return;
      for (var i = 0; i < s.x.length; i++) {
        var x = s.x[i], y = s.y[i];
        if (!isFinite(x) || !isFinite(y) || (self.xlog && x <= 0) || (self.ylog && y <= 0)) continue;
        if (x < xl) xl = x; if (x > xh) xh = x; if (y < yl) yl = y; if (y > yh) yh = y;
      }
    });
    if (!isFinite(xl)) return { x0: 1, x1: 10, y0: 1, y1: 10 };
    if (this.ylog && this.o.maxDecades && yh / yl > Math.pow(10, this.o.maxDecades)) yl = yh / Math.pow(10, this.o.maxDecades);
    var v = { x0: xl, x1: xh, y0: yl, y1: yh };
    if (this.ylog) { v.y0 /= 1.5; v.y1 *= 1.5; } else { var p = 0.05 * (yh - yl || 1); v.y0 -= p; v.y1 += p; if (this.o.yZero && yl >= 0) v.y0 = 0; }
    if (!this.xlog) { if (this.o.xZero && xl >= 0) v.x0 = 0; }
    if (v.x0 === v.x1) { v.x0 *= 0.9; v.x1 *= 1.1; }
    return v;
  };
  LinePlot.prototype.frame = function () {
    var m = { l: 64, r: 14, t: 12, b: 42 };
    return { l: m.l, t: m.t, w: Math.max(10, this.c.clientWidth - m.l - m.r), h: Math.max(10, this.c.clientHeight - m.t - m.b) };
  };
  LinePlot.prototype.tx = function (v, F, V) {
    return this.xlog ? F.l + F.w * Math.log(v / V.x0) / Math.log(V.x1 / V.x0) : F.l + F.w * (v - V.x0) / (V.x1 - V.x0);
  };
  LinePlot.prototype.ty = function (v, F, V) {
    return this.ylog ? F.t + F.h * (1 - Math.log(v / V.y0) / Math.log(V.y1 / V.y0)) : F.t + F.h * (1 - (v - V.y0) / (V.y1 - V.y0));
  };
  LinePlot.prototype.inv = function (px, py, F, V) {
    var fx = (px - F.l) / F.w, fy = 1 - (py - F.t) / F.h;
    return {
      x: this.xlog ? V.x0 * Math.pow(V.x1 / V.x0, fx) : V.x0 + fx * (V.x1 - V.x0),
      y: this.ylog ? V.y0 * Math.pow(V.y1 / V.y0, fy) : V.y0 + fy * (V.y1 - V.y0)
    };
  };
  LinePlot.prototype.zoomTo = function (a, b) {
    var F = this.frame(), V = this.view || this.autoRange(), p = this.inv(a.x, a.y, F, V), q = this.inv(b.x, b.y, F, V);
    this.view = { x0: Math.min(p.x, q.x), x1: Math.max(p.x, q.x), y0: Math.min(p.y, q.y), y1: Math.max(p.y, q.y) };
    this.draw();
  };
  LinePlot.prototype.hover = function (p) { this.mouse = p; this.draw(); };
  LinePlot.prototype.draw = function () {
    var g = setupCanvas(this.c), ctx = g.ctx, F = this.frame(), V = this.view || this.autoRange(), self = this;
    var fg = css('--fg', '#222'), muted = css('--muted', '#777'), grid = css('--grid', '#e5e5e5'), bg = css('--panel', '#fff');
    ctx.fillStyle = bg; ctx.fillRect(0, 0, g.w, g.h);
    ctx.font = '11px system-ui, sans-serif';
    // grid + ticks
    var xt = this.xlog ? logTicks(V.x0, V.x1) : linTicks(V.x0, V.x1, F.w / 80);
    var yt = this.ylog ? logTicks(V.y0, V.y1) : linTicks(V.y0, V.y1, F.h / 50);
    ctx.lineWidth = 1;
    xt.forEach(function (t) {
      var x = Math.round(self.tx(t.v, F, V)) + 0.5;
      ctx.strokeStyle = grid; ctx.globalAlpha = t.major ? 1 : 0.45;
      ctx.beginPath(); ctx.moveTo(x, F.t); ctx.lineTo(x, F.t + F.h); ctx.stroke();
      ctx.globalAlpha = 1;
      if (t.major || t.label) { ctx.fillStyle = muted; ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(fmt(t.v), x, F.t + F.h + 5); }
    });
    yt.forEach(function (t) {
      var y = Math.round(self.ty(t.v, F, V)) + 0.5;
      ctx.strokeStyle = grid; ctx.globalAlpha = t.major ? 1 : 0.45;
      ctx.beginPath(); ctx.moveTo(F.l, y); ctx.lineTo(F.l + F.w, y); ctx.stroke();
      ctx.globalAlpha = 1;
      if (t.major || t.label) { ctx.fillStyle = muted; ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(fmt(t.v), F.l - 6, y); }
    });
    ctx.strokeStyle = muted; ctx.strokeRect(F.l + 0.5, F.t + 0.5, F.w, F.h);
    ctx.fillStyle = fg; ctx.font = '12px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
    ctx.fillText(this.o.xlabel || '', F.l + F.w / 2, g.h - 4);
    ctx.save(); ctx.translate(13, F.t + F.h / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'middle';
    ctx.fillText(this.o.ylabel || '', 0, 0); ctx.restore();
    // series
    ctx.save(); ctx.beginPath(); ctx.rect(F.l, F.t, F.w, F.h); ctx.clip();
    this.series.forEach(function (s) {
      ctx.strokeStyle = s.color; ctx.fillStyle = s.color; ctx.lineWidth = s.width || 1.5;
      ctx.setLineDash(s.dash || []);
      var pen = false;
      if (s.style !== 'points') {
        ctx.beginPath();
        for (var i = 0; i < s.x.length; i++) {
          var x = s.x[i], y = s.y[i];
          if (!isFinite(y) || y === null || (self.xlog && x <= 0) || (self.ylog && y <= 0)) { pen = false; continue; }
          var px = self.tx(x, F, V), py = self.ty(y, F, V);
          if (pen) ctx.lineTo(px, py); else ctx.moveTo(px, py);
          pen = true;
        }
        ctx.stroke();
      }
      if (s.style === 'points' || s.style === 'both') {
        for (var j = 0; j < s.x.length; j++) {
          var yy = s.y[j];
          if (!isFinite(yy) || (self.xlog && s.x[j] <= 0) || (self.ylog && yy <= 0)) continue;
          ctx.beginPath(); ctx.arc(self.tx(s.x[j], F, V), self.ty(yy, F, V), s.r || 1.8, 0, 2 * Math.PI); ctx.fill();
        }
      }
    });
    ctx.setLineDash([]);
    // hover crosshair
    var m = this.mouse;
    if (m && m.x >= F.l && m.x <= F.l + F.w && m.y >= F.t && m.y <= F.t + F.h) {
      ctx.strokeStyle = muted; ctx.globalAlpha = 0.6; ctx.beginPath();
      ctx.moveTo(m.x + 0.5, F.t); ctx.lineTo(m.x + 0.5, F.t + F.h); ctx.stroke(); ctx.globalAlpha = 1;
      var p = this.inv(m.x, m.y, F, V);
      if (this.o.onHover) this.o.onHover(p.x, this.series);
    } else if (this.o.onHover) this.o.onHover(null);
    ctx.restore();
    // legend
    var leg = this.series.filter(function (s) { return s.label; });
    if (leg.length) {
      ctx.font = '11px system-ui, sans-serif';
      var wmax = 0; leg.forEach(function (s) { wmax = Math.max(wmax, ctx.measureText(s.label).width); });
      // put the legend in the corner that covers the fewest data points
      var lw = wmax + 40, lh = leg.length * 16 + 4, best = Infinity, lx, ly;
      [[F.l + F.w - lw + 6, F.t + 8], [F.l + 14, F.t + F.h - lh - 4], [F.l + 14, F.t + 8], [F.l + F.w - lw + 6, F.t + F.h - lh - 4]].forEach(function (c) {
        var n = 0;
        self.series.forEach(function (s) {
          for (var i = 0; i < s.x.length; i++) {
            if (!isFinite(s.y[i]) || (self.xlog && s.x[i] <= 0) || (self.ylog && s.y[i] <= 0)) continue;
            var px = self.tx(s.x[i], F, V), py = self.ty(s.y[i], F, V);
            if (px > c[0] - 6 && px < c[0] - 6 + lw && py > c[1] - 4 && py < c[1] - 4 + lh) n++;
          }
        });
        if (n < best) { best = n; lx = c[0]; ly = c[1]; }
      });
      ctx.fillStyle = bg; ctx.globalAlpha = 0.85; ctx.fillRect(lx - 6, ly - 4, wmax + 40, leg.length * 16 + 4); ctx.globalAlpha = 1;
      leg.forEach(function (s, i) {
        var y = ly + 8 + i * 16;
        ctx.strokeStyle = s.color; ctx.fillStyle = s.color; ctx.lineWidth = s.width || 1.5; ctx.setLineDash(s.dash || []);
        if (s.style === 'points') { ctx.beginPath(); ctx.arc(lx + 11, y, 2.5, 0, 7); ctx.fill(); }
        else { ctx.beginPath(); ctx.moveTo(lx, y); ctx.lineTo(lx + 22, y); ctx.stroke(); }
        ctx.setLineDash([]); ctx.fillStyle = fg; ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
        ctx.fillText(s.label, lx + 28, y);
      });
    }
    if (this.band) {
      ctx.strokeStyle = css('--accent', '#36c'); ctx.setLineDash([4, 3]);
      ctx.strokeRect(this.band.a.x, this.band.a.y, this.band.b.x - this.band.a.x, this.band.b.y - this.band.a.y);
      ctx.setLineDash([]);
    }
  };

  /* ---------------- ImageView ---------------- */

  /* Image with physical extent {x0, x1, y0, y1}; row 0 of the source is
     drawn at the top (y1). Zoom keeps a square aspect ratio. The source may
     be an array of canvases, finest first (an image pyramid). */
  function ImageView(canvas, opts) {
    this.c = canvas; this.o = opts || {}; this.src = null; this.ext = null; this.view = null; this.mouse = null;
    attachMouse(canvas, this);
  }
  ImageView.prototype.setImage = function (src, ext, keepView) {
    this.src = src;
    var same = this.ext && ext && this.ext.x0 === ext.x0 && this.ext.x1 === ext.x1;
    this.ext = ext;
    if (!keepView || !same) this.view = null;
    this.draw();
  };
  ImageView.prototype.resetView = function () { this.view = null; this.draw(); };
  ImageView.prototype.frame = function () {
    var ml = 56, mr = 12, mt = 10, mb = 40, s = Math.max(10, Math.min(this.c.clientWidth - ml - mr, this.c.clientHeight - mt - mb));
    return { l: ml + Math.max(0, (this.c.clientWidth - ml - mr - s) / 2), t: mt, w: s, h: s };
  };
  ImageView.prototype.toPx = function (x, y, F, V) {
    return { x: F.l + F.w * (x - V.x0) / (V.x1 - V.x0), y: F.t + F.h * (1 - (y - V.y0) / (V.y1 - V.y0)) };
  };
  ImageView.prototype.toData = function (px, py, F, V) {
    return { x: V.x0 + (px - F.l) / F.w * (V.x1 - V.x0), y: V.y0 + (1 - (py - F.t) / F.h) * (V.y1 - V.y0) };
  };
  ImageView.prototype.zoomTo = function (a, b) {
    if (!this.ext) return;
    var F = this.frame(), V = this.view || this.ext, p = this.toData(a.x, a.y, F, V), q = this.toData(b.x, b.y, F, V);
    var cx = (p.x + q.x) / 2, cy = (p.y + q.y) / 2, h = Math.max(Math.abs(p.x - q.x), Math.abs(p.y - q.y)) / 2;
    this.view = { x0: cx - h, x1: cx + h, y0: cy - h, y1: cy + h };
    this.draw();
  };
  ImageView.prototype.hover = function (p) { this.mouse = p; this.draw(); };
  ImageView.prototype.draw = function () {
    var g = setupCanvas(this.c), ctx = g.ctx, F = this.frame(), self = this;
    var fg = css('--fg', '#222'), muted = css('--muted', '#777'), bg = css('--panel', '#fff');
    ctx.fillStyle = bg; ctx.fillRect(0, 0, g.w, g.h);
    if (!this.src || !this.ext) return;
    var E = this.ext, V = this.view || E, src = this.src;
    if (Array.isArray(src)) {                       // pyramid: finest level not larger than the
      var need = F.w * (window.devicePixelRatio || 1) * (E.x1 - E.x0) / (V.x1 - V.x0);   // screen, so it is never averaged down
      var lvl = src[src.length - 1];
      for (var k = src.length - 1; k >= 0; k--) if (src[k].width <= need * 1.05) lvl = src[k];
      src = lvl;
    }
    var sw = src.width, sh = src.height;
    // source rectangle of the current view
    var sx0 = (V.x0 - E.x0) / (E.x1 - E.x0) * sw, sx1 = (V.x1 - E.x0) / (E.x1 - E.x0) * sw;
    var sy0 = (E.y1 - V.y1) / (E.y1 - E.y0) * sh, sy1 = (E.y1 - V.y0) / (E.y1 - E.y0) * sh;
    ctx.save(); ctx.beginPath(); ctx.rect(F.l, F.t, F.w, F.h); ctx.clip();
    ctx.fillStyle = '#000'; ctx.fillRect(F.l, F.t, F.w, F.h);
    ctx.imageSmoothingEnabled = (sx1 - sx0) > F.w * 1.2;   // show pixels when zoomed in
    ctx.drawImage(src, sx0, sy0, sx1 - sx0, sy1 - sy0, F.l, F.t, F.w, F.h);
    if (this.o.overlay) this.o.overlay(ctx, F, V, function (x, y) { return self.toPx(x, y, F, V); });
    ctx.restore();
    ctx.strokeStyle = muted; ctx.lineWidth = 1; ctx.strokeRect(F.l + 0.5, F.t + 0.5, F.w, F.h);
    ctx.font = '11px system-ui, sans-serif'; ctx.fillStyle = muted;
    var n = Math.max(2, F.w / 70);
    linTicks(V.x0, V.x1, n).forEach(function (t) {
      var x = self.toPx(t.v, 0, F, V).x;
      ctx.beginPath(); ctx.moveTo(x, F.t + F.h); ctx.lineTo(x, F.t + F.h + 4); ctx.stroke();
      ctx.textAlign = 'center'; ctx.textBaseline = 'top'; ctx.fillText(fmt(t.v), x, F.t + F.h + 6);
    });
    linTicks(V.y0, V.y1, n).forEach(function (t) {
      var y = self.toPx(0, t.v, F, V).y;
      ctx.beginPath(); ctx.moveTo(F.l - 4, y); ctx.lineTo(F.l, y); ctx.stroke();
      ctx.textAlign = 'right'; ctx.textBaseline = 'middle'; ctx.fillText(fmt(t.v), F.l - 6, y);
    });
    ctx.fillStyle = fg; ctx.font = '12px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
    ctx.fillText(this.o.xlabel || '', F.l + F.w / 2, g.h - 3);
    ctx.save(); ctx.translate(F.l - 44, F.t + F.h / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'middle';
    ctx.fillText(this.o.ylabel || '', 0, 0); ctx.restore();
    var m = this.mouse;
    if (this.o.onHover) {
      if (m && m.x >= F.l && m.x <= F.l + F.w && m.y >= F.t && m.y <= F.t + F.h) { var d = this.toData(m.x, m.y, F, V); this.o.onHover(d.x, d.y); }
      else this.o.onHover(null);
    }
    if (this.band) {
      ctx.strokeStyle = css('--accent', '#36c'); ctx.setLineDash([4, 3]);
      ctx.strokeRect(this.band.a.x, this.band.a.y, this.band.b.x - this.band.a.x, this.band.b.y - this.band.a.y);
      ctx.setLineDash([]);
    }
  };

  root.SAXS_plot = { LinePlot: LinePlot, ImageView: ImageView, fmt: fmt };
})(typeof self !== 'undefined' ? self : globalThis);
