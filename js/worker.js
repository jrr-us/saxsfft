/* SAXS-FFT engine: runs generation + scattering in a Web Worker built from a
 * Blob (the core parts' source is stitched together), so heavy work never
 * blocks the page and it works from file://. Falls back to the main thread.
 *
 * Request:  {id, gen: params | null, scat: options}; gen = null reuses the
 *           last structure.
 * Replies:  {id, progress: text} while working, then
 *           {id, ok: true, rho?, N, info?, scat} or {id, ok: false, error}.
 */
(function (root) {
  /* Body of the worker; ctx is the worker global (or a main-thread shim). */
  function workerMain(ctx, S) {
    var last = null;
    ctx.onmessage = function (e) {
      var m = e.data;
      function progress(t) { ctx.postMessage({ id: m.id, progress: t }); }
      try {
        var out = { id: m.id, ok: true, timings: {} }, transfer = [], t0 = Date.now();
        if (m.gen) {
          progress('Generating structure…');
          var gen = S.generate(m.gen, progress);
          out.timings.generate = Date.now() - t0;
          progress('Rasterising…');
          var rho = S.rasterize(gen, m.gen);
          last = { gen: gen, rho: rho, P: m.gen };
          out.rho = Float32Array.from(rho);
          out.N = m.gen.N;
          out.info = gen.info;
          transfer.push(out.rho.buffer);
        }
        if (!last) throw new Error('No structure generated yet');
        progress('Computing FFT…');
        var t1 = Date.now();
        out.scat = S.scatter(last.rho, last.gen, last.P, m.scat);
        out.timings.scatter = Date.now() - t1;
        transfer.push(out.scat.I2D.buffer);
        ctx.postMessage(out, transfer);
      } catch (err) {
        ctx.postMessage({ id: m.id, ok: false, error: String(err && err.message || err) });
      }
    };
  }

  /* Create the engine. onReply(msg) receives every worker message. */
  function createEngine(onReply) {
    var worker = null, shim = null, url = null;
    function start() {
      try {
        var src = 'self.SAXS = {};\n' +
          root.SAXS_PARTS.map(function (p) { return '(' + p.toString() + ')(self.SAXS);'; }).join('\n') +
          '\n(' + workerMain.toString() + ')(self, self.SAXS);';
        url = URL.createObjectURL(new Blob([src], { type: 'text/javascript' }));
        worker = new Worker(url);
        worker.onmessage = function (e) { onReply(e.data); };
        worker.onerror = function (e) { onReply({ id: -1, ok: false, error: 'Worker error: ' + e.message }); };
      } catch (err) {
        worker = null;
        var S = {};
        root.SAXS_PARTS.forEach(function (p) { p(S); });
        shim = { postMessage: function (msg) { setTimeout(function () { onReply(msg); }, 0); } };
        workerMain(shim, S);
      }
    }
    start();
    return {
      usesWorker: function () { return !!worker; },
      send: function (msg) {
        if (worker) worker.postMessage(msg);
        else setTimeout(function () { shim.onmessage({ data: msg }); }, 30);
      },
      /* Abort a running job (worker only); the next request must regenerate. */
      restart: function () {
        if (!worker) return false;
        worker.terminate(); URL.revokeObjectURL(url); start();
        return true;
      }
    };
  }

  root.SAXS_createEngine = createEngine;
})(typeof self !== 'undefined' ? self : globalThis);
