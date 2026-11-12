# NimBoltz in a web browser

A small page that runs NimBoltz (Magboltz 11.21 in Nim) as WebAssembly in
the browser: no server-side computation.

## Build

Needs Nim and Emscripten ([emsdk](https://emscripten.org/docs/getting_started/downloads.html)):

```sh
source /path/to/emsdk/emsdk_env.sh
web/build.sh             # → web/nimboltz.{js,wasm}, web/nimboltz-mt.{js,wasm}
```

`web/build.sh cli DIR` builds the command-line program instead, to run input
decks under Node: `node DIR/magboltz.js < deck.in`.

With the GPU: the drift run can also run on independent electrons on the
GPU through WebGPU. That build needs [tattletale](https://github.com/mratsim/tattletale)
(its Nim-to-GPU compiler Crucible, with the WGSL changes of the local
branch `nimboltz`):

```sh
TATTLETALE=/path/to/tattletale web/build.sh
```

## Run

The page needs to be served over HTTP (browsers don't start workers from
`file://`), with `.wasm` as `application/wasm`; Python's server does that:

```sh
python3 -m http.server -d web 8000    # then open http://localhost:8000
```

Any static file host works (GitHub Pages, ...): the page is `index.html`,
`style.css`, `app.js`, `worker.js`, `nimboltz.js` and `nimboltz.wasm`,
and for threads `nimboltz-mt.js` and `nimboltz-mt.wasm`.

Threads: `nimboltz-mt` runs on one thread per core (Emscripten's pthreads,
Web Workers that share the module's memory). Browsers allow that shared
memory (`SharedArrayBuffer`) only on pages with cross-origin isolation,
which the server turns on with two headers:

```
Cross-Origin-Opener-Policy: same-origin
Cross-Origin-Embedder-Policy: require-corp
```

`worker.js` loads `nimboltz-mt` when the page has them
(`crossOriginIsolated`), else `nimboltz` on one thread. `web/serve.py`
is Python's server with the headers:

```sh
python3 web/serve.py 8000             # then open http://localhost:8000
```

With nginx:

```nginx
location /nimboltz/ {
    add_header Cross-Origin-Opener-Policy  same-origin   always;
    add_header Cross-Origin-Embedder-Policy require-corp always;
}
```

(A `location` with `add_header` of its own inherits none from the
`server` block: repeat those there.) The page loads only its own files, so
the headers block nothing. GitHub Pages can't set headers.

## The page

- Up to 6 gases with their shares, temperature, pressure, the electric
  field (several fields separated by spaces run one after the other), the
  magnetic field and its angle to E.
- Collisions: real collisions of the drift run, in millions. Magboltz's
  NMAX = 2 is 80 million; a few million give rough numbers in seconds.
  Upper energy 0 finds the energy range as Magboltz does.
- The latest run: drift velocity, mean energy, diffusion, √(2D/v), D/µ,
  Townsend and attachment coefficients, all with Magboltz's errors; the
  energy distribution of the collisions; the collision frequencies by
  process; Magboltz's printout (also while running).
- All runs: a table (CSV download) and a plot of any result against E,
  one line per setup.

- Drift run: independent electrons (walkers, the default: Magboltz's
  results within errors, faster), on the CPU or on the GPU, or Magboltz's
  one electron, bit for bit as Magboltz 11.21. The GPU needs a browser with WebGPU and WGSL pointer parameters
  (Chrome or Edge 123 and later) and the build with it; it runs in single
  precision, 16384 walkers, and its runs count more collisions than asked
  (each walker warms up first). When the browser can't use the graphics
  card it hands out a software adapter (SwiftShader), which runs on the CPU;
  the page says so under the selector, and chrome://gpu tells why.
- Threads (with the headers above): up to one per core, for many electrons
  on the CPU (16 threads: about 9 times faster, the same numbers as one
  thread, since each electron has its own random numbers) and for the
  avalanches of the Townsend analyses (their own random numbers: not
  Magboltz's numbers, the same within errors). Magboltz's one electron
  stays on one thread.

Runs happen in a Web Worker; Stop ends the worker (and the run) and starts
a new one. Runs of the Townsend analyses (high fields, ionisation or
attachment) take much longer than the drift run itself.

## The module

`nimboltz.js` defines `createNimBoltz(options)` (Emscripten `MODULARIZE`)
with three C functions (`src/nimboltz/wasm.nim`):

- `nimboltzRun(setupJson)` → results as JSON (or `{"error": ...}`), with
  Magboltz's printout to `options.print` if the setup has `"report": true`
- `nimboltzGases()` → `[{"id", "name"}, ...]`
- `nimboltzCollisionScale()` → collisions per unit of Magboltz's NMAX
- `nimboltzWebGpu()` → 1 in the build with WebGPU; then set
  `m.gpuDevice` to a `GPUDevice` before a run with `"walkers": "gpu"`,
  and call `nimboltzRun` with `ccall(..., {async: true})` (it waits for the
  GPU), as `worker.js` does

```js
const m = await createNimBoltz({print: console.log});
const r = JSON.parse(m.ccall('nimboltzRun', 'string', ['string'],
  [JSON.stringify({mixture: [[2, 90], [12, 10]], efield: 500, collisions: 4e6})]));
console.log(r.velocity.z);   // {value: 41.2, error: 0.097, unit: "μm•ns⁻¹"}
```

`tests/wasm/trun.js` runs it under Node; `tests/wasm/regress.sh` runs the
regression decks through the command-line build; `tests/wasm/webgpu.sh`
compares WebGPU walkers with CPU walkers collision by collision in headless
Chromium (WebGPU on SwiftShader, no GPU needed; `WEBGPU=vulkan` for the
machine's GPU through Vulkan).
