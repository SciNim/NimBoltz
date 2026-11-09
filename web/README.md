# NimBoltz in a web browser

A small page that runs NimBoltz (Magboltz 11.21 in Nim) as WebAssembly in
the browser: no server-side computation.

## Build

Needs Nim and Emscripten ([emsdk](https://emscripten.org/docs/getting_started/downloads.html)):

```sh
source /path/to/emsdk/emsdk_env.sh
web/build.sh             # → web/nimboltz.js, web/nimboltz.wasm
```

`web/build.sh cli DIR` builds the command-line program instead, to run input
decks under Node: `node DIR/magboltz.js < deck.in`.

## Run

The page needs to be served over HTTP (browsers don't start workers from
`file://`), with `.wasm` as `application/wasm`; Python's server does that:

```sh
python3 -m http.server -d web 8000    # then open http://localhost:8000
```

Any static file host works (GitHub Pages, ...): the page is `index.html`,
`style.css`, `app.js`, `worker.js`, `nimboltz.js` and `nimboltz.wasm`.

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

```js
const m = await createNimBoltz({print: console.log});
const r = JSON.parse(m.ccall('nimboltzRun', 'string', ['string'],
  [JSON.stringify({mixture: [[2, 90], [12, 10]], efield: 500, collisions: 4e6})]));
console.log(r.velocity.z);   // {value: 41.2, error: 0.097, unit: "μm•ns⁻¹"}
```

`tests/wasm/trun.js` runs it under Node; `tests/wasm/regress.sh` runs the
regression decks through the command-line build.
