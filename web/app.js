// The page of NimBoltz's web build: the form, the runs (in worker.js), the results.
'use strict';

// Formulas and names of Magboltz's gases, by gas number (worker gives the numbers).
const GAS_NAMES = {
  1: 'CF₄', 2: 'Argon', 3: 'Helium-4', 4: 'Helium-3', 5: 'Neon', 6: 'Krypton', 7: 'Xenon',
  8: 'CH₄ methane', 9: 'C₂H₆ ethane', 10: 'C₃H₈ propane', 11: 'iC₄H₁₀ isobutane',
  12: 'CO₂', 13: 'neoC₅H₁₂ neopentane', 14: 'H₂O water', 15: 'O₂ oxygen', 16: 'N₂ nitrogen',
  17: 'NO nitric oxide', 18: 'N₂O nitrous oxide', 19: 'C₂H₄ ethylene', 20: 'C₂H₂ acetylene',
  21: 'H₂ hydrogen', 22: 'D₂ deuterium', 23: 'CO carbon monoxide', 24: 'Methylal',
  25: 'DME dimethyl ether', 26: 'Reid step model', 27: 'Maxwell model', 28: 'Reid ramp model',
  29: 'C₂F₆', 30: 'SF₆', 31: 'NH₃ ammonia', 32: 'C₃H₆ propene', 33: 'cC₃H₆ cyclopropane',
  34: 'CH₃OH methanol', 35: 'C₂H₅OH ethanol', 36: 'iC₃H₇OH isopropanol', 37: 'Cesium',
  38: 'F₂ fluorine', 39: 'CS₂', 40: 'COS', 41: 'CD₄', 42: 'BF₃', 43: 'C₂H₂F₄',
  44: 'TMA trimethylamine', 45: 'C₃F₄H₂', 46: 'nC₃H₇OH n-propanol', 47: 'para-H₂',
  48: 'ortho-D₂', 49: 'C₄F₇N', 50: 'CHF₃', 51: 'CF₃Br', 52: 'C₃F₈', 53: 'O₃ ozone',
  54: 'Hg₂ mercury dimer', 55: 'H₂S', 56: 'nC₄H₁₀ n-butane', 57: 'nC₅H₁₂ n-pentane',
  58: 'N₂ (Phelps)', 59: 'GeH₄ germane', 60: 'SiH₄ silane', 61: 'CCl₄',
};
const shortName = (id) => GAS_NAMES[id].split(' ')[0];

const $ = (id) => document.getElementById(id);
let gases = [];          // [{id, name}] from the module
let worker = null;
let queue = [];          // setups still to run
let running = null;      // {setup, started, lines}
let runs = [];           // finished runs: {setup, result, ms, printout}
let latest = null;

// --- worker -------------------------------------------------------------------
function startWorker() {
  worker = new Worker('worker.js');
  worker.onmessage = (e) => {
    const m = e.data;
    if (m.type === 'ready') {
      gases = m.gases;
      $('gpuOption').disabled = !m.gpu;
      $('gpuHint').textContent = m.gpu ? `GPU: ${m.gpu}, in single precision.` :
        'The GPU needs WebGPU with WGSL pointer parameters (Chrome or Edge 123+), and a build with it.';
      if (!$('gases').children.length) { addGas(2, 90); addGas(12, 10); }
      setIdle('Ready.');
    } else if (m.type === 'line') {
      if (running) running.lines.push(m.text);
    } else if (m.type === 'result') {
      finished(m.result, m.ms);
    }
  };
  worker.onerror = (e) => { setIdle('The simulation failed to load: ' + e.message); };
}

function setIdle(text) {
  $('run').disabled = false; $('stop').disabled = true;
  status(text);
}
function status(text) { $('status').textContent = text; }

function runNext() {
  if (!queue.length) return;
  const setup = queue.shift();
  running = {setup, started: performance.now(), lines: []};
  $('run').disabled = true; $('stop').disabled = false;
  tick();
  worker.postMessage({setup});
}

function tick() {
  if (!running) return;
  const s = ((performance.now() - running.started) / 1000).toFixed(0);
  const left = queue.length;
  status(`Running E = ${running.setup.efield} V/cm${left ? ` (${left} more queued)` : ''}: ${s} s`);
  if (String(running.lines.length) !== $('printout').dataset.lines) {  // the printout so far
    showPrintout(running.lines.join('\n'));
    $('printout').dataset.lines = running.lines.length;
    $('printout').scrollTop = $('printout').scrollHeight;
  }
  setTimeout(tick, 1000);
}

function finished(result, ms) {
  const run = {setup: running.setup, result, ms, printout: running.lines.join('\n')};
  running = null;
  if (result.error) {
    queue = [];
    setIdle(`Stopped: ${result.error}`);
    showPrintout(run.printout);
    if (result.fatal) { worker.terminate(); startWorker(); }
    return;
  }
  runs.push(run);
  latest = run;
  showLatest(run);
  showRuns();
  if (queue.length) runNext();
  else setIdle(`Done in ${(ms / 1000).toFixed(1)} s.`);
}

function stop() {
  if (!running) return;
  worker.terminate();
  queue = []; running = null;
  status('Stopped.');
  $('run').disabled = true; $('stop').disabled = true;
  startWorker();
}

// --- form -----------------------------------------------------------------------
function addGas(id = 16, percent = 0) {
  const rows = $('gases');
  if (rows.children.length >= 6) return;
  const row = document.createElement('div');
  row.className = 'gasrow';
  const sel = document.createElement('select');
  sel.setAttribute('aria-label', 'Gas');
  for (const g of gases) {
    if (!GAS_NAMES[g.id]) continue;
    const o = document.createElement('option');
    o.value = g.id; o.textContent = `${GAS_NAMES[g.id]} (${g.id})`;
    sel.appendChild(o);
  }
  sel.value = id;
  const pct = document.createElement('input');
  pct.type = 'number'; pct.step = 'any'; pct.min = '0'; pct.max = '100'; pct.value = percent;
  pct.setAttribute('aria-label', 'Share (%)');
  const del = document.createElement('button');
  del.type = 'button'; del.textContent = '×'; del.title = 'Remove';
  del.onclick = () => { if (rows.children.length > 1) { row.remove(); updateTotal(); } };
  pct.oninput = updateTotal;
  row.append(sel, pct, del);
  rows.appendChild(row);
  updateTotal();
}

function mixture() {
  return [...$('gases').children].map((r) => [Number(r.children[0].value), Number(r.children[1].value)]);
}

function updateTotal() {
  const total = mixture().reduce((s, [, p]) => s + p, 0);
  $('total').textContent = `Total ${+total.toFixed(6)} %`;
  $('total').classList.toggle('bad', Math.abs(total - 100) > 1e-9);
  $('addGas').disabled = $('gases').children.length >= 6;
}

function setups() {
  const fields = $('efield').value.split(/[\s,;]+/).filter(Boolean).map(Number);
  if (!fields.length || fields.some((e) => !(e > 0))) throw new Error('electric fields must be positive numbers');
  const mix = mixture();
  const total = mix.reduce((s, [, p]) => s + p, 0);
  if (Math.abs(total - 100) > 1e-6) throw new Error(`the shares add up to ${total} %, not 100 %`);
  if (new Set(mix.map(([g]) => g)).size !== mix.length) throw new Error('a gas appears twice');
  const base = {
    mixture: mix.filter(([, p]) => p > 0),
    temperature: Number($('temperature').value) + 273.15,
    pressure: Number($('pressure').value),
    bfield: Number($('bfield').value),
    angle: Number($('angle').value),
    collisions: Math.round(Number($('collisions').value) * 1e6),
    efinal: Number($('efinal').value),
    thermal: $('thermal').checked,
    penning: $('penning').checked,
    report: true,
  };
  if ($('seed').value) base.seed = Number($('seed').value);
  if ($('walkers').value) base.walkers = $('walkers').value;
  if (!(base.pressure > 0)) throw new Error('the pressure must be positive');
  if (!(base.collisions >= 1e5)) throw new Error('at least 0.1 million collisions');
  return fields.map((efield) => ({...base, efield}));
}

$('setup').onsubmit = (e) => {
  e.preventDefault();
  try { queue = setups(); } catch (err) { status(`Not run: ${err.message}.`); return; }
  runNext();
};
$('stop').onclick = stop;
$('addGas').onclick = () => addGas();

// --- numbers ----------------------------------------------------------------------
const unitText = (u) => u.replaceAll('•', '·');

const SUP = {'-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹'};
const power = (k) => `×10${String(k).split('').map((c) => SUP[c]).join('')}`;

function fmt(m) {
  // value ± error, the error with 2 significant digits; ×10ⁿ outside 10⁻³..10⁶
  const v = m.value, e = m.error;
  if (!isFinite(v)) return String(v);
  const k = v === 0 ? 0 : Math.floor(Math.log10(Math.abs(v)));
  const scaled = k < -3 || k >= 6;
  const f = scaled ? 10 ** k : 1;
  if (!e || !isFinite(e)) return Number((v / f).toPrecision(5)) + (scaled ? power(k) : '');
  const digits = Math.min(Math.max(0, 1 - Math.floor(Math.log10(e / f))), 10);
  const text = `${(v / f).toFixed(digits)} ± ${(e / f).toFixed(digits)}`;
  return scaled ? `(${text})${power(k)}` : text;
}

function mixtureLabel(s) {
  return s.mixture.map(([g, p]) => `${shortName(g)} ${+p.toFixed(4)}`).join(' / ');
}
function conditionsLabel(s) {
  let t = `${mixtureLabel(s)}, ${+(s.temperature - 273.15).toFixed(2)} °C, ${s.pressure} Torr`;
  if (s.bfield) t += `, B = ${s.bfield} T at ${s.angle}°`;
  return t;
}

// --- latest run -------------------------------------------------------------------
function quantities(r, s) {
  // [label, measurement] of a run, as Magboltz reports them
  const q = [['Drift velocity (along E)', r.velocity.z]];
  if (s.bfield) q.push(['Drift velocity x', r.velocity.x], ['Drift velocity y (E × B)', r.velocity.y]);
  q.push(['Mean energy', r.meanEnergy],
         ['Longitudinal diffusion', r.longitudinal], ['Transverse diffusion', r.transverse]);
  if (!s.bfield || s.angle === 0) {
    q.push(['Longitudinal spread √(2D/v)', r.spread.longitudinal],
           ['Transverse spread √(2D/v)', r.spread.transverse],
           ['Characteristic energy D_L/µ', r.characteristicEnergy.longitudinal],
           ['Characteristic energy D_T/µ', r.characteristicEnergy.transverse]);
  }
  q.push(['Townsend coefficient α', r.ionisation], ['Attachment coefficient η', r.attachment]);
  return q;
}

function showLatest(run) {
  const r = run.result, s = run.setup;
  $('latest').hidden = false;
  $('latestTitle').textContent = `E = ${s.efield} V/cm: ${conditionsLabel(s)}`;
  const tiles = $('tiles');
  tiles.replaceChildren();
  for (const [label, m] of quantities(r, s)) {
    const t = document.createElement('div');
    t.className = 'tile';
    t.innerHTML = `<div class="label"></div><div class="value"></div>`;
    t.children[0].textContent = label;
    t.children[1].textContent = fmt(m) + ' ';
    const u = document.createElement('span');
    u.className = 'unit'; u.textContent = unitText(m.unit);
    t.children[1].appendChild(u);
    tiles.appendChild(t);
  }
  const w = $('warnings');
  w.replaceChildren();
  for (const text of r.warnings) {
    const li = document.createElement('li');
    li.textContent = text;
    w.appendChild(li);
  }
  $('runinfo').textContent = `Upper energy ${+r.efinal.value.toPrecision(6)} eV; ` +
    `${r.collisions.toLocaleString()} collisions (${r.nullCollisions.toLocaleString()} null) ` +
    `in ${(run.ms / 1000).toFixed(1)} s.`;
  showSpectrum(r);
  showProcesses(r);
  showPrintout(run.printout);
}

function showPrintout(text) {
  $('printPanel').hidden = false;
  $('printout').textContent = text;
}
$('savePrint').onclick = (e) => {
  e.preventDefault();
  download('magboltz.txt', $('printout').textContent, 'text/plain');
};

let showAllProcesses = false;
$('allProcesses').onclick = () => { showAllProcesses = !showAllProcesses; showProcesses(latest.result); };

function showProcesses(r) {
  $('processPanel').hidden = false;
  const f = r.frequencies;
  const u = unitText(f.total.unit);
  const g = (x) => Number(x.value.toPrecision(4));
  $('kinds').textContent = `Total ${g(f.total)} ${u}: elastic ${g(f.elastic)}, inelastic ${g(f.inelastic)}, ` +
    `ionisation ${g(f.ionisation)}, attachment ${g(f.attachment)} ${u}`;
  const rows = [...r.processes].filter((p) => p.frequency.value > 0)
    .sort((a, b) => b.frequency.value - a.frequency.value);
  const t = $('processes');
  t.replaceChildren();
  const head = t.createTHead().insertRow();
  for (const h of ['Gas', 'Process', `Frequency (${u})`, 'Share']) {
    const th = document.createElement('th'); th.textContent = h; head.appendChild(th);
  }
  const body = t.createTBody();
  const shown = showAllProcesses ? rows : rows.slice(0, 15);
  $('allProcesses').hidden = rows.length <= 15;
  $('allProcesses').textContent = showAllProcesses ? 'Show the first 15' : `Show all ${rows.length}`;
  for (const p of shown) {
    const tr = body.insertRow();
    tr.insertCell().textContent = shortName(latest.setup.mixture[p.gas - 1][0]);
    const d = tr.insertCell(); d.textContent = p.description; d.className = 'text';
    tr.insertCell().textContent = fmt(p.frequency);
    tr.insertCell().textContent = (100 * p.frequency.value / f.total.value).toFixed(2) + ' %';
  }
}

// --- all runs -------------------------------------------------------------------------
const COLUMNS = [
  ['vz', 'Drift velocity', (r) => r.velocity.z],
  ['vx', 'Drift velocity x', (r) => r.velocity.x],
  ['vy', 'Drift velocity y', (r) => r.velocity.y],
  ['e', 'Mean energy', (r) => r.meanEnergy],
  ['dl', 'D_L', (r) => r.longitudinal],
  ['dt', 'D_T', (r) => r.transverse],
  ['sl', 'σ_L √(2D_L/v)', (r) => r.spread.longitudinal],
  ['st', 'σ_T √(2D_T/v)', (r) => r.spread.transverse],
  ['a', 'α', (r) => r.ionisation],
  ['eta', 'η', (r) => r.attachment],
];
for (const [key, label] of COLUMNS) {
  const o = document.createElement('option');
  o.value = key; o.textContent = label;
  $('scanY').appendChild(o);
}
$('scanY').onchange = $('scanLogX').onchange = () => showRuns();

function showRuns() {
  $('scanPanel').hidden = !runs.length;
  if (!runs.length) return;
  const t = $('runs');
  t.replaceChildren();
  const head = t.createTHead().insertRow();
  const unitOf = (c) => unitText(c[2](runs[0].result).unit);
  const cols = COLUMNS.filter((c) => runs.some((r) => c[2](r.result).value !== 0));
  for (const h of ['Setup', 'E (V/cm)', ...cols.map((c) => `${c[1]} (${unitOf(c)})`), 'Time (s)']) {
    const th = document.createElement('th'); th.textContent = h; head.appendChild(th);
  }
  const body = t.createTBody();
  const series = seriesOf(runs);
  for (const run of runs) {
    const tr = body.insertRow();
    const c = tr.insertCell();
    const key = document.createElement('span');
    key.className = 'key';
    key.style.background = series.color(run);
    c.append(key, conditionsLabel(run.setup));
    tr.insertCell().textContent = run.setup.efield;
    const oblique = run.setup.bfield && run.setup.angle !== 0;  // no √(2D/v) then
    for (const col of cols) {
      tr.insertCell().textContent = oblique && ['sl', 'st'].includes(col[0]) ? '–' : fmt(col[2](run.result));
    }
    tr.insertCell().textContent = (run.ms / 1000).toFixed(1);
  }
  const col = COLUMNS.find((c) => c[0] === $('scanY').value);
  const groups = series.groups.map((g) => ({
    name: g.label, color: g.color,
    points: g.runs.filter((r) => !(['sl', 'st'].includes(col[0]) && r.setup.bfield && r.setup.angle !== 0))
      .map((r) => ({x: r.setup.efield, y: col[2](r.result).value, err: col[2](r.result).error, run: r}))
      .sort((a, b) => a.x - b.x),
  }));
  plot($('scanChart'), groups, {
    xlog: $('scanLogX').checked, xlabel: 'E (V/cm)',
    ylabel: `${col[1]} (${unitOf(col)})`, lines: true, errors: true,
    tip: (p) => `E = ${p.x} V/cm<br>${col[1]} = ${fmt(col[2](p.run.result))} ${unitOf(col)}`,
  });
}

function seriesOf(rs) {
  // one series per setup apart from E, in fixed colour order (3 slots, then grey)
  const slots = ['var(--series-1)', 'var(--series-2)', 'var(--series-3)'];
  const groups = [];
  const keyOf = (r) => conditionsLabel(r.setup) + JSON.stringify([r.setup.collisions, r.setup.thermal, r.setup.penning, r.setup.walkers]);
  for (const r of rs) {
    let g = groups.find((g) => g.key === keyOf(r));
    if (!g) {
      g = {key: keyOf(r), label: conditionsLabel(r.setup), runs: [],
           color: slots[groups.length] || 'var(--text-muted)'};
      groups.push(g);
    }
    g.runs.push(r);
  }
  return {groups, color: (r) => groups.find((g) => g.key === keyOf(r)).color};
}

$('clear').onclick = () => { runs = []; showRuns(); };
$('csv').onclick = () => {
  const head = ['mixture', 'temperature_K', 'pressure_Torr', 'E_V/cm', 'B_T', 'angle_deg', 'collisions'];
  for (const [key, label, get] of COLUMNS) {
    const u = get(runs[0].result).unit;
    head.push(`${key}_${u}`, `${key}_err`);
  }
  const lines = [head.join(',')];
  for (const {setup: s, result: r} of runs) {
    const row = [`"${mixtureLabel(s)}"`, s.temperature, s.pressure, s.efield, s.bfield, s.angle, r.collisions];
    for (const c of COLUMNS) { const m = c[2](r); row.push(m.value, m.error); }
    lines.push(row.join(','));
  }
  download('nimboltz-runs.csv', lines.join('\n') + '\n', 'text/csv');
};

function download(name, text, type) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([text], {type}));
  a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

// --- energy distribution ----------------------------------------------------------
function showSpectrum(r) {
  // the energy bins merged to ≤ 200, as a density: fraction of the collisions per eV
  $('spectrumPanel').hidden = false;
  const {energy, fraction} = r.spectrum;
  const n = energy.length, merge = Math.max(1, Math.ceil(n / 200));
  const width = n > 1 ? energy[1] - energy[0] : 1;
  const pts = [];
  for (let i = 0; i < n; i += merge) {
    const k = Math.min(merge, n - i);
    let sum = 0;
    for (let j = i; j < i + k; j++) sum += fraction[j];
    pts.push({x: energy[i] + (k - 1) * width / 2, y: sum / (k * width), share: sum, lo: energy[i] - width / 2, hi: energy[i + k - 1] + width / 2});
  }
  const ylog = $('specLogY').checked;
  plot($('spectrum'), [{name: 'collisions', color: 'var(--series-1)', points: ylog ? pts.filter((p) => p.y > 0) : pts}], {
    xlabel: 'Energy (eV)', ylabel: 'Fraction of the collisions per eV', ylog, area: !ylog, lines: true,
    tip: (p) => `${p.lo.toPrecision(3)}–${p.hi.toPrecision(3)} eV<br>${(100 * p.share).toPrecision(3)} % of the collisions`,
  });
}
$('specLogY').onchange = () => latest && showSpectrum(latest.result);

// --- plots (SVG) -------------------------------------------------------------------
function niceTicks(lo, hi, n = 6) {
  const span = hi - lo || Math.abs(hi) || 1;
  const step0 = span / n, mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map((f) => f * mag).find((s) => s >= step0);
  const ticks = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step) ticks.push(+t.toPrecision(12));
  return ticks;
}
function logTicks(lo, hi) {
  const ticks = [];
  for (let k = Math.floor(Math.log10(lo)); k <= Math.ceil(Math.log10(hi)); k++) {
    const t = 10 ** k;
    if (t >= lo * 0.999 && t <= hi * 1.001) ticks.push(t);
  }
  if (ticks.length < 3) for (const f of [2, 3, 5]) {
    for (let k = Math.floor(Math.log10(lo)); k <= Math.ceil(Math.log10(hi)); k++) {
      const t = f * 10 ** k;
      if (t >= lo && t <= hi) ticks.push(t);
    }
  }
  return ticks.sort((a, b) => a - b);
}
const tickText = (t) => Math.abs(t) >= 1e5 || (Math.abs(t) < 1e-3 && t !== 0) ? t.toExponential(0) : String(+t.toPrecision(6));

function plot(svg, series, o) {
  const NS = 'http://www.w3.org/2000/svg';
  const W = svg.clientWidth || 800, H = svg.clientHeight || 300;
  const m = {l: 64, r: 16, t: series.length > 1 ? 30 : 12, b: 40};
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  svg.replaceChildren();
  $('tooltip').hidden = true;
  const el = (tag, attrs, parent = svg) => {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    parent.appendChild(e);
    return e;
  };
  const all = series.flatMap((s) => s.points);
  if (!all.length) return;
  let [x0, x1] = [Math.min(...all.map((p) => p.x)), Math.max(...all.map((p) => p.x))];
  let [y0, y1] = [Math.min(...all.map((p) => p.y - (o.errors ? p.err : 0))), Math.max(...all.map((p) => p.y + (o.errors ? p.err : 0)))];
  if (o.xlog) { x0 = Math.max(x0, 1e-300); if (x0 === x1) { x0 /= 2; x1 *= 2; } }
  else { if (x0 === x1) { x0 -= Math.abs(x0) * 0.5 || 1; x1 += Math.abs(x1) * 0.5 || 1; } const p = (x1 - x0) * 0.03; x0 -= p; x1 += p; if (all.every((q) => q.x >= 0)) x0 = Math.max(x0, 0); }
  if (o.ylog) { y0 = Math.min(...all.filter((p) => p.y > 0).map((p) => p.y)); if (y0 === y1) { y0 /= 2; y1 *= 2; } }
  else { if (y0 > 0 && y0 < (y1 - y0)) y0 = 0; if (y0 === y1) { y0 -= Math.abs(y0) * 0.5 || 1; y1 += Math.abs(y1) * 0.5 || 1; } y1 += (y1 - y0) * 0.05; }
  const fx = o.xlog ? (x) => m.l + (Math.log10(x) - Math.log10(x0)) / (Math.log10(x1) - Math.log10(x0)) * (W - m.l - m.r)
                    : (x) => m.l + (x - x0) / (x1 - x0) * (W - m.l - m.r);
  const fy = o.ylog ? (y) => H - m.b - (Math.log10(y) - Math.log10(y0)) / (Math.log10(y1) - Math.log10(y0)) * (H - m.t - m.b)
                    : (y) => H - m.b - (y - y0) / (y1 - y0) * (H - m.t - m.b);
  // grid and axes
  const xt = o.xlog ? logTicks(x0, x1) : niceTicks(x0, x1);
  const yt = o.ylog ? logTicks(y0, y1) : niceTicks(y0, y1, 5);
  for (const t of yt) {
    el('line', {x1: m.l, x2: W - m.r, y1: fy(t), y2: fy(t), class: 'gridline'});
    el('text', {x: m.l - 6, y: fy(t) + 4, 'text-anchor': 'end'}).textContent = tickText(t);
  }
  for (const t of xt) {
    el('line', {x1: fx(t), x2: fx(t), y1: H - m.b, y2: H - m.b + 4, class: 'axis'});
    el('text', {x: fx(t), y: H - m.b + 16, 'text-anchor': 'middle'}).textContent = tickText(t);
  }
  el('line', {x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: 'axis'});
  el('text', {x: (m.l + W - m.r) / 2, y: H - 6, 'text-anchor': 'middle'}).textContent = o.xlabel;
  el('text', {x: 12, y: (m.t + H - m.b) / 2, 'text-anchor': 'middle', transform: `rotate(-90 12 ${(m.t + H - m.b) / 2})`}).textContent = o.ylabel;
  // legend (two or more series)
  if (series.length > 1) {
    let lx = m.l;
    for (const s of series) {
      el('circle', {cx: lx + 5, cy: 12, r: 5, fill: s.color});
      const t = el('text', {x: lx + 14, y: 16});
      t.textContent = s.name;
      lx += 24 + t.getComputedTextLength();
    }
  }
  // marks
  for (const s of series) {
    const pts = s.points;
    const d = pts.map((p, i) => `${i ? 'L' : 'M'}${fx(p.x).toFixed(1)},${fy(p.y).toFixed(1)}`).join('');
    if (o.area && pts.length > 1) {
      el('path', {d: `${d}L${fx(pts.at(-1).x).toFixed(1)},${fy(Math.max(y0, 0))}L${fx(pts[0].x).toFixed(1)},${fy(Math.max(y0, 0))}Z`,
                  fill: s.color, 'fill-opacity': 0.1, stroke: 'none'});
    }
    if (o.lines && pts.length > 1) {
      el('path', {d, fill: 'none', stroke: s.color, 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round'});
    }
    if (o.errors) {
      for (const p of pts) {
        if (p.err > 0) el('line', {x1: fx(p.x), x2: fx(p.x), y1: fy(p.y - p.err), y2: fy(p.y + p.err), stroke: s.color, 'stroke-width': 1.5});
        el('circle', {cx: fx(p.x), cy: fy(p.y), r: 4, fill: s.color, stroke: 'var(--surface-1)', 'stroke-width': 2});
      }
    }
  }
  // hover: crosshair and tooltip at the nearest point
  const cross = el('line', {y1: m.t, y2: H - m.b, class: 'crosshair', visibility: 'hidden'});
  const dot = el('circle', {r: 5, stroke: 'var(--surface-1)', 'stroke-width': 2, visibility: 'hidden'});
  const tip = $('tooltip');
  svg.onmousemove = (e) => {
    const b = svg.getBoundingClientRect();
    const mx = (e.clientX - b.left) * W / b.width, my = (e.clientY - b.top) * H / b.height;
    let best = null, bd = Infinity;
    for (const s of series) for (const p of s.points) {
      const dx = fx(p.x) - mx, dy = series.length > 1 || o.errors ? fy(p.y) - my : 0;
      const dd = dx * dx + dy * dy;
      if (dd < bd) { bd = dd; best = {p, s}; }
    }
    if (!best) return;
    cross.setAttribute('x1', fx(best.p.x)); cross.setAttribute('x2', fx(best.p.x));
    cross.setAttribute('visibility', 'visible');
    dot.setAttribute('cx', fx(best.p.x)); dot.setAttribute('cy', fy(best.p.y));
    dot.setAttribute('fill', best.s.color); dot.setAttribute('visibility', 'visible');
    tip.innerHTML = o.tip(best.p);
    tip.hidden = false;
    tip.style.left = `${e.clientX + 14}px`; tip.style.top = `${e.clientY + 14}px`;
  };
  svg.onmouseleave = () => { tip.hidden = true; cross.setAttribute('visibility', 'hidden'); dot.setAttribute('visibility', 'hidden'); };
}
window.addEventListener('resize', () => { showRuns(); if (latest) showSpectrum(latest.result); });

startWorker();
