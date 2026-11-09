// Runs the web module (web/nimboltz.js) under Node, as web/worker.js does:
//   node tests/wasm/trun.js '{"mixture": [[2, 90], [12, 10]], "efield": 500}'
const path = require('path');
const lines = [];
const createNimBoltz = require(path.join(__dirname, '../../web/nimboltz.js'));
createNimBoltz({
  print: (s) => lines.push(s),
  printErr: (s) => console.error(s),
}).then((Module) => {
  const gases = JSON.parse(Module.UTF8ToString(Module.ccall('nimboltzGases', 'number', [], [])));
  console.log(`${gases.length} gases, ${gases[0].id}: ${gases[0].name} .. ${gases.at(-1).id}: ${gases.at(-1).name}`);
  console.log('collisions per NMAX', Module.ccall('nimboltzCollisionScale', 'number', [], []));
  const setup = process.argv[2] || '{"mixture": [[2, 90], [12, 10]], "efield": 500, "report": true}';
  const t = Date.now();
  const r = JSON.parse(Module.ccall('nimboltzRun', 'string', ['string'], [setup]));
  console.log(`${Date.now() - t} ms, ${lines.length} lines printed`);
  if (r.error) { console.log('error:', r.error); return; }
  const m = (x) => `${x.value.toPrecision(5)} ± ${x.error.toPrecision(2)} ${x.unit}`;
  console.log('vz', m(r.velocity.z), ' meanE', m(r.meanEnergy));
  console.log('DL', m(r.longitudinal), ' DT', m(r.transverse));
  console.log('alpha', m(r.ionisation), ' eta', m(r.attachment));
  console.log('spectrum bins', r.spectrum.energy.length, 'processes', r.processes.length,
              'warnings', r.warnings);
});
