// Runs NimBoltz (nimboltz.js/.wasm) off the page's thread. Messages:
//   in:  {setup}                        one run (see src/nimboltz/wasm.nim)
//   out: {type: 'ready', gases, scale}  once the module is loaded
//        {type: 'line', text}           a line of Magboltz's printout
//        {type: 'result', result, ms}   the run's results (or {error})
// A run can't be interrupted from inside: the page terminates the worker.
importScripts('nimboltz.js');

const ready = createNimBoltz({
  print: (text) => postMessage({type: 'line', text}),
  printErr: (text) => postMessage({type: 'line', text}),
}).then((m) => {
  const gases = JSON.parse(m.UTF8ToString(m.ccall('nimboltzGases', 'number', [], [])));
  postMessage({type: 'ready', gases, scale: m.ccall('nimboltzCollisionScale', 'number', [], [])});
  return m;
});

onmessage = async (e) => {
  const m = await ready;
  const t = performance.now();
  let result;
  try {
    result = JSON.parse(m.ccall('nimboltzRun', 'string', ['string'], [JSON.stringify(e.data.setup)]));
  } catch (err) {  // a Nim defect or a trap: the module is unusable now
    result = {error: String(err), fatal: true};
  }
  postMessage({type: 'result', result, ms: performance.now() - t});
};
