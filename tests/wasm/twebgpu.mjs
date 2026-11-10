// Runs tests/wasm/twebgpu.nim (built by webgpu.sh) in headless Chromium with
// WebGPU on SwiftShader: node twebgpu.mjs http://127.0.0.1:PORT/twebgpu.html
// PW: the playwright package (default ~/opt/pw/node_modules/playwright).
// WEBGPU=vulkan: the machine's GPU through Vulkan instead of SwiftShader
// (headless Chromium needs ANGLE on Vulkan for it, and the driver's libEGL).
import os from 'os';
const pw = process.env.PW || os.homedir() + '/opt/pw/node_modules/playwright';
const { chromium } = (await import(pw + '/index.js')).default;
const adapter = process.env.WEBGPU === 'vulkan' ?
  ['--use-vulkan=native', '--use-gl=angle', '--use-angle=vulkan', '--enable-gpu', '--ozone-platform=headless'] :
  ['--use-vulkan=swiftshader', '--use-webgpu-adapter=swiftshader'];
const browser = await chromium.launch({args: ['--no-zygote', '--enable-unsafe-webgpu', '--enable-features=Vulkan',
  ...adapter, '--disable-vulkan-surface', '--ignore-gpu-blocklist']});
const page = await browser.newPage();
const done = new Promise((resolve) => page.on('console', (m) => {
  console.log(m.text());
  if (m.text() === 'DONE' || m.text().startsWith('FAILED')) resolve();
}));
page.on('pageerror', (e) => console.log('pageerror:', e.message));
await page.goto(process.argv[2]);
await done;
await browser.close();
