// Runs tests/wasm/twebgpu.nim (built by webgpu.sh) in headless Chromium with
// WebGPU on SwiftShader: node twebgpu.mjs http://127.0.0.1:PORT/twebgpu.html
// PW: the playwright package (default ~/opt/pw/node_modules/playwright).
import os from 'os';
const pw = process.env.PW || os.homedir() + '/opt/pw/node_modules/playwright';
const { chromium } = (await import(pw + '/index.js')).default;
const browser = await chromium.launch({args: ['--no-zygote', '--enable-unsafe-webgpu', '--enable-features=Vulkan',
  '--use-vulkan=swiftshader', '--use-webgpu-adapter=swiftshader', '--disable-vulkan-surface', '--ignore-gpu-blocklist']});
const page = await browser.newPage();
const done = new Promise((resolve) => page.on('console', (m) => {
  console.log(m.text());
  if (m.text() === 'DONE' || m.text().startsWith('FAILED')) resolve();
}));
page.on('pageerror', (e) => console.log('pageerror:', e.message));
await page.goto(process.argv[2]);
await done;
await browser.close();
