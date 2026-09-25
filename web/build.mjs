import { build } from "esbuild";
import { access, copyFile, mkdir } from "node:fs/promises";

const outDir = "../app/src/main/assets/www";
const vendorDir = outDir + "/vendor";
await mkdir(vendorDir, { recursive: true });

await build({
  entryPoints: ["src/main.js"],
  outfile: outDir + "/game.js",
  bundle: true,
  minify: false,
  sourcemap: false,
  format: "esm",
  platform: "browser",
  target: ["chrome120"],
  external: ["@dimforge/rapier3d-compat"],
  define: {
    "process.env.NODE_ENV": "\"production\""
  },
  logLevel: "info"
});

const rapierRoot = "node_modules/@dimforge/rapier3d-compat";
let rapierModule = null;
for (const candidate of ["rapier.mjs", "rapier.es.js", "rapier.js"]) {
  try {
    await access(rapierRoot + "/" + candidate);
    rapierModule = candidate;
    break;
  } catch {}
}
if (!rapierModule) {
  throw new Error("Rapier ESM runtime not found in installed package");
}

await copyFile("src/index.html", outDir + "/index.html");
await copyFile("src/style.css", outDir + "/style.css");
await copyFile(
  rapierRoot + "/" + rapierModule,
  vendorDir + "/rapier.mjs"
);
await copyFile(
  rapierRoot + "/rapier_wasm3d_bg.wasm",
  vendorDir + "/rapier_wasm3d_bg.wasm"
);
