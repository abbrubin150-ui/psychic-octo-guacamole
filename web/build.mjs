import { build } from "esbuild";
import { mkdir, copyFile } from "node:fs/promises";

const outDir = "../app/src/main/assets/www";
await mkdir(outDir, { recursive: true });

await build({
  entryPoints: ["src/main.js"],
  outfile: outDir + "/game.js",
  bundle: true,
  minify: true,
  sourcemap: true,
  format: "iife",
  platform: "browser",
  target: ["chrome120"],
  define: {
    "process.env.NODE_ENV": "\"production\""
  },
  logLevel: "info"
});

await copyFile("src/index.html", outDir + "/index.html");
await copyFile("src/style.css", outDir + "/style.css");
await copyFile(
  "node_modules/@dimforge/rapier3d-compat/rapier_wasm3d_bg.wasm",
  outDir + "/rapier_wasm3d_bg.wasm"
);
