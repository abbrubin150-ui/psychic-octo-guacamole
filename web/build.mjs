import { build } from "esbuild";
import { mkdir, copyFile } from "node:fs/promises";

await mkdir("app/src/main/assets/www", { recursive: true });

await build({
  entryPoints: ["web/src/main.js"],
  outfile: "app/src/main/assets/www/game.js",
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

await copyFile("web/src/index.html", "app/src/main/assets/www/index.html");
await copyFile("web/src/style.css", "app/src/main/assets/www/style.css");
