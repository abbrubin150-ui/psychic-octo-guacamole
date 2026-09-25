import * as THREE from "three";
import { PhysicsKernel, ContactPointGrabber, RAPIER } from "./physics.js";
import { PointerFusion } from "./input.js";
import { Telemetry } from "./telemetry.js";
import { RENDER, SAVE_SCHEMA } from "./config.js";
import { OrthoRig, PixelPipeline } from "./view.js";

const bootEl = document.querySelector("#boot");

function bootMessage(message) {
  if (!bootEl) return;
  bootEl.hidden = false;
  const msg = bootEl.querySelector("[data-role=message]");
  if (msg) msg.textContent = message;
}

function bootReady() {
  if (bootEl) bootEl.hidden = true;
  globalThis.__PIXEL_PHYSICS_LAB_READY__ = true;
  globalThis.AndroidBridge?.log?.("info", "APP_READY Pixel Physics Lab 0.8.2");
}

function bootFatal(error) {
  const message = String(error?.stack || error?.message || error || "Unknown startup failure");
  globalThis.__PIXEL_PHYSICS_LAB_READY__ = false;
  globalThis.AndroidBridge?.log?.("error", "BOOT_FATAL " + message);
  if (!bootEl) return;
  bootEl.hidden = false;
  bootEl.classList.add("fatal");
  const title = bootEl.querySelector("[data-role=title]");
  const msg = bootEl.querySelector("[data-role=message]");
  if (title) title.textContent = "STARTUP ERROR";
  if (msg) msg.textContent = message.slice(0, 420);
}

async function bootstrap() {
  bootMessage("INITIALIZING RENDERER…");
const canvas = document.querySelector("#game");
const statusEl = document.querySelector("#status");
const contextEl = document.querySelector("#context");
const telemetry = new Telemetry(globalThis.AndroidBridge);
telemetry.installGlobalFaultHooks();

bootMessage("STARTING WEBGL…");
const renderer = new THREE.WebGLRenderer({
  canvas,
  antialias: false,
  alpha: false,
  powerPreference: "high-performance",
  premultipliedAlpha: false
});
renderer.setPixelRatio(1);
renderer.setSize(RENDER.width, RENDER.height, false);
renderer.outputColorSpace = THREE.SRGBColorSpace;
renderer.shadowMap.enabled = true;
renderer.shadowMap.type = THREE.PCFShadowMap;

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x101721);

const camera = new THREE.OrthographicCamera(-7.2, 7.2, 4.05, -4.05, 0.05, 80);
const rig = new OrthoRig(camera);

const hemi = new THREE.HemisphereLight(0xffddb0, 0x273247, 1.35);
scene.add(hemi);

const sun = new THREE.DirectionalLight(0xffc27a, 2.1);
sun.position.set(-4, 8, 5);
sun.castShadow = true;
sun.shadow.mapSize.set(1024, 1024);
sun.shadow.camera.left = -7;
sun.shadow.camera.right = 7;
sun.shadow.camera.top = 7;
sun.shadow.camera.bottom = -7;
sun.shadow.bias = -0.0008;
scene.add(sun);

const pipeline = new PixelPipeline(renderer, scene, camera);
bootMessage("LOADING PHYSICS CORE…");
const physics = await new PhysicsKernel(telemetry).init();
const grabber = new ContactPointGrabber(physics, telemetry);
const raycaster = new THREE.Raycaster();
const dynamicRoots = [];
const spawnRecords = new Map();

bootMessage("BUILDING WORKSHOP…");
buildWorkshop();
await restoreOrCreateWorld();

const input = new PointerFusion(canvas, {
  pick: (x, y, radius) => pickFatFinger(x, y, radius),
  beginGrab: (hit, x, y, timeSec) => {
    selected = hit.entity;
    const ndc = clientToNdc(x, y);
    grabber.start(hit.entity, hit.point, camera, ndc.x, ndc.y, timeSec);
    haptic(12, 105);
    telemetry.count("input.pickup");
  },
  moveGrab: (x, y, timeSec) => {
    const ndc = clientToNdc(x, y);
    grabber.setScreen(ndc.x, ndc.y, camera, timeSec);
  },
  gestureGrab: ({ x, y, depthDeltaPx, twistDelta, dt, timeSec }) => {
    const ndc = clientToNdc(x, y);
    grabber.adjustDepth(depthDeltaPx, camera, ndc.x, ndc.y, timeSec);
    grabber.addTwist(-twistDelta, dt);
  },
  endGrab: () => {
    grabber.release();
    scheduleSave();
  },
  isGrabbing: () => !!grabber.active,
  orbit: (dx, dy) => rig.orbit(dx, dy),
  cameraGesture: ({ panX, panY, pinchDeltaPx, twistDelta }) => {
    rig.pan(panX, panY);
    rig.zoomByPixels(pinchDeltaPx);
    rig.orbit(-twistDelta * 130, 0);
  },
  tap: hit => {
    selected = hit?.entity || null;
  },
  doubleTap: hit => {
    if (hit?.entity) rig.focusEntity(hit.entity);
    else rig.reset();
    haptic(14, 95);
  },
  longPress: (hit, x, y) => {
    if (!hit?.entity) return;
    selected = hit.entity;
    openContext(hit.entity, x, y);
    haptic(22, 135);
  }
}, telemetry);

let selected = null;
let running = true;
let lastFrame = performance.now();
let saveTimer = null;
let fpsAccum = 0;
let fpsFrames = 0;
let fps = 0;

rig.update(0);
pipeline.render();
bootReady();
requestAnimationFrame(frame);

function frame(now) {
  requestAnimationFrame(frame);
  if (!running) {
    lastFrame = now;
    return;
  }

  const start = performance.now();
  const dt = Math.min(0.1, Math.max(0, (now - lastFrame) / 1000));
  lastFrame = now;

  const alpha = physics.stepFrame(dt, stepDt => {
    grabber.preStep(stepDt);
  });
  physics.syncVisuals(alpha);

  rig.update(dt);
  pipeline.render();

  fpsAccum += dt;
  fpsFrames++;
  if (fpsAccum >= 0.5) {
    fps = Math.round(fpsFrames / fpsAccum);
    fpsAccum = 0;
    fpsFrames = 0;
  }

  updateStatus();
  telemetry.frame(performance.now() - start);
  invariantSweep();
}

function buildWorkshop() {
  const floorMat = pixelMaterial(0x9b6a47);
  const wallMat = pixelMaterial(0x57433d);
  const trimMat = pixelMaterial(0x6a2f1d);

  const floor = meshBox(9.6, 0.2, 9.6, floorMat);
  floor.position.set(0, -0.1, 0);
  floor.receiveShadow = true;
  scene.add(floor);
  physics.createStaticFloor(floor, [4.8, 0.10, 4.8], [0, -0.10, 0]);

  const wallBack = meshBox(9.6, 3.6, 0.22, wallMat);
  wallBack.position.set(0, 1.8, -4.7);
  wallBack.receiveShadow = true;
  scene.add(wallBack);
  physics.createStaticBox({
    mesh: wallBack,
    position: wallBack.position,
    halfExtents: new THREE.Vector3(4.8, 1.8, 0.11)
  });

  const wallSide = meshBox(0.22, 3.6, 9.6, wallMat);
  wallSide.position.set(-4.7, 1.8, 0);
  wallSide.receiveShadow = true;
  scene.add(wallSide);
  physics.createStaticBox({
    mesh: wallSide,
    position: wallSide.position,
    halfExtents: new THREE.Vector3(0.11, 1.8, 4.8)
  });

  for (const [pos, size] of [
    [[0, 0.14, -4.52], [9.45, 0.18, 0.18]],
    [[-4.52, 0.14, 0], [0.18, 0.18, 9.45]],
    [[0, 3.42, -4.52], [9.45, 0.22, 0.22]],
    [[-4.52, 3.42, 0], [0.22, 0.22, 9.45]]
  ]) {
    const trim = meshBox(size[0], size[1], size[2], trimMat);
    trim.position.set(...pos);
    scene.add(trim);
  }

  const bench = meshBox(2.8, 0.22, 1.35, pixelMaterial(0xa85d31));
  bench.position.set(-1.8, 0.95, -2.6);
  bench.castShadow = true;
  bench.receiveShadow = true;
  scene.add(bench);
  physics.createStaticBox({
    mesh: bench,
    position: bench.position,
    halfExtents: new THREE.Vector3(1.4, 0.11, 0.675),
    friction: 0.74
  });

  const shelf = meshBox(1.5, 0.14, 0.55, pixelMaterial(0xb86a36));
  shelf.position.set(2.5, 2.1, -4.22);
  shelf.castShadow = true;
  scene.add(shelf);
  physics.createStaticBox({
    mesh: shelf,
    position: shelf.position,
    halfExtents: new THREE.Vector3(0.75, 0.07, 0.275)
  });
}

async function restoreOrCreateWorld() {
  const raw = globalThis.AndroidBridge?.loadState?.("autosave") || "";
  if (raw) {
    try {
      const save = JSON.parse(raw);
      if (save.schema === SAVE_SCHEMA && Array.isArray(save.objects)) {
        for (const spec of save.objects) spawnFromSpec(spec);
        if (save.camera) rig.restore(save.camera);
        return;
      }
    } catch (err) {
      telemetry.fault("save_restore_failed", { message: String(err) });
    }
  }

  spawnFromSpec({
    kind: "crate",
    material: "wood",
    size: [0.62, 0.62, 0.62],
    position: [-1.9, 1.4, -0.6],
    rotation: [0, 0, 0, 1]
  });
  spawnFromSpec({
    kind: "beam",
    material: "wood",
    size: [1.65, 0.24, 0.30],
    position: [-0.1, 1.2, 0.4],
    rotation: [0, 0, 0, 1]
  });
  spawnFromSpec({
    kind: "ball",
    material: "rubber",
    radius: 0.28,
    position: [1.25, 1.5, 0.2],
    rotation: [0, 0, 0, 1]
  });
  spawnFromSpec({
    kind: "cube",
    material: "metal",
    size: [0.38, 0.38, 0.38],
    position: [2.1, 1.2, -0.6],
    rotation: [0, 0, 0, 1]
  });
  spawnFromSpec({
    kind: "cube",
    material: "wood",
    size: [0.46, 0.46, 0.46],
    position: [2.4, 0.4, 1.3],
    rotation: [0, 0, 0, 1]
  });
}

function spawnFromSpec(spec) {
  const material = spec.material || "wood";
  const position = new THREE.Vector3(...(spec.position || [0, 1, 0]));
  const qArr = spec.rotation || [0, 0, 0, 1];
  const rotation = new THREE.Quaternion(qArr[0], qArr[1], qArr[2], qArr[3]);

  let mesh;
  let entity;
  if (spec.kind === "ball") {
    const radius = spec.radius || 0.25;
    mesh = new THREE.Mesh(
      new THREE.SphereGeometry(radius, 14, 9),
      pixelMaterial(materialColor(material), true)
    );
    mesh.position.copy(position);
    mesh.quaternion.copy(rotation);
    mesh.castShadow = true;
    mesh.receiveShadow = true;
    scene.add(mesh);

    entity = physics.spawnBall({ mesh, position, radius, material });
  } else {
    const size = spec.size || [0.5, 0.5, 0.5];
    mesh = makeVoxelLikeBox(size, material, spec.kind === "crate");
    mesh.position.copy(position);
    mesh.quaternion.copy(rotation);
    scene.add(mesh);

    entity = physics.spawnBox({
      mesh,
      position,
      rotation,
      halfExtents: new THREE.Vector3(size[0] / 2, size[1] / 2, size[2] / 2),
      material
    });
  }

  entity.spec = structuredCloneSafe(spec);
  entity.spec.rotation = qArr.slice();
  spawnRecords.set(entity.id, entity.spec);
  dynamicRoots.push(mesh);
  return entity;
}

function makeVoxelLikeBox(size, material, crate = false) {
  const group = new THREE.Group();
  const base = new THREE.Mesh(
    new THREE.BoxGeometry(size[0], size[1], size[2]),
    pixelMaterial(materialColor(material))
  );
  base.castShadow = true;
  base.receiveShadow = true;
  group.add(base);

  if (crate) {
    const slatMat = pixelMaterial(darken(materialColor(material), 0.64));
    const t = Math.min(size[0], size[1], size[2]) * 0.08;
    for (const y of [-size[1] * 0.34, size[1] * 0.34]) {
      const slat = new THREE.Mesh(new THREE.BoxGeometry(size[0] * 1.01, t, t), slatMat);
      slat.position.set(0, y, size[2] * 0.505);
      group.add(slat);
    }
    for (const x of [-size[0] * 0.34, size[0] * 0.34]) {
      const slat = new THREE.Mesh(new THREE.BoxGeometry(t, size[1] * 1.01, t), slatMat);
      slat.position.set(x, 0, size[2] * 0.505);
      group.add(slat);
    }
  }
  return group;
}

function pixelMaterial(color, emissive = false) {
  return new THREE.MeshStandardMaterial({
    color,
    roughness: 0.86,
    metalness: color === 0x6f7d88 ? 0.58 : 0.05,
    flatShading: true,
    emissive: emissive ? new THREE.Color(color).multiplyScalar(0.04) : new THREE.Color(0x000000)
  });
}

function materialColor(material) {
  if (material === "metal") return 0x6f7d88;
  if (material === "rubber") return 0x44624d;
  return 0x9a4a28;
}

function darken(color, factor) {
  const c = new THREE.Color(color).multiplyScalar(factor);
  return c.getHex();
}

function meshBox(x, y, z, mat) {
  const mesh = new THREE.Mesh(new THREE.BoxGeometry(x, y, z), mat);
  mesh.castShadow = true;
  mesh.receiveShadow = true;
  return mesh;
}

function pickFatFinger(clientX, clientY, radius) {
  const offsets = [
    [0, 0],
    [radius, 0], [-radius, 0], [0, radius], [0, -radius],
    [radius * 0.7, radius * 0.7], [-radius * 0.7, radius * 0.7],
    [radius * 0.7, -radius * 0.7], [-radius * 0.7, -radius * 0.7]
  ];

  let best = null;
  for (const [ox, oy] of offsets) {
    const ndc = clientToNdc(clientX + ox, clientY + oy);
    raycaster.setFromCamera(ndc, camera);
    const hits = raycaster.intersectObjects(dynamicRoots, true);
    if (!hits.length) continue;

    const hit = hits[0];
    const entity = physics.entityFromObject(hit.object);
    if (!entity?.dynamic) continue;

    const radial = Math.hypot(ox, oy);
    const score = hit.distance + radial * 0.018;
    if (!best || score < best.score) {
      best = { entity, point: hit.point.clone(), object: hit.object, score };
    }
  }

  telemetry.sample("input.pick_samples", offsets.length);
  return best;
}

function clientToNdc(x, y) {
  const r = canvas.getBoundingClientRect();
  return {
    x: ((x - r.left) / r.width) * 2 - 1,
    y: -(((y - r.top) / r.height) * 2 - 1)
  };
}

function openContext(entity, clientX, clientY) {
  const r = canvas.getBoundingClientRect();
  contextEl.style.left = `${Math.min(r.width - 170, Math.max(8, clientX - r.left + 12))}px`;
  contextEl.style.top = `${Math.min(r.height - 150, Math.max(8, clientY - r.top - 20))}px`;
  contextEl.hidden = false;
  contextEl.innerHTML = "";

  const make = (label, fn) => {
    const b = document.createElement("button");
    b.textContent = label;
    b.addEventListener("pointerup", e => {
      e.stopPropagation();
      contextEl.hidden = true;
      fn();
      scheduleSave();
    });
    contextEl.appendChild(b);
  };

  make(entity.frozen ? "UNFREEZE" : "FREEZE", () => {
    physics.setFrozen(entity, !entity.frozen);
    haptic(14, 100);
  });
  make("DUPLICATE", () => {
    const s = snapshotSpec(entity);
    s.position[0] += 0.28;
    s.position[2] += 0.18;
    spawnFromSpec(s);
    haptic(14, 90);
  });
  make("DELETE", () => {
    dynamicRoots.splice(dynamicRoots.indexOf(entity.mesh), 1);
    spawnRecords.delete(entity.id);
    physics.removeEntity(entity);
    if (selected === entity) selected = null;
    haptic(18, 125);
  });
}

document.addEventListener("pointerdown", e => {
  if (!contextEl.hidden && !contextEl.contains(e.target)) contextEl.hidden = true;
}, true);

function snapshotSpec(entity) {
  const t = entity.body.translation();
  const q = entity.body.rotation();
  const base = structuredCloneSafe(entity.spec || {});
  base.position = [t.x, t.y, t.z];
  base.rotation = [q.x, q.y, q.z, q.w];
  return base;
}

function captureState() {
  const objects = [];
  for (const entity of physics.entities.values()) {
    if (!entity.dynamic) continue;
    objects.push(snapshotSpec(entity));
  }
  return {
    schema: SAVE_SCHEMA,
    version: "0.8.2-pixel-physics-lab",
    savedAt: Date.now(),
    camera: rig.snapshot(),
    objects
  };
}

function scheduleSave() {
  clearTimeout(saveTimer);
  saveTimer = setTimeout(() => {
    const json = JSON.stringify(captureState());
    globalThis.AndroidBridge?.saveState?.("autosave", json);
  }, 350);
}

function invariantSweep() {
  let nonFinite = 0;
  let maxSpeed = 0;
  for (const entity of physics.entities.values()) {
    if (!entity.dynamic) continue;
    const v = entity.body.linvel();
    const s = Math.hypot(v.x, v.y, v.z);
    maxSpeed = Math.max(maxSpeed, s);
    const t = entity.body.translation();
    if (!Number.isFinite(t.x + t.y + t.z + s)) nonFinite++;
  }

  telemetry.gauge("physics.max_speed", maxSpeed);
  if (nonFinite) telemetry.fault("invariant_non_finite_body", { count: nonFinite });
  if (maxSpeed > 120) telemetry.fault("invariant_runaway_speed", { maxSpeed });
}

function updateStatus() {
  const selectedText = selected
    ? ` · #${selected.id} ${selected.material}`
    : "";
  statusEl.textContent =
    `${fps} FPS · 120Hz physics · 16 solver · CCD×4${selectedText}`;
}

function haptic(ms, amp) {
  globalThis.AndroidBridge?.vibrate?.(ms, amp);
}

window.__onAppPause = () => {
  running = false;
  input.reset();
  grabber.release();
  globalThis.AndroidBridge?.saveState?.("autosave", JSON.stringify(captureState()));
};

window.__onAppResume = () => {
  running = true;
  lastFrame = performance.now();
};

function structuredCloneSafe(value) {
  if (globalThis.structuredClone) return structuredClone(value);
  return JSON.parse(JSON.stringify(value));
}

}

bootstrap().catch(bootFatal);
