import test from "node:test";
import assert from "node:assert/strict";
import { PHYSICS, INPUT, RENDER } from "../src/config.js";

test("physics baseline exceeds reference solver cadence", () => {
  assert.equal(PHYSICS.fixedDt, 1 / 120);
  assert.ok(PHYSICS.solverIterations >= 16);
  assert.ok(PHYSICS.maxCcdSubsteps >= 4);
  assert.equal(PHYSICS.releaseVelocityInjection, 0);
});

test("touch thresholds remain finger-scale and internally ordered", () => {
  assert.ok(INPUT.dragSlopCssPx < INPUT.longPressSlopCssPx);
  assert.ok(INPUT.fatFingerRadiusCssPx >= INPUT.longPressSlopCssPx);
  assert.ok(INPUT.longPressMs >= 450 && INPUT.longPressMs <= 650);
  assert.ok(INPUT.sampleWindowMs >= 60 && INPUT.sampleWindowMs <= 120);
});

test("pixel renderer is fixed logical resolution", () => {
  assert.equal(RENDER.width, 480);
  assert.equal(RENDER.height, 270);
  assert.ok(RENDER.paletteLevels >= 6);
});
