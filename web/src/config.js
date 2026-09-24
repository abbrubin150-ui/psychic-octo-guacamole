export const PHYSICS = Object.freeze({
  fixedDt: 1 / 120,
  maxFrameDt: 0.08,
  maxCatchupSteps: 10,
  solverIterations: 16,
  internalPgsIterations: 4,
  maxCcdSubsteps: 4,
  baseAdditionalSolverIterations: 2,
  heldAdditionalSolverIterations: 8,
  gravity: Object.freeze({ x: 0, y: -9.81, z: 0 }),
  ccdTravelFraction: 0.35,
  maxGrabAcceleration: 95,
  grabFrequencyHz: 9.5,
  grabDampingRatio: 0.96,
  maxAngularImpulse: 1.8,
  releaseVelocityInjection: 0
});

export const INPUT = Object.freeze({
  dragSlopCssPx: 7,
  longPressSlopCssPx: 11,
  longPressMs: 500,
  doubleTapMs: 320,
  doubleTapRadiusCssPx: 28,
  fatFingerRadiusCssPx: 13,
  fatFingerSamples: 9,
  assistLiftCssPx: 10,
  assistRampMs: 90,
  sampleWindowMs: 90,
  velocitySmoothing: 0.32
});

export const RENDER = Object.freeze({
  width: 480,
  height: 270,
  paletteLevels: 8,
  outlineDepthThreshold: 0.0022,
  cameraExtent: 7.2
});

export const SAVE_SCHEMA = 1;
