# Enterprise v8 Architecture

## Runtime layers

1. **Android host**
   - secure WebView
   - local asset origin
   - haptics
   - bounded versioned save storage
   - diagnostic bridge

2. **Input fusion**
   - Pointer Events
   - contact geometry
   - coalesced samples
   - intent state machine
   - fat-finger acquisition

3. **Manipulation controller**
   - exact local hit point
   - camera-parallel target plane
   - near-critical contact-point servo
   - physical torque from off-centre impulse
   - torque impulse for explicit two-finger twist

4. **Physics kernel**
   - Rapier 3D
   - 120 Hz
   - 16 solver iterations
   - adaptive CCD
   - interpolation
   - invariant monitoring

5. **Pixel renderer**
   - Three.js scene
   - orthographic camera
   - 480×270 render target
   - palette quantization
   - depth outlines
   - nearest-neighbour presentation

6. **Product state**
   - schema-versioned autosave
   - stable spawn specifications
   - deterministic restore boundary

7. **Telemetry**
   - local rolling metrics
   - bounded Android logging
   - no user conversation/content telemetry

## Core manipulation invariant

The visible object point under the finger is the physical point receiving the grab impulse.

There is no hidden COM-only translational controller and no release velocity injection.

## Failure policy

- invalid save: reject and preserve previous state;
- non-finite body transform: emit fault, never silently persist it;
- simulation backlog: drop excess wall-clock catch-up and count it;
- app pause: cancel gesture, release body, save stable state;
- app resume: reset frame clock before continuing.

## Next enterprise milestones

- command-sourced undo/redo;
- deterministic input replay;
- golden physics scenes with numeric tolerances;
- device matrix: Samsung/Pixel/tablet, 60/90/120 Hz displays;
- benchmark gates for 25/50/100 dynamic bodies;
- accessibility profile for one-handed depth/rotation;
- save migrations;
- crash reporting connector suitable for production deployment.
