# Enterprise v8 Review Board

This document records the multidisciplinary review used to define Enterprise v8. The "team" is implemented as explicit expert roles, each with falsifiable acceptance criteria.

## Reference baseline observed in the supplied APK

The reference build uses Three.js + Rapier 3D. Relevant measured implementation details:

- physics timestep: 1/60 s;
- solver iterations: 8;
- max CCD substeps: 2;
- CCD enabled selectively for small bodies / balls;
- direct pointer capture with one- and two-finger gesture modes;
- long press around 450 ms;
- drag thresholds around 8–10 CSS px;
- grab control applies a COM impulse servo and separately drives angular velocity;
- release adds a bounded correction toward measured finger velocity;
- low-resolution rendering plus pixel/palette/outline post-processing;
- undo/redo, replay, autosave and product UI already exist.

Enterprise v8 must beat this baseline measurably, not rhetorically.

---

## 1. Rigid-body / solver reviewer

### Required changes
- 120 Hz fixed physics step.
- 16 solver iterations.
- 4 CCD substeps.
- 4 internal PGS iterations where supported.
- additional solver iterations on actively manipulated bodies.
- material density drives mass; mass is never a cosmetic constant.
- no position teleport for normal grab interaction.

### v8 implementation
- Rapier integration parameters are configured at startup.
- held bodies receive 8 additional solver iterations.
- adaptive CCD considers predicted travel relative to the smallest body extent.
- small bodies remain CCD-enabled even at lower speed.

### Acceptance
- no body tunnels through a 5 cm obstacle at the defined high-speed test envelope;
- a 10-body stack remains bounded without increasing kinetic energy while untouched;
- physics backlog is observable and dropped catch-up time is counted.

---

## 2. Numerical stability reviewer

### Required changes
- fixed-step accumulator independent of render FPS;
- capped catch-up work;
- transform finite-value guard;
- explicit runaway-velocity fault;
- interpolation from previous to current simulation state.

### v8 implementation
- simulation runs at 1/120 s;
- at most 10 catch-up steps per render frame;
- excess time is dropped and telemetered rather than creating a death spiral;
- visual transforms interpolate between the final two physics states;
- non-finite body state is logged through AndroidBridge.

### Acceptance
- 30 s runtime soak without process death;
- no NaN/Infinity transform;
- p95 physics step time is exported.

---

## 3. Collision / CCD reviewer

### Required changes
- CCD not tied only to object type;
- collision precision follows velocity and geometric scale;
- no render-only approximation controls collision.

### v8 implementation
- adaptive CCD uses speed × dt versus minimum body extent;
- soft CCD prediction is configured when available;
- Rapier colliders use the same physical dimensions as the visible meshes.

### Acceptance
- fast small objects remain inside room bounds;
- collision geometry and visible dimensions differ by less than the asset tolerance.

---

## 4. Touch / motor-control UX reviewer

### Reference weakness
A single pick ray is precise for a mouse but unnecessarily unforgiving for a fingertip.

### Required changes
- fat-finger target acquisition;
- contact-width-aware hit radius;
- pointer capture;
- coalesced pointer samples;
- explicit tap/drag/long-press hysteresis;
- gradual visual finger offset;
- two-finger transition with no target discontinuity.

### v8 implementation
- nine-ray acquisition pattern over a minimum 13 CSS-pixel radius;
- radius expands with PointerEvent width/height;
- coalesced pointer samples drive target updates;
- drag slop: 7 px;
- long-press slop: 11 px;
- long press: 500 ms;
- 10 px assist lift ramps over 90 ms;
- midpoint motion continues to steer the object while pinching/twisting.

### Acceptance
- a visible object whose projected body intersects the acquisition circle is selectable without pixel-perfect aim;
- transition one finger → two fingers → one finger does not release the body;
- context gesture cannot accidentally trigger after movement exceeds long-press slop.

---

## 5. Manipulation-physics reviewer

### Reference weakness
The reference grab servo applies translational impulse at COM, then controls rotation through a second angular-velocity servo. Release also injects a bounded throw correction.

### v8 improvement
The grab servo acts at the exact hit point with applyImpulseAtPoint.

This means off-centre grabs create torque through the rigid-body equations themselves.

The controller is defined from:
- exact local grab point;
- world-space target plane;
- filtered target velocity;
- near-critical PD acceleration;
- mass-aware acceleration/force caps.

Two-finger twist adds bounded torque impulse. It does not overwrite orientation.

Release invariant:

> no artificial throw velocity is injected on release.

The body's release momentum is exactly the momentum produced while it was physically coupled to the moving finger target.

### Acceptance
- grabbing a beam near one end produces more angular response than grabbing its center;
- release with a stationary finger does not create extra velocity;
- heavy objects visibly lag light objects under the same grip limit.

---

## 6. Pixel / rendering reviewer

### Required changes
- real 3D simulation remains visually constrained by a pixel pipeline;
- no antialiasing;
- fixed logical render resolution;
- nearest-neighbour presentation;
- depth-derived outline;
- palette quantization.

### v8 implementation
- 480×270 render target;
- WebGL antialias disabled;
- nearest filtering;
- palette snap in post shader;
- depth discontinuity outline;
- orthographic camera defaults to 45° elevation.

### Acceptance
- one logical pixel is never blurred by presentation scaling;
- render resolution is independent of device resolution.

---

## 7. Android / platform reviewer

### Required changes
- no remote runtime dependencies;
- no cleartext traffic;
- WebView content served from the Android asset origin;
- external navigation blocked;
- file/content access disabled;
- JS/Android bridge bounded and versioned.

### v8 implementation
- WebViewAssetLoader serves app assets through appassets.androidplatform.net;
- allowFileAccess=false;
- allowContentAccess=false;
- mixed content disabled;
- save payload capped at 2 MB and JSON-validated before persistence;
- bridge exposes only haptics, versioned state storage, build info and bounded diagnostics.

### Acceptance
- app starts offline;
- no runtime CDN fetch is needed;
- malformed saves cannot replace valid state.

---

## 8. QA / enterprise operations reviewer

### Required gates
1. JS contract tests.
2. Web bundle build.
3. Android unit tests.
4. signed APK verification.
5. Android 16 cold-start smoke.
6. 30 s interaction soak.
7. screenshot artifact.
8. fault log on any process death.

### Observability
- rolling p50/p95 frame and physics-step metrics;
- catch-up drop counter;
- body count;
- max speed;
- non-finite transform faults;
- JS window errors and unhandled promise rejections.

### Release rule

A build is not "verified" merely because Gradle succeeded.

Release status must distinguish:

- BUILD PASS
- CONTRACT TEST PASS
- RUNTIME PASS
- PERFORMANCE PASS
- DEVICE PLAYTEST UNKNOWN/PASS

---

# Enterprise target

Enterprise v8 is accepted only if it is simultaneously:

1. more physically faithful than the supplied reference baseline;
2. easier to manipulate with one or two fingers;
3. diagnosable after a failure;
4. reproducible in CI;
5. offline-capable and self-contained;
6. explicit about unverified device behavior.
