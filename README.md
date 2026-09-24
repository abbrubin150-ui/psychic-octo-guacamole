# Pixel Physics Sandbox — Voxel Orthographic v6

This branch rebuilds the visual pipeline around **true voxel art under a fixed orthographic 45° camera**.

## Core raster invariant

- Simulation units are meters.
- **0.01 m = 1 voxel unit = 1 viewport pixel unit.**
- The camera is orthographic at 45°.
- The 45° loss factor is compensated by **sqrt(2)** in the Y/Z projection terms.
- After compensation and integer snap, the renderer uses the exact raster basis:

```
pixelX = X - Y
pixelY = (X + Y) / 2 - Z
```

where X/Y/Z are voxel coordinates.

Every exposed voxel top/front face is emitted as a single viewport pixel. There is no texture filtering and no pseudo-pixel post effect.

## Asset model

Voxel assets are stored as one 3D voxel model, not four hand-painted sprites.

At runtime the same model is rotated into the four cardinal orientations and rasterized through the compensated orthographic projection. This means a staircase, crate, beam, wheel, barrel, spring, etc. remain visually coherent when rotated.

Current voxel models:

- Cube
- Beam Short
- Beam Long
- Plank
- Weight
- Wheel
- Barrel
- Crate
- Staircase
- Spring

The three balls are intentionally **non-voxel hybrid objects**, matching the reference technique where smooth/non-voxel elements can be mixed into the voxel scene.

## Rendering implementation

The renderer is implemented in pure Android Canvas/Java. It uses:

- fixed 480×270 logical canvas
- integer nearest-neighbor presentation scaling
- software orthographic projection
- per-face single-pixel emission
- a tiny per-sprite depth buffer
- deterministic material shading
- cardinal model rotation from one voxel dataset

No Bullet, Box2D, libGDX, OpenGL engine, JNI, or native `.so` libraries are required.

## Physics / interaction

The v5 custom XYZ physics remains:

- gravity on Z
- floor/object support
- stacking
- horizontal collisions
- material friction/restitution
- direct grab via spring-damper coupling
- two-finger height control
- twist rotation, snapped to the nearest cardinal direction on release
- freeze / delete / duplicate / material cycle
- undo
- autosave

## Android

Package: `com.pixelphysics.sandbox`

Version: `0.5.0-voxel-ortho`

Version code: `6`

Minimum Android: API 26

Target SDK: API 35
