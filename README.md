# Pixel Physics Sandbox — Pure 2D Pixel Edition

A direct-touch physics sandbox rebuilt as a **pure 2D pixel-art game**.

## Visual invariant

There is no 3D renderer, no mesh pipeline, no perspective camera, and no Bullet physics.

The game renders to a fixed 480×270 canvas and scales with nearest-neighbor integer scaling.

Every prop is represented by **four separately authored pixel sprites**:

- 0°
- 90°
- 180°
- 270°

The renderer never rotates a sprite image. Box2D may simulate continuous rotation, but the visual representation selects the nearest authored directional frame.

The Spawn Drawer uses the same four-frame sprite sets and cycles through them as a visible verification that each icon has multiple authored directions.

## Physics

- Box2D 2D rigid bodies
- direct touch grab at the actual contact point
- spring-damper force coupling
- partial mass compensation
- natural release/throw
- two-finger torque
- freeze/unfreeze
- delete / duplicate
- discrete undo
- autosave

## Pixel assets

The MVP contains four-direction sprite sets for:

Cube, Beam Short, Beam Long, Plank, Wood Ball, Metal Ball, Weight, Wheel, Rubber Ball, Barrel, Crate, and Ramp.

Each material has a restricted hand-authored palette: mahogany, metal, and rubber.

## Android

Package: `com.pixelphysics.sandbox`

Version: `0.2.0-pure-2d-pixel`

Minimum Android: API 26

Target SDK: API 35
