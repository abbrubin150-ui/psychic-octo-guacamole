# Pixel Physics Sandbox — Isometric 2.5D Pixel Edition

A physics sandbox rendered as **true 2D pixel art in an isometric 2.5D room**.

## Rendering model

This build deliberately does **not** use a 3D renderer, meshes, a perspective camera, or Bullet.

The simulation is a 2D Box2D floor plane plus an independent height channel:

`(x, y)` = floor physics  
`z` = vertical height, gravity and bounce  
`screen = isometric(x, y) + vertical(z)`

The game renders to a fixed **480×320** pixel canvas and scales with nearest-neighbor sampling.

## Art model

Every spawnable prop has four separately authored isometric pixel frames:

- NE
- SE
- SW
- NW

The renderer selects the nearest authored frame from the physical angle. It does not rotate a bitmap.

The workshop itself is drawn as an isometric pixel-art room: diamond floor, two rear walls, timber frame, shelves, window, hanging spring/weight, rails and fixtures.

## Interaction

- one finger: grab and move on the isometric floor plane
- two-finger pinch while grabbed: raise/lower the prop on the height channel
- two-finger twist: physical torque
- release: height gravity and natural Box2D momentum continue
- long press: Freeze / Delete / Duplicate / Material
- Spawn and Undo
- autosave

## Physics

- Box2D handles floor-plane rigid-body collision and momentum
- custom vertical gravity handles the 2.5D height channel
- material-dependent vertical bounce
- spring-damper direct-touch manipulation
- fixed 60 Hz simulation

## Android

Package: `com.pixelphysics.sandbox`  
Version: `0.3.0-isometric-25d`  
Minimum Android: API 26  
Target SDK: API 35
