# Physics v7 — Expert Review & Validation Protocol

This document records the multidisciplinary review used for the v7 physics rewrite.

## 1. Rigid-body dynamics review

### Problem in v6
The old solver reduced every horizontal body interaction to a circular radius. Long beams, crates, stairs and thin planks therefore collided with geometry that did not match their visible footprint.

### v7 change
- Rectangular voxel props use oriented bounding boxes (OBB).
- Balls and barrels use circular footprints.
- OBB-vs-OBB uses a 4-axis separating axis test (SAT).
- Circle-vs-circle and circle-vs-OBB have dedicated narrow-phase solvers.
- Contact points are retained for impulse application.
- Yaw angular inertia is derived from body footprint and mass.

## 2. Contact solver review

v7 uses sequential impulses instead of direct velocity reflection.

Contact constraints contain:
- normal
- contact point
- penetration
- friction
- restitution
- body A / body B

Solver:
- 10 velocity iterations
- 4 position iterations
- Baumgarte-style positional correction with slop
- Coulomb friction clamped by the normal impulse
- restitution disabled below a velocity threshold to suppress micro-bounce

## 3. Vertical / stacking review

The system is 2.5D:
- X/Y = planar rigid-body dynamics
- Z = vertical position and velocity
- yaw = rotational degree of freedom

For body-body contacts, horizontal penetration and Z penetration are compared. The minimum separating axis decides whether the contact is a side collision or a stacking contact.

Stairs use a local height-field with 8 discrete treads rather than one flat top surface.

## 4. Numerical integration review

- fixed outer timestep: 1/60 s
- semi-implicit velocity integration
- adaptive collision substeps: 1–6
- substep count is driven by travel distance relative to the smallest active feature
- weak air drag only; contact friction is responsible for surface damping

This is intended to reduce fast-body tunneling without requiring a native CCD library.

## 5. Manipulation review

Touch grabbing is now applied at the actual local grab offset:
- point velocity includes yaw angular velocity
- spring force acts at that point
- off-center force therefore creates torque naturally
- two-finger twist drives an angular spring target rather than directly overwriting the body angle

This removes the previous non-physical direct rotation injection.

## 6. Sleeping review

Dynamic bodies can sleep only when:
- horizontal speed < 0.035 m/s
- vertical speed < 0.045 m/s
- yaw speed < 0.08 rad/s
- the body has continuous support
- the condition persists for 0.70 s

A sufficiently strong collision or active grab wakes the body.

## 7. Material model

Working coefficients:

| Material | restitution | friction |
|---|---:|---:|
| Mahogany | 0.14 | 0.64 |
| Metal | 0.08 | 0.36 |
| Rubber | 0.72 | 0.96 |

Body-body coefficients use geometric mixing.

## 8. Required regression scenarios

### STACK-10
Drop ten identical cubes into a vertical stack.
Pass target:
- no explosive separation
- no perpetual micro-bounce
- no visible lateral drift after settling
- sleep reached after settling

### BEAM-CORNER
Strike one end of a beam against a crate.
Pass target:
- impulse produces yaw torque
- beam does not behave like a circle
- collision point is visibly consistent with the contacted corner

### RUBBER-DROP
Drop rubber, wood and metal balls from equal height.
Pass target:
- rubber rebound is visibly highest
- metal/wood settle without repeated micro-bounces

### FRICTION-SLIDE
Give equal horizontal velocity to wood, metal and rubber bodies.
Pass target:
- metal travels furthest
- rubber travels least
- no instant velocity deletion

### FAST-BODY
Launch the smallest body at high velocity toward another body and a room wall.
Pass target:
- adaptive substeps prevent ordinary gameplay-speed tunneling
- no NaN / infinity state

### GRAB-OFFCENTER
Grab the end of a long beam and drag laterally.
Pass target:
- translation and rotation appear together
- center of mass is not teleported to the finger
- release preserves physical velocity

### STAIRS
Drop a cube onto different staircase treads.
Pass target:
- cube rests on the local tread height
- cube does not float on one flat invisible lid

### SLEEP-WAKE
Allow a stack to settle, then strike the bottom body.
Pass target:
- settled bodies stop consuming integration movement
- impact wakes affected bodies

## 9. Current model limits

v7 is still a custom 2.5D solver, not a general 6-DOF rigid-body engine.

Not modeled:
- pitch/roll angular dynamics
- arbitrary convex 3D hull collision
- gyroscopic effects
- true rolling constraints for upright wheels
- continuous-time exact TOI collision

These are deliberate boundaries. Within the game’s fixed isometric voxel presentation, the target is stable, coherent 2.5D rigid-body behavior rather than pretending to provide full 3D mechanics.
