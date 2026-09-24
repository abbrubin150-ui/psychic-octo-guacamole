# Physics v7 — multidisciplinary review

This branch treats the simulation as a constrained 2.5D rigid-body system: continuous XYZ translation, yaw rotation, and cardinal voxel rendering.

## 1. Rigid-body dynamics review

Changes:
- fixed simulation rate raised from 60 Hz to 120 Hz;
- semi-implicit integration;
- mass-aware normal impulses;
- yaw moment of inertia for boxes/cylinders;
- impulse transfer at an actual XY contact point, so off-centre hits generate spin;
- restitution threshold prevents low-speed micro-bounces;
- material friction is handled as a Coulomb impulse rather than frame damping.

Invariant:
> A collision may change linear and angular momentum only through a contact impulse or a declared world constraint.

## 2. Collision-detection review

Old model:
- one approximate radius per object.

v7:
- oriented bounding boxes (OBB) for beams, planks, cubes, crates, stairs, wheels, weights and springs;
- circles for balls and barrel footprints;
- SAT for OBB/OBB;
- exact closest-point circle/OBB manifold;
- circle/circle manifold;
- broadphase AABB rejection before narrowphase;
- Z interval separation decides side contact vs support contact.

The visual cardinal model and collision footprint share the same metric dimensions: one voxel is 0.01 m.

## 3. Numerical-stability review

- maximum travel per substep: 0.018 m;
- up to four conservative substeps;
- nine sequential velocity iterations;
- four positional iterations;
- contact slop: 1.2 mm;
- Baumgarte-style positional correction: 0.68;
- sleeping only when grounded and linear/vertical/angular speeds are all below thresholds;
- sleeping bodies wake on impulse or grab.

These controls target the common sandbox failures: tunnelling, stack jitter, energy injection, interpenetration and endlessly vibrating resting bodies.

## 4. Interaction-physics review

The grab manipulator is now a near-critically damped mass-aware spring:
- stiffness scales sublinearly with mass;
- damping is derived from sqrt(k*m);
- maximum force is capped;
- no artificial release velocity is injected.

This preserves the "touch = physical manipulator" design while reducing oscillation and slingshot energy.

## Stairs

Stairs retain an OBB side collider but their vertical support height is sampled from the actual eight-step voxel profile. Objects therefore rest on the current step height rather than on the staircase bounding-box roof.

## Deterministic geometry tests

The CI suite checks:
- separated OBBs;
- metric penetration depth;
- rotated OBB SAT;
- circle/circle depth;
- circle/OBB normal direction and symmetry;
- rotated AABB expansion;
- point-in-rotated-box behavior.

## Known model boundary

v7 is still a 2.5D solver, not a six-degree-of-freedom 3D rigid-body engine. It models XYZ translation and yaw, not pitch/roll. Irregular voxel assets use convex XY footprints except for the staircase's explicit stepped support surface. That limitation is intentional so the simulation stays visually consistent with the four-cardinal voxel renderer.
