import * as THREE from "three";
import RAPIER from "@dimforge/rapier3d-compat";
import { PHYSICS } from "./config.js";

const tmpV1 = new THREE.Vector3();
const tmpV2 = new THREE.Vector3();
const tmpV3 = new THREE.Vector3();
const tmpQ1 = new THREE.Quaternion();
const tmpQ2 = new THREE.Quaternion();

export class PhysicsKernel {
  constructor(telemetry) {
    this.telemetry = telemetry;
    this.world = null;
    this.eventQueue = null;
    this.entities = new Map();
    this.nextId = 1;
    this.accumulator = 0;
    this.droppedTime = 0;
  }

  async init() {
    await RAPIER.init();

    this.world = new RAPIER.World(PHYSICS.gravity);
    this.eventQueue = new RAPIER.EventQueue(true);

    const ip = this.world.integrationParameters;
    ip.dt = PHYSICS.fixedDt;
    ip.numSolverIterations = PHYSICS.solverIterations;
    ip.maxCcdSubsteps = PHYSICS.maxCcdSubsteps;

    if ("numInternalPgsIterations" in ip) {
      ip.numInternalPgsIterations = PHYSICS.internalPgsIterations;
    }
    if ("contact_natural_frequency" in ip) {
      ip.contact_natural_frequency = 45;
    }

    return this;
  }

  createStaticFloor(mesh, halfExtents = [4.8, 0.10, 4.8], center = [0, -0.10, 0]) {
    const body = this.world.createRigidBody(
      RAPIER.RigidBodyDesc.fixed().setTranslation(center[0], center[1], center[2])
    );
    const collider = this.world.createCollider(
      RAPIER.ColliderDesc.cuboid(halfExtents[0], halfExtents[1], halfExtents[2])
        .setFriction(0.82)
        .setRestitution(0.02),
      body
    );

    const entity = this.#registerEntity({
      body,
      collider,
      mesh,
      dynamic: false,
      minExtent: Math.min(...halfExtents) * 2,
      material: "floor"
    });
    return entity;
  }

  spawnBox({ mesh, position, halfExtents, material = "wood", rotation = null }) {
    const density = materialDensity(material);
    let desc = RAPIER.RigidBodyDesc.dynamic()
      .setTranslation(position.x, position.y, position.z)
      .setLinearDamping(0.015)
      .setAngularDamping(0.035)
      .setCanSleep(true);

    if (rotation) {
      desc = desc.setRotation({ x: rotation.x, y: rotation.y, z: rotation.z, w: rotation.w });
    }

    const minExtent = Math.min(halfExtents.x, halfExtents.y, halfExtents.z) * 2;
    if (minExtent <= 0.45) desc = desc.setCcdEnabled(true);

    const body = this.world.createRigidBody(desc);
    if (body.setAdditionalSolverIterations) {
      body.setAdditionalSolverIterations(PHYSICS.baseAdditionalSolverIterations);
    }

    const colliderDesc = RAPIER.ColliderDesc.cuboid(halfExtents.x, halfExtents.y, halfExtents.z)
      .setDensity(density)
      .setFriction(materialFriction(material))
      .setRestitution(materialRestitution(material));

    if (RAPIER.CoefficientCombineRule) {
      colliderDesc.setFrictionCombineRule(RAPIER.CoefficientCombineRule.Average);
      colliderDesc.setRestitutionCombineRule(RAPIER.CoefficientCombineRule.Max);
    }

    const collider = this.world.createCollider(colliderDesc, body);
    return this.#registerEntity({ body, collider, mesh, dynamic: true, minExtent, material });
  }

  spawnBall({ mesh, position, radius, material = "rubber" }) {
    let desc = RAPIER.RigidBodyDesc.dynamic()
      .setTranslation(position.x, position.y, position.z)
      .setLinearDamping(0.01)
      .setAngularDamping(0.025)
      .setCanSleep(true)
      .setCcdEnabled(true);

    const body = this.world.createRigidBody(desc);
    if (body.setAdditionalSolverIterations) {
      body.setAdditionalSolverIterations(PHYSICS.baseAdditionalSolverIterations);
    }

    const colliderDesc = RAPIER.ColliderDesc.ball(radius)
      .setDensity(materialDensity(material))
      .setFriction(materialFriction(material))
      .setRestitution(materialRestitution(material));

    if (RAPIER.CoefficientCombineRule) {
      colliderDesc.setFrictionCombineRule(RAPIER.CoefficientCombineRule.Average);
      colliderDesc.setRestitutionCombineRule(RAPIER.CoefficientCombineRule.Max);
    }

    const collider = this.world.createCollider(colliderDesc, body);
    return this.#registerEntity({
      body,
      collider,
      mesh,
      dynamic: true,
      minExtent: radius * 2,
      material
    });
  }

  #registerEntity({ body, collider, mesh, dynamic, minExtent, material }) {
    const id = this.nextId++;
    const t = body.translation();
    const r = body.rotation();
    const entity = {
      id,
      body,
      collider,
      mesh,
      dynamic,
      minExtent,
      material,
      prevPos: new THREE.Vector3(t.x, t.y, t.z),
      currPos: new THREE.Vector3(t.x, t.y, t.z),
      prevQuat: new THREE.Quaternion(r.x, r.y, r.z, r.w),
      currQuat: new THREE.Quaternion(r.x, r.y, r.z, r.w)
    };

    mesh.userData.entityId = id;
    mesh.traverse?.(child => {
      child.userData.entityId = id;
    });

    this.entities.set(id, entity);
    return entity;
  }

  entityFromObject(object) {
    let o = object;
    while (o) {
      const id = o.userData?.entityId;
      if (id && this.entities.has(id)) return this.entities.get(id);
      o = o.parent;
    }
    return null;
  }

  stepFrame(realDt, beforeStep) {
    if (!this.world) return 0;

    const frameDt = Math.min(Math.max(realDt, 0), PHYSICS.maxFrameDt);
    this.accumulator += frameDt;

    let steps = 0;
    const start = performance.now();

    while (this.accumulator >= PHYSICS.fixedDt && steps < PHYSICS.maxCatchupSteps) {
      this.#capturePrevious();
      beforeStep?.(PHYSICS.fixedDt);
      this.#adaptiveCcd();
      this.world.timestep = PHYSICS.fixedDt;
      this.world.step(this.eventQueue);
      this.#captureCurrent();

      this.accumulator -= PHYSICS.fixedDt;
      steps++;
    }

    if (this.accumulator >= PHYSICS.fixedDt) {
      this.droppedTime += this.accumulator;
      this.telemetry?.count("physics.catchup_drop");
      this.accumulator = this.accumulator % PHYSICS.fixedDt;
    }

    this.telemetry?.sample("physics.step_ms", performance.now() - start);
    this.telemetry?.sample("physics.steps_per_frame", steps);
    this.telemetry?.gauge("physics.body_count", this.entities.size);

    return Math.min(1, this.accumulator / PHYSICS.fixedDt);
  }

  syncVisuals(alpha) {
    for (const e of this.entities.values()) {
      if (!e.mesh || !e.dynamic) continue;

      e.mesh.position.lerpVectors(e.prevPos, e.currPos, alpha);
      e.mesh.quaternion.slerpQuaternions(e.prevQuat, e.currQuat, alpha);
    }
  }

  #capturePrevious() {
    for (const e of this.entities.values()) {
      if (!e.dynamic) continue;
      e.prevPos.copy(e.currPos);
      e.prevQuat.copy(e.currQuat);
    }
  }

  #captureCurrent() {
    for (const e of this.entities.values()) {
      if (!e.dynamic) continue;
      const t = e.body.translation();
      const q = e.body.rotation();

      if (!Number.isFinite(t.x + t.y + t.z + q.x + q.y + q.z + q.w)) {
        this.telemetry?.count("physics.non_finite_transform");
        this.telemetry?.fault("non_finite_transform", { id: e.id });
        continue;
      }

      e.currPos.set(t.x, t.y, t.z);
      e.currQuat.set(q.x, q.y, q.z, q.w);
    }
  }

  #adaptiveCcd() {
    for (const e of this.entities.values()) {
      if (!e.dynamic) continue;
      const v = e.body.linvel();
      const speed = Math.hypot(v.x, v.y, v.z);
      const travel = speed * PHYSICS.fixedDt;
      const need = travel > e.minExtent * PHYSICS.ccdTravelFraction || e.minExtent <= 0.30;

      if (e.body.enableCcd) e.body.enableCcd(need);
      if (need && e.body.setSoftCcdPrediction) {
        e.body.setSoftCcdPrediction(Math.max(0.002, e.minExtent * 0.5));
      }
    }
  }
}

export class ContactPointGrabber {
  constructor(kernel, telemetry) {
    this.kernel = kernel;
    this.telemetry = telemetry;
    this.active = null;
    this.raycaster = new THREE.Raycaster();
    this.plane = new THREE.Plane();
    this.planePoint = new THREE.Vector3();
    this.target = new THREE.Vector3();
    this.targetVelocity = new THREE.Vector3();
    this.lastTarget = new THREE.Vector3();
    this.axis = new THREE.Vector3();
    this.twistRate = 0;
    this.depthOffset = 0;
    this.lastSampleTime = 0;
  }

  start(entity, hitPoint, camera, ndcX, ndcY, timeSec = performance.now() / 1000) {
    this.release();

    const body = entity.body;
    const p = body.translation();
    const q = body.rotation();

    const pos = tmpV1.set(p.x, p.y, p.z);
    const quat = tmpQ1.set(q.x, q.y, q.z, q.w);

    const localGrab = hitPoint.clone().sub(pos).applyQuaternion(tmpQ2.copy(quat).invert());
    const cameraDir = camera.getWorldDirection(new THREE.Vector3()).normalize();

    this.planePoint.copy(hitPoint);
    this.plane.setFromNormalAndCoplanarPoint(cameraDir, this.planePoint);
    this.axis.copy(cameraDir);

    this.active = {
      entity,
      localGrab: localGrab.clone(),
      startedAt: timeSec,
      heldSolverIterations: PHYSICS.heldAdditionalSolverIterations
    };

    this.target.copy(hitPoint);
    this.lastTarget.copy(hitPoint);
    this.targetVelocity.set(0, 0, 0);
    this.lastSampleTime = timeSec;
    this.depthOffset = 0;
    this.twistRate = 0;

    if (body.setAdditionalSolverIterations) {
      body.setAdditionalSolverIterations(PHYSICS.heldAdditionalSolverIterations);
    }
    if (body.enableCcd) body.enableCcd(true);
    if (body.setSoftCcdPrediction) body.setSoftCcdPrediction(Math.max(0.002, entity.minExtent * 0.7));
    body.wakeUp?.();

    this.setScreen(ndcX, ndcY, camera, timeSec);
    return this.active;
  }

  setScreen(ndcX, ndcY, camera, timeSec = performance.now() / 1000) {
    if (!this.active) return;

    this.raycaster.setFromCamera({ x: ndcX, y: ndcY }, camera);
    const movedPlanePoint = tmpV1.copy(this.planePoint).addScaledVector(this.axis, this.depthOffset);
    this.plane.setFromNormalAndCoplanarPoint(this.axis, movedPlanePoint);

    const hit = this.raycaster.ray.intersectPlane(this.plane, tmpV2);
    if (!hit) return;

    const dt = Math.max(1 / 240, timeSec - this.lastSampleTime);
    const rawVelocity = tmpV3.copy(hit).sub(this.lastTarget).multiplyScalar(1 / dt);

    this.targetVelocity.lerp(rawVelocity, 0.32);
    this.lastTarget.copy(hit);
    this.target.copy(hit);
    this.lastSampleTime = timeSec;
  }

  adjustDepth(deltaCssPx, camera, ndcX, ndcY, timeSec) {
    if (!this.active) return;
    this.depthOffset = THREE.MathUtils.clamp(this.depthOffset - deltaCssPx * 0.006, -4, 4);
    this.setScreen(ndcX, ndcY, camera, timeSec);
  }

  addTwist(deltaRadians, dt) {
    if (!this.active) return;
    const rate = deltaRadians / Math.max(dt, 1 / 240);
    this.twistRate = THREE.MathUtils.lerp(this.twistRate, rate, 0.45);
  }

  preStep(dt) {
    const a = this.active;
    if (!a) return;

    const body = a.entity.body;
    const t = body.translation();
    const r = body.rotation();
    const lin = body.linvel();
    const ang = body.angvel();

    const pos = tmpV1.set(t.x, t.y, t.z);
    const quat = tmpQ1.set(r.x, r.y, r.z, r.w);
    const worldGrab = tmpV2.copy(a.localGrab).applyQuaternion(quat).add(pos);
    const arm = tmpV3.copy(worldGrab).sub(pos);

    const pointVelocity = new THREE.Vector3(lin.x, lin.y, lin.z)
      .add(new THREE.Vector3(ang.x, ang.y, ang.z).cross(arm));

    const error = this.target.clone().sub(worldGrab);
    const omega = 2 * Math.PI * PHYSICS.grabFrequencyHz;
    const accel = error.multiplyScalar(omega * omega)
      .add(
        this.targetVelocity.clone()
          .sub(pointVelocity)
          .multiplyScalar(2 * PHYSICS.grabDampingRatio * omega)
      );

    const mass = Math.max(0.05, body.mass());
    const maxAccel = Math.max(24, PHYSICS.maxGrabAcceleration / (1 + 0.14 * Math.sqrt(mass)));
    if (accel.length() > maxAccel) accel.setLength(maxAccel);

    const impulse = accel.multiplyScalar(mass * dt);
    if (Number.isFinite(impulse.x + impulse.y + impulse.z)) {
      body.applyImpulseAtPoint(
        { x: impulse.x, y: impulse.y, z: impulse.z },
        { x: worldGrab.x, y: worldGrab.y, z: worldGrab.z },
        true
      );
    }

    if (Math.abs(this.twistRate) > 0.01 && body.applyTorqueImpulse) {
      const av = new THREE.Vector3(ang.x, ang.y, ang.z);
      const along = av.dot(this.axis);
      const desired = THREE.MathUtils.clamp(this.twistRate, -7, 7);
      const delta = desired - along;

      const magnitude = THREE.MathUtils.clamp(
        delta * Math.max(0.02, mass * 0.035),
        -PHYSICS.maxAngularImpulse,
        PHYSICS.maxAngularImpulse
      );

      body.applyTorqueImpulse(
        {
          x: this.axis.x * magnitude,
          y: this.axis.y * magnitude,
          z: this.axis.z * magnitude
        },
        true
      );
    }

    this.twistRate *= 0.72;
  }

  release() {
    if (!this.active) return;
    const body = this.active.entity.body;

    // Enterprise invariant: release never injects an artificial throw velocity.
    // The body's momentum is whatever the contact-point servo physically produced.
    if (body.setAdditionalSolverIterations) {
      body.setAdditionalSolverIterations(PHYSICS.baseAdditionalSolverIterations);
    }

    this.telemetry?.sample(
      "input.grab_duration_ms",
      (performance.now() / 1000 - this.active.startedAt) * 1000
    );

    this.active = null;
    this.twistRate = 0;
    this.targetVelocity.set(0, 0, 0);
  }
}

function materialDensity(material) {
  switch (material) {
    case "metal": return 7800;
    case "rubber": return 1100;
    case "wood": return 700;
    default: return 1000;
  }
}

function materialFriction(material) {
  switch (material) {
    case "metal": return 0.42;
    case "rubber": return 0.92;
    case "wood": return 0.68;
    default: return 0.65;
  }
}

function materialRestitution(material) {
  switch (material) {
    case "metal": return 0.08;
    case "rubber": return 0.72;
    case "wood": return 0.12;
    default: return 0.1;
  }
}

export { RAPIER };
