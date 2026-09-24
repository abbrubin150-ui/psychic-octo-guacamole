import { INPUT } from "./config.js";

export class PointerFusion {
  constructor(element, handlers, telemetry) {
    this.el = element;
    this.h = handlers;
    this.telemetry = telemetry;
    this.ptrs = new Map();
    this.mode = "none";
    this.longTimer = null;
    this.pressHit = null;
    this.gesture = null;
    this.lastTap = null;

    element.style.touchAction = "none";
    element.addEventListener("pointerdown", e => this.down(e), { passive: false });
    element.addEventListener("pointermove", e => this.move(e), { passive: false });
    element.addEventListener("pointerup", e => this.up(e, false), { passive: false });
    element.addEventListener("pointercancel", e => this.up(e, true), { passive: false });
    element.addEventListener("contextmenu", e => e.preventDefault());
  }

  down(e) {
    e.preventDefault();
    try { this.el.setPointerCapture(e.pointerId); } catch {}

    const p = {
      id: e.pointerId,
      x: e.clientX,
      y: e.clientY,
      sx: e.clientX,
      sy: e.clientY,
      t0: e.timeStamp * 0.001,
      lastT: e.timeStamp * 0.001,
      type: e.pointerType,
      width: e.width || 0,
      height: e.height || 0,
      pressure: e.pressure || 0,
      followOffsetX: 0,
      followOffsetY: 0
    };
    this.ptrs.set(e.pointerId, p);

    if (this.ptrs.size === 1) {
      const radius = Math.max(
        INPUT.fatFingerRadiusCssPx,
        0.5 * Math.max(p.width, p.height)
      );
      this.pressHit = this.h.pick?.(p.x, p.y, radius) || null;

      if (this.pressHit?.entity) {
        this.mode = "objPending";
        this.h.beginGrab?.(this.pressHit, p.x, p.y, p.t0);
        this.#startLongPress();
      } else {
        this.mode = "camOrbit";
      }
      return;
    }

    if (this.ptrs.size === 2) {
      this.#cancelLongPress();
      const [a, b] = [...this.ptrs.values()];
      const mid = midpoint(a, b);
      const dist = distance(a, b);
      const ang = angle(a, b);

      if ((this.mode === "objPending" || this.mode === "objDrag") && this.h.isGrabbing?.()) {
        this.mode = "objGesture";
      } else {
        this.mode = "camGesture";
      }

      const first = [...this.ptrs.values()][0];
      const objectHandoff = this.mode === "objGesture";
      this.gesture = {
        lastDist: dist,
        lastAngle: ang,
        lastMidX: mid.x,
        lastMidY: mid.y,
        lastT: e.timeStamp * 0.001,
        targetOffsetX: objectHandoff ? first.x - mid.x : 0,
        targetOffsetY: objectHandoff ? first.y - mid.y : 0
      };
    }
  }

  move(e) {
    const p = this.ptrs.get(e.pointerId);
    if (!p) return;
    e.preventDefault();

    const samples = e.getCoalescedEvents?.() || [e];
    for (const sample of samples) {
      const oldX = p.x;
      const oldY = p.y;
      const oldT = p.lastT;

      p.x = sample.clientX;
      p.y = sample.clientY;
      p.lastT = sample.timeStamp * 0.001;
      p.pressure = sample.pressure || p.pressure;

      const travel = Math.hypot(p.x - p.sx, p.y - p.sy);
      const sampleDt = Math.max(1 / 240, p.lastT - oldT);

      if (this.mode === "objPending" && travel > INPUT.dragSlopCssPx) {
        this.mode = "objDrag";
        this.#cancelLongPress();
        this.h.onDragStart?.();
      }

      if (this.mode === "objDrag" && this.ptrs.size === 1) {
        const heldMs = Math.max(0, (p.lastT - p.t0) * 1000);
        const assist = INPUT.assistLiftCssPx * Math.min(1, heldMs / INPUT.assistRampMs);

        const decay = Math.exp(-sampleDt / 0.11);
        p.followOffsetX *= decay;
        p.followOffsetY *= decay;

        this.h.moveGrab?.(
          p.x + p.followOffsetX,
          p.y + p.followOffsetY - assist,
          p.lastT,
          p.pressure
        );
      } else if (this.mode === "camOrbit" && this.ptrs.size === 1) {
        const dx = p.x - oldX;
        const dy = p.y - oldY;
        if (travel > 3) this.h.orbit?.(dx, dy, p.lastT);
      }
    }

    if (this.mode === "objGesture" || this.mode === "camGesture") {
      this.#updateGesture(e.timeStamp * 0.001);
    }
  }

  #updateGesture(timeSec) {
    if (this.ptrs.size < 2 || !this.gesture) return;
    const [a, b] = [...this.ptrs.values()];
    const g = this.gesture;
    const mid = midpoint(a, b);
    const dist = Math.max(1, distance(a, b));
    let ang = angle(a, b);
    let da = ang - g.lastAngle;
    if (da > Math.PI) da -= Math.PI * 2;
    if (da < -Math.PI) da += Math.PI * 2;

    const dd = dist - g.lastDist;
    const dt = Math.max(1 / 240, timeSec - g.lastT);

    if (this.mode === "objGesture") {
      const decay = Math.exp(-dt / 0.12);
      g.targetOffsetX *= decay;
      g.targetOffsetY *= decay;

      this.h.gestureGrab?.({
        x: mid.x + g.targetOffsetX,
        y: mid.y + g.targetOffsetY,
        depthDeltaPx: dd,
        twistDelta: da,
        dt,
        timeSec
      });
    } else {
      this.h.cameraGesture?.({
        x: mid.x,
        y: mid.y,
        panX: mid.x - g.lastMidX,
        panY: mid.y - g.lastMidY,
        pinchDeltaPx: dd,
        twistDelta: da,
        dt
      });
    }

    g.lastDist = dist;
    g.lastAngle = ang;
    g.lastMidX = mid.x;
    g.lastMidY = mid.y;
    g.lastT = timeSec;
  }

  up(e, cancelled) {
    const p = this.ptrs.get(e.pointerId);
    if (!p) return;
    e.preventDefault();

    this.ptrs.delete(e.pointerId);
    const durationMs = Math.max(0, (e.timeStamp * 0.001 - p.t0) * 1000);
    const travel = Math.hypot(p.x - p.sx, p.y - p.sy);

    if (this.ptrs.size === 1) {
      const remaining = [...this.ptrs.values()][0];
      const oldGesture = this.gesture;

      if (this.mode === "objGesture" && oldGesture) {
        remaining.followOffsetX = oldGesture.lastMidX - remaining.x;
        remaining.followOffsetY = oldGesture.lastMidY - remaining.y;
      } else {
        remaining.followOffsetX = 0;
        remaining.followOffsetY = 0;
      }

      remaining.sx = remaining.x;
      remaining.sy = remaining.y;
      remaining.t0 = e.timeStamp * 0.001;
      this.gesture = null;

      if (this.mode === "objGesture") this.mode = "objDrag";
      else if (this.mode === "camGesture") this.mode = "camOrbit";
      return;
    }

    if (this.ptrs.size > 1) return;

    this.#cancelLongPress();
    const oldMode = this.mode;
    this.mode = "none";
    this.gesture = null;

    if (cancelled) {
      if (oldMode.startsWith("obj")) this.h.endGrab?.({ cancelled: true });
      return;
    }

    if (oldMode === "objPending") {
      this.h.endGrab?.({ cancelled: false, tapLike: true });
      this.#handleTap(p, this.pressHit, durationMs, travel);
    } else if (oldMode === "objDrag" || oldMode === "objGesture") {
      this.h.endGrab?.({ cancelled: false, tapLike: false });
      this.telemetry?.sample("input.drag_distance_px", travel);
    } else if (oldMode === "camOrbit" && durationMs < 400 && travel < 9) {
      const radius = Math.max(INPUT.fatFingerRadiusCssPx, 0.5 * Math.max(p.width, p.height));
      const hit = this.h.pick?.(p.x, p.y, radius) || null;
      this.#handleTap(p, hit, durationMs, travel);
    }
  }

  #handleTap(p, hit) {
    const now = performance.now();
    const entId = hit?.entity?.id ?? null;
    const prev = this.lastTap;

    if (
      prev &&
      now - prev.t < INPUT.doubleTapMs &&
      Math.hypot(p.x - prev.x, p.y - prev.y) < INPUT.doubleTapRadiusCssPx &&
      prev.entId === entId
    ) {
      this.lastTap = null;
      this.h.doubleTap?.(hit, p.x, p.y);
      return;
    }

    this.lastTap = { t: now, x: p.x, y: p.y, entId };
    this.h.tap?.(hit, p.x, p.y);
  }

  #startLongPress() {
    this.#cancelLongPress();
    this.longTimer = setTimeout(() => {
      this.longTimer = null;
      if (this.mode !== "objPending") return;
      const p = [...this.ptrs.values()][0];
      if (!p) return;
      if (Math.hypot(p.x - p.sx, p.y - p.sy) > INPUT.longPressSlopCssPx) return;

      this.h.endGrab?.({ cancelled: true, longPress: true });
      this.mode = "context";
      this.h.longPress?.(this.pressHit, p.x, p.y);
    }, INPUT.longPressMs);
  }

  #cancelLongPress() {
    if (this.longTimer) clearTimeout(this.longTimer);
    this.longTimer = null;
  }

  reset() {
    this.#cancelLongPress();
    this.ptrs.clear();
    this.mode = "none";
    this.gesture = null;
    this.h.endGrab?.({ cancelled: true });
  }
}

function midpoint(a, b) {
  return { x: (a.x + b.x) * 0.5, y: (a.y + b.y) * 0.5 };
}

function distance(a, b) {
  return Math.hypot(b.x - a.x, b.y - a.y);
}

function angle(a, b) {
  return Math.atan2(b.y - a.y, b.x - a.x);
}
