export class Telemetry {
  constructor(bridge = globalThis.AndroidBridge) {
    this.bridge = bridge;
    this.samples = new Map();
    this.counters = new Map();
    this.gauges = new Map();
    this.faults = [];
    this.lastFlush = performance.now();
    this.flushEveryMs = 5000;
    this.maxSamples = 240;
  }

  sample(name, value) {
    if (!Number.isFinite(value)) return;
    let a = this.samples.get(name);
    if (!a) this.samples.set(name, (a = []));
    a.push(value);
    if (a.length > this.maxSamples) a.splice(0, a.length - this.maxSamples);
  }

  count(name, amount = 1) {
    this.counters.set(name, (this.counters.get(name) || 0) + amount);
  }

  gauge(name, value) {
    if (Number.isFinite(value)) this.gauges.set(name, value);
  }

  fault(kind, detail = {}) {
    const fault = {
      t: Math.round(performance.now()),
      kind,
      detail
    };
    this.faults.push(fault);
    if (this.faults.length > 40) this.faults.shift();
    this.bridge?.log?.("error", JSON.stringify(fault));
  }

  frame(frameMs) {
    this.sample("render.frame_ms", frameMs);
    const now = performance.now();
    if (now - this.lastFlush >= this.flushEveryMs) {
      this.flush();
      this.lastFlush = now;
    }
  }

  flush() {
    if (!this.bridge?.reportMetric) return;

    for (const [name, values] of this.samples) {
      if (!values.length) continue;
      const sorted = [...values].sort((a, b) => a - b);
      this.bridge.reportMetric(name + ".p50", percentile(sorted, 0.50));
      this.bridge.reportMetric(name + ".p95", percentile(sorted, 0.95));
    }

    for (const [name, value] of this.counters) {
      this.bridge.reportMetric(name, value);
    }
    for (const [name, value] of this.gauges) {
      this.bridge.reportMetric(name, value);
    }

    this.counters.clear();
  }

  installGlobalFaultHooks() {
    addEventListener("error", e => {
      this.fault("window_error", {
        message: String(e.message || "unknown"),
        file: String(e.filename || ""),
        line: e.lineno || 0,
        col: e.colno || 0
      });
    });

    addEventListener("unhandledrejection", e => {
      this.fault("unhandled_rejection", {
        reason: String(e.reason?.stack || e.reason || "unknown")
      });
    });
  }
}

function percentile(sorted, p) {
  if (!sorted.length) return 0;
  const i = Math.min(sorted.length - 1, Math.max(0, Math.round((sorted.length - 1) * p)));
  return sorted[i];
}
