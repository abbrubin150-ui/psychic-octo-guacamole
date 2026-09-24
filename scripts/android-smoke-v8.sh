#!/usr/bin/env bash
set -euo pipefail

APK="${1:-PixelPhysicsLab-v8.1.apk}"
PKG="com.pixelphysics.sandbox"
ACT="$PKG/.AndroidLauncher"
LOG_TAG="PixelPhysicsLab"

adb install -r "$APK"
adb logcat -c
adb shell am force-stop "$PKG"
adb shell am start -n "$ACT" || true

PID=""
for i in $(seq 1 30); do
  PID="$(adb shell pidof "$PKG" 2>/dev/null | tr -d '\r' || true)"
  if [[ -n "$PID" ]]; then
    break
  fi
  sleep 1
done

if [[ -z "$PID" ]]; then
  echo "APP_PROCESS_MISSING"
  adb logcat -d -b crash || true
  adb logcat -d -v threadtime | tail -1400 || true
  exit 1
fi

READY=0
for i in $(seq 1 50); do
  LOG="$(adb logcat -d -s "$LOG_TAG:V" 2>/dev/null || true)"

  if grep -q "APP_READY Pixel Physics Lab 0.8.1" <<<"$LOG"; then
    READY=1
    break
  fi

  if grep -qE "BOOT_FATAL|WEB_ERROR|WEBVIEW_RENDERER_GONE" <<<"$LOG"; then
    echo "STARTUP_FAULT_DETECTED"
    printf '%s\n' "$LOG"
    exit 1
  fi

  PID="$(adb shell pidof "$PKG" 2>/dev/null | tr -d '\r' || true)"
  if [[ -z "$PID" ]]; then
    echo "APP_DIED_BEFORE_READY"
    adb logcat -d -b crash || true
    adb logcat -d -v threadtime | tail -1400 || true
    exit 1
  fi

  sleep 1
done

if [[ "$READY" != "1" ]]; then
  echo "APP_READY_TIMEOUT"
  adb logcat -d -s "$LOG_TAG:V" || true
  adb logcat -d -v threadtime | tail -1400 || true
  exit 1
fi

# Exercise camera and broad object-manipulation paths.
adb shell input swipe 450 450 760 390 650 || true
adb shell input swipe 820 500 1040 360 650 || true
adb shell input tap 600 420 || true
sleep 15

PID="$(adb shell pidof "$PKG" 2>/dev/null | tr -d '\r' || true)"
if [[ -z "$PID" ]]; then
  echo "APP_DIED_DURING_SOAK"
  adb logcat -d -b crash || true
  adb logcat -d -v threadtime | tail -1400 || true
  exit 1
fi

if adb logcat -d -b crash | grep -q "$PKG"; then
  echo "CRASH_BUFFER_MATCH"
  adb logcat -d -b crash
  exit 1
fi

if adb logcat -d -s "$LOG_TAG:E" | grep -qE "BOOT_FATAL|fault|WEB_ERROR|WEBVIEW_RENDERER_GONE"; then
  echo "RUNTIME_FAULT"
  adb logcat -d -s "$LOG_TAG:V"
  exit 1
fi

adb logcat -d -s "$LOG_TAG:V" > pixel-physics-lab-v81-metrics.txt || true
adb exec-out screencap -p > pixel-physics-lab-v81-smoke.png
test -s pixel-physics-lab-v81-smoke.png

echo "RUNTIME_PASS PID=$PID"
