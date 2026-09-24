# Pixel Physics Sandbox

Android-first 3D rigid-body sandbox MVP.

## Core interaction

- Touch a dynamic object to grab at the exact hit point.
- Drag moves a physical grab target on a camera-parallel plane.
- Pinch while grabbed changes depth.
- Two-finger twist applies torque.
- Release simply removes the grab force; existing rigid-body velocity produces the throw.
- Touch empty space to orbit; pinch empty space zooms/pans the camera.

## MVP systems

Spring-damper manipulator with force cap and partial mass compensation; Bullet rigid-body physics; 12 spawnable props across mahogany, metal and rubber profiles; freeze/unfreeze; delete; duplicate; material cycling; discrete undo; autosave/restore; procedural impact audio; haptics; low-resolution nearest-neighbor world rendering with native-resolution UI.

## Build

GitHub Actions builds and verifies `PixelPhysicsSandbox-MVP.apk` from the `pixel-physics-mvp` branch/PR. The app is offline after installation.

Package: `com.pixelphysics.sandbox`
Minimum Android: API 26
Target/compile SDK: API 35

CI additionally installs the APK on an Android 15 x86_64 emulator, launches the main activity, checks that the app process remains alive, checks the crash buffer, and captures a smoke-test screenshot using the SDK's platform-tools directly.
