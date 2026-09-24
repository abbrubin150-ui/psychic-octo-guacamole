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

CI additionally installs the APK on an Android 15 x86_64 emulator, launches the main activity, checks that the app process remains alive, checks the crash buffer, and captures a smoke-test screenshot using a managed Android 15 emulator runner.

The verified APK artifact is produced independently of the optional emulator smoke job, so CI infrastructure around AVD boot cannot block delivery of a structurally verified install package.


## True pixel-art renderer

The game now renders the entire world and HUD into a fixed 480x270 pixel canvas before nearest-neighbor integer upscaling. Mahogany, metal, rubber, workshop wood, floor and dark-metal surfaces are authored as deterministic 16x16/32x32 pixel textures rather than smooth PBR materials. UI panels, object markers, menus and text are rendered on the same pixel canvas so the visual language is consistent end-to-end.
