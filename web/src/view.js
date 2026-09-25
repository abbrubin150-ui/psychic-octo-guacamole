import * as THREE from "three";
import { RENDER } from "./config.js";

export class OrthoRig {
  constructor(camera) {
    this.camera = camera;
    this.target = new THREE.Vector3(0, 0.8, 0);
    this.yaw = Math.PI / 4;
    this.pitch = Math.PI / 4;
    this.distance = 12;
    this.zoom = 1;
    this.update();
  }

  orbit(dx, dy) {
    this.yaw -= dx * 0.0036;
    this.pitch = THREE.MathUtils.clamp(
      this.pitch + dy * 0.0027,
      THREE.MathUtils.degToRad(34),
      THREE.MathUtils.degToRad(56)
    );
  }

  pan(dx, dy) {
    this.camera.updateMatrixWorld();
    const right = new THREE.Vector3().setFromMatrixColumn(this.camera.matrixWorld, 0);
    const up = new THREE.Vector3().setFromMatrixColumn(this.camera.matrixWorld, 1);
    const scale = 0.0065 / this.zoom;
    this.target.addScaledVector(right, -dx * scale);
    this.target.addScaledVector(up, dy * scale);
  }

  zoomByPixels(delta) {
    this.zoom *= Math.exp(delta * 0.006);
    this.zoom = THREE.MathUtils.clamp(this.zoom, 0.72, 2.6);
  }

  focusEntity(entity) {
    const t = entity.body.translation();
    this.target.set(t.x, t.y, t.z);
    this.zoom = Math.max(this.zoom, 1.25);
  }

  reset() {
    this.target.set(0, 0.8, 0);
    this.yaw = Math.PI / 4;
    this.pitch = Math.PI / 4;
    this.zoom = 1;
  }

  update() {
    const horizontal = Math.cos(this.pitch) * this.distance;
    this.camera.position.set(
      this.target.x + Math.sin(this.yaw) * horizontal,
      this.target.y + Math.sin(this.pitch) * this.distance,
      this.target.z + Math.cos(this.yaw) * horizontal
    );
    this.camera.lookAt(this.target);
    this.camera.zoom = this.zoom;
    this.camera.updateProjectionMatrix();
    this.camera.updateMatrixWorld();
  }

  snapshot() {
    return {
      target: this.target.toArray(),
      yaw: this.yaw,
      pitch: this.pitch,
      zoom: this.zoom
    };
  }

  restore(data) {
    if (!data) return;
    if (Array.isArray(data.target) && data.target.length === 3) this.target.fromArray(data.target);
    if (Number.isFinite(data.yaw)) this.yaw = data.yaw;
    if (Number.isFinite(data.pitch)) this.pitch = data.pitch;
    if (Number.isFinite(data.zoom)) this.zoom = data.zoom;
    this.update();
  }
}

export class PixelPipeline {
  constructor(renderer, scene, camera) {
    this.renderer = renderer;
    this.scene = scene;
    this.camera = camera;

    this.target = new THREE.WebGLRenderTarget(RENDER.width, RENDER.height, {
      minFilter: THREE.NearestFilter,
      magFilter: THREE.NearestFilter,
      depthBuffer: true,
      stencilBuffer: false
    });
    this.target.depthTexture = new THREE.DepthTexture(
      RENDER.width,
      RENDER.height,
      THREE.UnsignedIntType
    );

    this.postScene = new THREE.Scene();
    this.postCamera = new THREE.Camera();
    this.postMaterial = new THREE.ShaderMaterial({
      depthTest: false,
      depthWrite: false,
      uniforms: {
        tColor: { value: this.target.texture },
        tDepth: { value: this.target.depthTexture },
        resolution: { value: new THREE.Vector2(RENDER.width, RENDER.height) },
        levels: { value: RENDER.paletteLevels },
        edgeThreshold: { value: RENDER.outlineDepthThreshold }
      },
      vertexShader: `
        varying vec2 vUv;
        void main() {
          vUv = uv;
          gl_Position = vec4(position.xy, 0.0, 1.0);
        }
      `,
      fragmentShader: `
        precision highp float;
        uniform sampler2D tColor;
        uniform sampler2D tDepth;
        uniform vec2 resolution;
        uniform float levels;
        uniform float edgeThreshold;
        varying vec2 vUv;

        void main() {
          vec2 px = 1.0 / resolution;
          vec3 c = texture2D(tColor, vUv).rgb;
          float d = texture2D(tDepth, vUv).r;
          float d1 = texture2D(tDepth, vUv + vec2(px.x, 0.0)).r;
          float d2 = texture2D(tDepth, vUv + vec2(0.0, px.y)).r;
          float d3 = texture2D(tDepth, vUv - vec2(px.x, 0.0)).r;
          float d4 = texture2D(tDepth, vUv - vec2(0.0, px.y)).r;

          float edge = max(max(abs(d-d1), abs(d-d2)), max(abs(d-d3), abs(d-d4)));
          vec3 q = floor(c * (levels - 1.0) + 0.5) / (levels - 1.0);

          if (d < 0.99999 && edge > edgeThreshold) {
            q *= 0.38;
          }

          gl_FragColor = vec4(q, 1.0);
        }
      `
    });

    const quad = new THREE.Mesh(new THREE.PlaneGeometry(2, 2), this.postMaterial);
    quad.frustumCulled = false;
    this.postScene.add(quad);
  }

  render() {
    this.renderer.setRenderTarget(this.target);
    this.renderer.render(this.scene, this.camera);
    this.renderer.setRenderTarget(null);
    this.renderer.render(this.postScene, this.postCamera);
  }
}
