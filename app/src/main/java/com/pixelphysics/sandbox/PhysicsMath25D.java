package com.pixelphysics.sandbox;

/**
 * Allocation-free 2D footprint collision math for the 2.5D solver.
 * Shapes live in the world XY plane; Z is handled by the vertical contact layer.
 */
final class PhysicsMath25D {
    static final int BOX = 0;
    static final int CIRCLE = 1;
    static final float EPS = 1e-7f;

    static final class Shape {
        int kind = BOX;
        float x, y, yaw;
        float halfW, halfD, radius;

        Shape setBox(float x, float y, float yaw, float halfW, float halfD) {
            this.kind = BOX;
            this.x = x; this.y = y; this.yaw = yaw;
            this.halfW = halfW; this.halfD = halfD; this.radius = 0f;
            return this;
        }

        Shape setCircle(float x, float y, float radius) {
            this.kind = CIRCLE;
            this.x = x; this.y = y; this.yaw = 0f;
            this.halfW = this.halfD = 0f; this.radius = radius;
            return this;
        }
    }

    static final class Manifold {
        boolean hit;
        float nx, ny;
        float penetration;
        float cx, cy;

        void clear() {
            hit = false;
            nx = ny = penetration = cx = cy = 0f;
        }
    }

    static boolean collide(Shape a, Shape b, Manifold out) {
        out.clear();
        if (a.kind == CIRCLE && b.kind == CIRCLE) return circleCircle(a, b, out);
        if (a.kind == CIRCLE && b.kind == BOX) return circleBox(a, b, true, out);
        if (a.kind == BOX && b.kind == CIRCLE) return circleBox(b, a, false, out);
        return boxBox(a, b, out);
    }

    static float aabbHalfX(Shape s) {
        if (s.kind == CIRCLE) return s.radius;
        float c = Math.abs((float)Math.cos(s.yaw));
        float sn = Math.abs((float)Math.sin(s.yaw));
        return c * s.halfW + sn * s.halfD;
    }

    static float aabbHalfY(Shape s) {
        if (s.kind == CIRCLE) return s.radius;
        float c = Math.abs((float)Math.cos(s.yaw));
        float sn = Math.abs((float)Math.sin(s.yaw));
        return sn * s.halfW + c * s.halfD;
    }

    static boolean pointInside(Shape s, float px, float py, float margin) {
        if (s.kind == CIRCLE) {
            float dx = px - s.x, dy = py - s.y;
            float r = s.radius + margin;
            return dx * dx + dy * dy <= r * r;
        }
        float c = (float)Math.cos(s.yaw), sn = (float)Math.sin(s.yaw);
        float dx = px - s.x, dy = py - s.y;
        float lx = dx * c + dy * sn;
        float ly = -dx * sn + dy * c;
        return Math.abs(lx) <= s.halfW + margin && Math.abs(ly) <= s.halfD + margin;
    }

    private static boolean circleCircle(Shape a, Shape b, Manifold out) {
        float dx = b.x - a.x, dy = b.y - a.y;
        float rr = a.radius + b.radius;
        float d2 = dx * dx + dy * dy;
        if (d2 >= rr * rr) return false;

        float d = (float)Math.sqrt(Math.max(d2, EPS));
        float nx, ny;
        if (d2 < EPS) { nx = 1f; ny = 0f; }
        else { nx = dx / d; ny = dy / d; }

        out.hit = true;
        out.nx = nx; out.ny = ny;
        out.penetration = rr - d;
        float ax = a.x + nx * a.radius;
        float ay = a.y + ny * a.radius;
        float bx = b.x - nx * b.radius;
        float by = b.y - ny * b.radius;
        out.cx = (ax + bx) * 0.5f;
        out.cy = (ay + by) * 0.5f;
        return true;
    }

    /**
     * circle is the actual circle shape. circleIsA tells whether output normal should
     * point from circle->box (true) or box->circle (false).
     */
    private static boolean circleBox(Shape circle, Shape box, boolean circleIsA, Manifold out) {
        float c = (float)Math.cos(box.yaw), sn = (float)Math.sin(box.yaw);
        float dx = circle.x - box.x, dy = circle.y - box.y;

        float localX = dx * c + dy * sn;
        float localY = -dx * sn + dy * c;

        float qx = clamp(localX, -box.halfW, box.halfW);
        float qy = clamp(localY, -box.halfD, box.halfD);
        float ex = localX - qx, ey = localY - qy;
        float d2 = ex * ex + ey * ey;

        float nxLocal, nyLocal, penetration;
        float contactLocalX = qx, contactLocalY = qy;

        if (d2 > EPS) {
            float d = (float)Math.sqrt(d2);
            if (d >= circle.radius) return false;
            nxLocal = ex / d; nyLocal = ey / d;
            penetration = circle.radius - d;
        } else {
            // Circle center lies inside the box. Choose the nearest box face.
            float px = box.halfW - Math.abs(localX);
            float py = box.halfD - Math.abs(localY);
            if (px < py) {
                nxLocal = localX >= 0f ? 1f : -1f;
                nyLocal = 0f;
                contactLocalX = nxLocal * box.halfW;
                contactLocalY = localY;
                penetration = circle.radius + px;
            } else {
                nxLocal = 0f;
                nyLocal = localY >= 0f ? 1f : -1f;
                contactLocalX = localX;
                contactLocalY = nyLocal * box.halfD;
                penetration = circle.radius + py;
            }
        }

        // nxWorld currently points box -> circle.
        float nxWorld = nxLocal * c - nyLocal * sn;
        float nyWorld = nxLocal * sn + nyLocal * c;

        float cx = box.x + contactLocalX * c - contactLocalY * sn;
        float cy = box.y + contactLocalX * sn + contactLocalY * c;

        out.hit = true;
        if (circleIsA) {
            // Required normal is A(circle) -> B(box), the opposite direction.
            out.nx = -nxWorld; out.ny = -nyWorld;
        } else {
            out.nx = nxWorld; out.ny = nyWorld;
        }
        out.penetration = penetration;
        out.cx = cx; out.cy = cy;
        return true;
    }

    private static boolean boxBox(Shape a, Shape b, Manifold out) {
        float ca = (float)Math.cos(a.yaw), sa = (float)Math.sin(a.yaw);
        float cb = (float)Math.cos(b.yaw), sb = (float)Math.sin(b.yaw);

        float[] axesX = {ca, -sa, cb, -sb};
        float[] axesY = {sa,  ca, sb,  cb};

        float dx = b.x - a.x, dy = b.y - a.y;
        float bestOverlap = Float.POSITIVE_INFINITY;
        float bestNx = 0f, bestNy = 0f;

        for (int i = 0; i < 4; i++) {
            float ax = axesX[i], ay = axesY[i];
            float center = Math.abs(dx * ax + dy * ay);
            float ra = projectionRadius(a, ax, ay);
            float rb = projectionRadius(b, ax, ay);
            float overlap = ra + rb - center;
            if (overlap <= 0f) return false;

            if (overlap < bestOverlap) {
                bestOverlap = overlap;
                float sign = (dx * ax + dy * ay) >= 0f ? 1f : -1f;
                bestNx = ax * sign;
                bestNy = ay * sign;
            }
        }

        float[] pa = support(a, bestNx, bestNy);
        float[] pb = support(b, -bestNx, -bestNy);

        out.hit = true;
        out.nx = bestNx; out.ny = bestNy;
        out.penetration = bestOverlap;
        out.cx = (pa[0] + pb[0]) * 0.5f;
        out.cy = (pa[1] + pb[1]) * 0.5f;
        return true;
    }

    private static float projectionRadius(Shape s, float ax, float ay) {
        if (s.kind == CIRCLE) return s.radius;
        float c = (float)Math.cos(s.yaw), sn = (float)Math.sin(s.yaw);
        float ux = c, uy = sn;
        float vx = -sn, vy = c;
        return Math.abs(ax * ux + ay * uy) * s.halfW
                + Math.abs(ax * vx + ay * vy) * s.halfD;
    }

    private static float[] support(Shape s, float nx, float ny) {
        if (s.kind == CIRCLE) {
            float len = (float)Math.sqrt(nx * nx + ny * ny);
            if (len < EPS) return new float[]{s.x, s.y};
            return new float[]{s.x + nx / len * s.radius, s.y + ny / len * s.radius};
        }

        float c = (float)Math.cos(s.yaw), sn = (float)Math.sin(s.yaw);
        float ux = c, uy = sn;
        float vx = -sn, vy = c;
        float su = (nx * ux + ny * uy) >= 0f ? 1f : -1f;
        float sv = (nx * vx + ny * vy) >= 0f ? 1f : -1f;
        return new float[]{
                s.x + ux * s.halfW * su + vx * s.halfD * sv,
                s.y + uy * s.halfW * su + vy * s.halfD * sv
        };
    }

    private static float clamp(float v, float lo, float hi) {
        return Math.max(lo, Math.min(hi, v));
    }

    private PhysicsMath25D() {}
}
