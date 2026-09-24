package com.pixelphysics.sandbox;

import android.content.Context;
import android.content.SharedPreferences;
import android.graphics.Bitmap;
import android.graphics.Canvas;
import android.graphics.Color;
import android.graphics.Paint;
import android.graphics.Path;
import android.graphics.PointF;
import android.graphics.Rect;
import android.graphics.RectF;
import android.graphics.Typeface;
import android.os.SystemClock;
import android.view.HapticFeedbackConstants;
import android.view.MotionEvent;
import android.view.View;

import org.json.JSONArray;
import org.json.JSONObject;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

public final class PixelPhysics25DView extends View {
    private static final int W = 480;
    private static final int H = 270;
    private static final int FRAME = 80;

    private static final float ROOM = 4.3f;
    private static final float WALL_H = 3.6f;
    private static final float ORIGIN_X = 240f;
    private static final float ORIGIN_Y = 183f;
    private static final float ISO_X = 26f;
    private static final float ISO_Y = 10f;
    private static final float Z_PX = 24f;

    private static final float FIXED_DT = 1f / 60f;
    private static final float MAX_ACCUM = 0.12f;
    private static final float GRAVITY = 9.81f;
    private static final int MAX_PROPS = 100;

    private final Paint paint = new Paint();
    private final Paint spritePaint = new Paint();
    private final Bitmap logicalBitmap = Bitmap.createBitmap(W, H, Bitmap.Config.ARGB_8888);
    private final Canvas logical = new Canvas(logicalBitmap);
    private final Rect srcRect = new Rect(0, 0, W, H);
    private final RectF dstRect = new RectF();

    private float presentScale = 1f;
    private float presentX = 0f;
    private float presentY = 0f;
    private float presentW = W;
    private float presentH = H;

    private final ArrayList<Prop> props = new ArrayList<>();
    private final ArrayList<Prop> drawOrder = new ArrayList<>();
    private final HashMap<Integer, Prop> byId = new HashMap<>();
    private final HashMap<String, Bitmap[]> sprites = new HashMap<>();
    private final ArrayDeque<UndoAction> undo = new ArrayDeque<>();
    private final ArrayList<Platform> platforms = new ArrayList<>();
    private final ArrayList<RampSurface> ramps = new ArrayList<>();

    private int nextId = 1;
    private final SharedPreferences prefs;

    private long lastFrameNanos = 0L;
    private float accumulator = 0f;
    private long lastSaveMs = 0L;
    private boolean haptics = true;

    private enum Mode { IDLE, GRAB, CONTEXT, SPAWN, SETTINGS }
    private Mode mode = Mode.IDLE;

    private int primaryId = -1;
    private int secondaryId = -1;
    private float primaryX;
    private float primaryY;
    private float downX;
    private float downY;
    private long downNanos;
    private boolean contextTriggered;

    private Prop pressed;
    private Prop grabbed;
    private Prop contextProp;
    private float contextX;
    private float contextY;

    private float grabOffsetX;
    private float grabOffsetY;
    private float targetX;
    private float targetY;
    private float targetZ;
    private float lastTwoDistance;
    private float lastTwoAngle;

    private enum MaterialKind { MAHOGANY, METAL, RUBBER }

    private enum PropType {
        CUBE("CUBE", 1.0f, 1.0f, 1.0f, 1.0f, MaterialKind.MAHOGANY, true),
        BEAM_SHORT("BEAM S", 2.2f, 0.55f, 0.50f, 1.3f, MaterialKind.MAHOGANY, true),
        BEAM_LONG("BEAM L", 3.3f, 0.55f, 0.50f, 2.0f, MaterialKind.MAHOGANY, true),
        PLANK("PLANK", 2.7f, 1.0f, 0.28f, 1.2f, MaterialKind.MAHOGANY, true),
        WOOD_BALL("WOOD BALL", 0.9f, 0.9f, 0.9f, 0.7f, MaterialKind.MAHOGANY, false),
        METAL_BALL("METAL BALL", 0.9f, 0.9f, 0.9f, 4.0f, MaterialKind.METAL, false),
        WEIGHT("WEIGHT", 1.0f, 1.0f, 1.05f, 10.0f, MaterialKind.METAL, true),
        WHEEL("WHEEL", 1.25f, 0.42f, 1.25f, 1.1f, MaterialKind.MAHOGANY, false),
        RUBBER_BALL("RUBBER", 1.0f, 1.0f, 1.0f, 0.85f, MaterialKind.RUBBER, false),
        BARREL("BARREL", 1.1f, 1.1f, 1.5f, 3.6f, MaterialKind.METAL, false),
        CRATE("CRATE", 1.4f, 1.4f, 1.4f, 1.8f, MaterialKind.MAHOGANY, true),
        RAMP("RAMP", 2.4f, 1.4f, 0.8f, 1.8f, MaterialKind.MAHOGANY, false),
        SPRING("SPRING", 0.8f, 0.8f, 2.0f, 1.5f, MaterialKind.METAL, false);

        final String label;
        final float w;
        final float d;
        final float h;
        final float mass;
        final MaterialKind defaultMaterial;
        final boolean cuboid;

        PropType(String label, float w, float d, float h, float mass, MaterialKind defaultMaterial, boolean cuboid) {
            this.label = label;
            this.w = w;
            this.d = d;
            this.h = h;
            this.mass = mass;
            this.defaultMaterial = defaultMaterial;
            this.cuboid = cuboid;
        }
    }

    private static final class Palette {
        final int outline;
        final int dark;
        final int base;
        final int mid;
        final int light;
        final int hi;

        Palette(int outline, int dark, int base, int mid, int light, int hi) {
            this.outline = outline;
            this.dark = dark;
            this.base = base;
            this.mid = mid;
            this.light = light;
            this.hi = hi;
        }
    }

    private static final class Platform {
        final float cx;
        final float cy;
        final float w;
        final float d;
        final float top;

        Platform(float cx, float cy, float w, float d, float top) {
            this.cx = cx;
            this.cy = cy;
            this.w = w;
            this.d = d;
            this.top = top;
        }

        boolean contains(float x, float y, float margin) {
            return Math.abs(x - cx) <= w * 0.5f + margin && Math.abs(y - cy) <= d * 0.5f + margin;
        }
    }

    private static final class RampSurface {
        final float x0;
        final float x1;
        final float y0;
        final float y1;
        final float z0;
        final float z1;

        RampSurface(float x0, float x1, float y0, float y1, float z0, float z1) {
            this.x0 = x0;
            this.x1 = x1;
            this.y0 = y0;
            this.y1 = y1;
            this.z0 = z0;
            this.z1 = z1;
        }

        boolean contains(float x, float y, float margin) {
            return x >= Math.min(x0, x1) - margin && x <= Math.max(x0, x1) + margin
                    && y >= Math.min(y0, y1) - margin && y <= Math.max(y0, y1) + margin;
        }

        float heightAt(float x) {
            float t = clamp((x - x0) / (x1 - x0), 0f, 1f);
            return lerp(z0, z1, t);
        }

        float downhillSign() {
            return z1 > z0 ? -1f : 1f;
        }
    }

    private final class Prop {
        int id;
        PropType type;
        MaterialKind material;
        float x;
        float y;
        float z;
        float vx;
        float vy;
        float vz;
        float yaw;
        float spin;
        boolean frozen;

        float radius() {
            return Math.max(type.w, type.d) * 0.5f;
        }

        float bottom() {
            return z - type.h * 0.5f;
        }

        float top() {
            return z + type.h * 0.5f;
        }

        SaveState snapshot() {
            SaveState s = new SaveState();
            s.id = id;
            s.type = type.name();
            s.material = material.name();
            s.x = x;
            s.y = y;
            s.z = z;
            s.yaw = yaw;
            s.frozen = frozen;
            return s;
        }
    }

    private static final class SaveState {
        int id;
        String type;
        String material;
        float x;
        float y;
        float z;
        float yaw;
        boolean frozen;
    }

    private interface UndoAction {
        void undo();
    }

    public PixelPhysics25DView(Context context) {
        super(context);
        setFocusable(true);
        setFocusableInTouchMode(true);
        setKeepScreenOn(true);

        paint.setAntiAlias(false);
        paint.setDither(false);
        paint.setFilterBitmap(false);
        paint.setStrokeCap(Paint.Cap.SQUARE);
        paint.setStrokeJoin(Paint.Join.MITER);
        paint.setTypeface(Typeface.create(Typeface.MONOSPACE, Typeface.BOLD));

        spritePaint.setAntiAlias(false);
        spritePaint.setDither(false);
        spritePaint.setFilterBitmap(false);

        prefs = context.getSharedPreferences("pixel_physics_25d_canvas_v1", Context.MODE_PRIVATE);
        haptics = prefs.getBoolean("haptics", true);

        createSpriteSets();
        buildStaticPhysics();
        if (!restoreWorld()) {
            createStarterSet();
        }

        lastFrameNanos = System.nanoTime();
        lastSaveMs = SystemClock.uptimeMillis();
    }

    @Override
    protected void onSizeChanged(int w, int h, int oldw, int oldh) {
        updatePresentation(w, h);
    }

    private void updatePresentation(int viewW, int viewH) {
        if (viewW <= 0 || viewH <= 0) return;
        float raw = Math.min(viewW / (float) W, viewH / (float) H);
        float integer = (float) Math.floor(raw);
        presentScale = integer >= 1f ? integer : raw;
        presentW = W * presentScale;
        presentH = H * presentScale;
        presentX = (viewW - presentW) * 0.5f;
        presentY = (viewH - presentH) * 0.5f;
        dstRect.set(presentX, presentY, presentX + presentW, presentY + presentH);
    }

    @Override
    protected void onDraw(Canvas screen) {
        super.onDraw(screen);

        long now = System.nanoTime();
        float delta = Math.min(0.05f, (now - lastFrameNanos) / 1_000_000_000f);
        lastFrameNanos = now;

        checkLongPress(now);

        accumulator = Math.min(MAX_ACCUM, accumulator + delta);
        while (accumulator >= FIXED_DT) {
            updateGrab(FIXED_DT);
            physicsStep(FIXED_DT);
            accumulator -= FIXED_DT;
        }

        renderLogical();

        screen.drawColor(rgb("0B0E14"));
        spritePaint.setFilterBitmap(false);
        screen.drawBitmap(logicalBitmap, srcRect, dstRect, spritePaint);

        long ms = SystemClock.uptimeMillis();
        if (ms - lastSaveMs > 3000L) {
            lastSaveMs = ms;
            saveWorld();
        }

        postInvalidateOnAnimation();
    }

    private void renderLogical() {
        logical.drawColor(rgb("111822"));
        drawRoomBack(logical);
        drawPropShadows(logical);
        drawProps(logical);
        drawRoomFront(logical);
        drawHud(logical);
    }

    private Palette palette(MaterialKind kind) {
        switch (kind) {
            case METAL:
                return new Palette(rgb("101319"), rgb("272D34"), rgb("4A5660"), rgb("73818B"), rgb("AAB4BA"), rgb("E0E4E6"));
            case RUBBER:
                return new Palette(rgb("0F1511"), rgb("1A241C"), rgb("2B3B2E"), rgb("435746"), rgb("627A63"), rgb("A5B39D"));
            default:
                return new Palette(rgb("24100D"), rgb("421813"), rgb("6F271A"), rgb("963B21"), rgb("C25D30"), rgb("E88A49"));
        }
    }

    private static int rgb(String hex) {
        return Color.parseColor("#" + hex);
    }

    private void createSpriteSets() {
        for (PropType type : PropType.values()) {
            for (MaterialKind material : MaterialKind.values()) {
                Bitmap[] frames = new Bitmap[4];
                for (int dir = 0; dir < 4; dir++) {
                    frames[dir] = makeSprite(type, material, dir);
                }
                sprites.put(spriteKey(type, material), frames);
            }
        }
    }

    private String spriteKey(PropType type, MaterialKind material) {
        return type.name() + ":" + material.name();
    }

    private Bitmap makeSprite(PropType type, MaterialKind material, int dir) {
        Bitmap bmp = Bitmap.createBitmap(FRAME, FRAME, Bitmap.Config.ARGB_8888);
        Canvas c = new Canvas(bmp);
        Paint p = new Paint();
        p.setAntiAlias(false);
        p.setDither(false);
        p.setFilterBitmap(false);
        p.setStyle(Paint.Style.FILL);
        p.setStrokeCap(Paint.Cap.SQUARE);
        p.setStrokeJoin(Paint.Join.MITER);

        Palette pal = palette(material);

        if (type == PropType.WOOD_BALL || type == PropType.METAL_BALL || type == PropType.RUBBER_BALL) {
            drawBallSprite(c, p, type, pal, dir);
        } else if (type == PropType.WHEEL) {
            drawWheelSprite(c, p, pal, dir);
        } else if (type == PropType.RAMP) {
            drawRampSprite(c, p, pal, dir);
        } else if (type == PropType.SPRING) {
            drawSpringSprite(c, p, pal, dir);
        } else if (type == PropType.BARREL) {
            drawBarrelSprite(c, p, pal, dir);
        } else {
            drawCuboidSprite(c, p, type, pal, dir);
        }

        return bmp;
    }

    private void drawCuboidSprite(Canvas c, Paint p, PropType type, Palette pal, int dir) {
        float rw = dir % 2 == 0 ? type.w : type.d;
        float rd = dir % 2 == 0 ? type.d : type.w;
        int hx = Math.min(28, Math.max(5, Math.round(rw * 8.5f)));
        int hy = Math.min(14, Math.max(3, Math.round(rd * 4.0f)));
        int hh = Math.min(34, Math.max(5, Math.round(type.h * 15f)));

        int cx = FRAME / 2;
        int topY = 16;
        PointF a = new PointF(cx, topY);
        PointF b = new PointF(cx + hx, topY + hy);
        PointF cc = new PointF(cx, topY + hy * 2f);
        PointF d = new PointF(cx - hx, topY + hy);
        PointF b2 = new PointF(b.x, b.y + hh);
        PointF c2 = new PointF(cc.x, cc.y + hh);
        PointF d2 = new PointF(d.x, d.y + hh);

        int top = dir == 0 ? pal.light : dir == 1 ? pal.mid : dir == 2 ? pal.base : pal.hi;
        int left = dir == 3 ? pal.light : pal.base;
        int right = dir == 1 ? pal.light : pal.dark;

        fillPoly(c, p, left, d, cc, c2, d2);
        fillPoly(c, p, right, b, cc, c2, b2);
        fillPoly(c, p, top, a, b, cc, d);

        p.setColor(pal.outline);
        p.setStyle(Paint.Style.STROKE);
        p.setStrokeWidth(1f);
        drawPolyOutline(c, p, a, b, cc, d);
        c.drawLine(d.x, d.y, d2.x, d2.y, p);
        c.drawLine(cc.x, cc.y, c2.x, c2.y, p);
        c.drawLine(b.x, b.y, b2.x, b2.y, p);
        c.drawLine(d2.x, d2.y, c2.x, c2.y, p);
        c.drawLine(c2.x, c2.y, b2.x, b2.y, p);
        p.setStyle(Paint.Style.FILL);

        p.setColor(pal.mid);
        if (type == PropType.CRATE) {
            p.setStyle(Paint.Style.STROKE);
            c.drawLine(d.x + 3, d.y + 4, c2.x - 3, c2.y - 4, p);
            c.drawLine(d2.x + 3, d2.y - 4, cc.x - 3, cc.y + 4, p);
            c.drawLine(b.x - 3, b.y + 4, c2.x + 3, c2.y - 4, p);
            c.drawLine(b2.x - 3, b2.y - 4, cc.x + 3, cc.y + 4, p);
            p.setStyle(Paint.Style.FILL);
        } else {
            int sy = Math.min((int) c2.y - 4, (int) d.y + 7 + dir * 2);
            c.drawRect(d.x + 3, sy, cc.x - 2, sy + 1, p);
        }

        p.setColor(pal.hi);
        float markX = (dir == 1 || dir == 2) ? b.x - 4 : d.x + 3;
        float markY = dir >= 2 ? Math.min(d2.y - 5, d.y + hh - 5) : d.y + 5;
        c.drawRect(markX, markY, markX + 2, markY + 2, p);
    }

    private void drawBallSprite(Canvas c, Paint p, PropType type, Palette pal, int dir) {
        int r = Math.max(8, Math.round(type.w * 13f));
        int cx = FRAME / 2;
        int cy = FRAME / 2 + 5;

        p.setColor(pal.outline);
        c.drawCircle(cx, cy, r + 1, p);
        p.setColor(pal.dark);
        c.drawCircle(cx, cy, r, p);
        p.setColor(pal.base);
        c.drawCircle(cx - 1, cy - 1, r - 2, p);
        p.setColor(pal.mid);
        c.drawCircle(cx - 2, cy - 2, Math.max(3, r - 5), p);

        int[][] h = {{-5, -5}, {5, -5}, {5, 5}, {-5, 5}};
        p.setColor(pal.light);
        c.drawRect(cx + h[dir][0] - 2, cy + h[dir][1] - 2, cx + h[dir][0] + 3, cy + h[dir][1] + 2, p);
        p.setColor(pal.hi);
        c.drawRect(cx + h[dir][0] - 1, cy + h[dir][1] - 1, cx + h[dir][0] + 1, cy + h[dir][1] + 1, p);

        p.setColor(pal.outline);
        p.setStrokeWidth(1f);
        if (dir % 2 == 0) {
            c.drawLine(cx - r + 3, cy + 3, cx + r - 3, cy - 3, p);
        } else {
            c.drawLine(cx - 3, cy - r + 3, cx + 3, cy + r - 3, p);
        }
    }

    private void drawWheelSprite(Canvas c, Paint p, Palette pal, int dir) {
        int cx = FRAME / 2;
        int cy = FRAME / 2 + 4;
        int r = 17;

        p.setColor(pal.outline);
        c.drawCircle(cx, cy, r + 1, p);
        p.setColor(pal.dark);
        c.drawCircle(cx, cy, r, p);
        p.setColor(pal.base);
        c.drawCircle(cx, cy, r - 4, p);
        p.setColor(pal.dark);
        c.drawCircle(cx, cy, 5, p);
        p.setColor(pal.hi);
        c.drawCircle(cx, cy, 2, p);

        p.setColor(pal.light);
        p.setStrokeWidth(2f);
        float phase = dir * ((float) Math.PI / 8f);
        for (int i = 0; i < 8; i++) {
            float a = phase + i * ((float) Math.PI / 4f);
            float x1 = cx + (float) Math.cos(a) * 5f;
            float y1 = cy + (float) Math.sin(a) * 5f;
            float x2 = cx + (float) Math.cos(a) * (r - 5);
            float y2 = cy + (float) Math.sin(a) * (r - 5);
            c.drawLine(x1, y1, x2, y2, p);
        }
        p.setStrokeWidth(1f);

        p.setColor(pal.hi);
        float mx = cx - 7 + (dir % 2) * 10;
        c.drawRect(mx, cy - r + 4, mx + 3, cy - r + 6, p);
    }

    private void drawRampSprite(Canvas c, Paint p, Palette pal, int dir) {
        int cx = FRAME / 2;
        int baseY = 61;
        int left = cx - 28;
        int right = cx + 28;
        int topY = 24;
        boolean flip = dir == 2 || dir == 3;

        p.setColor(pal.outline);
        Path outline = new Path();
        if (!flip) {
            outline.moveTo(left - 1, baseY + 1);
            outline.lineTo(right + 1, baseY + 1);
            outline.lineTo(right + 1, topY - 1);
        } else {
            outline.moveTo(left - 1, baseY + 1);
            outline.lineTo(right + 1, baseY + 1);
            outline.lineTo(left - 1, topY - 1);
        }
        outline.close();
        c.drawPath(outline, p);

        p.setColor(pal.base);
        Path body = new Path();
        if (!flip) {
            body.moveTo(left, baseY);
            body.lineTo(right, baseY);
            body.lineTo(right, topY);
        } else {
            body.moveTo(left, baseY);
            body.lineTo(right, baseY);
            body.lineTo(left, topY);
        }
        body.close();
        c.drawPath(body, p);

        p.setColor(pal.light);
        p.setStrokeWidth(2f);
        if (!flip) c.drawLine(left + 4, baseY - 3, right - 3, topY + 4, p);
        else c.drawLine(left + 3, topY + 4, right - 4, baseY - 3, p);

        p.setColor(pal.dark);
        c.drawLine(left + 4, baseY - 5, right - 4, baseY - 5, p);
        p.setStrokeWidth(1f);

        p.setColor(pal.hi);
        float mx = flip ? left + 7 : right - 9;
        c.drawRect(mx, baseY - 10, mx + 3, baseY - 7, p);
    }

    private void drawSpringSprite(Canvas c, Paint p, Palette pal, int dir) {
        int cx = FRAME / 2;
        p.setColor(pal.outline);
        c.drawRect(cx - 9, 11, cx + 9, 16, p);
        c.drawRect(cx - 11, 64, cx + 11, 69, p);
        p.setColor(pal.base);
        c.drawRect(cx - 7, 12, cx + 7, 15, p);
        c.drawRect(cx - 9, 65, cx + 9, 68, p);

        int shift = dir % 2 == 0 ? 0 : 2;
        p.setStrokeWidth(2f);
        for (int y = 20; y < 60; y += 5) {
            p.setColor((y / 5) % 2 == 0 ? pal.light : pal.mid);
            if (((y / 5) + dir) % 2 == 0) c.drawLine(cx - 9 + shift, y, cx + 9 - shift, y + 3, p);
            else c.drawLine(cx + 9 - shift, y, cx - 9 + shift, y + 3, p);
        }
        p.setStrokeWidth(1f);

        p.setColor(pal.hi);
        c.drawRect(cx - 2 + (dir - 1), 14, cx + 1 + (dir - 1), 16, p);
    }

    private void drawBarrelSprite(Canvas c, Paint p, Palette pal, int dir) {
        int cx = FRAME / 2;
        int cy = 39;
        int w = 23;
        int h = 36;

        p.setColor(pal.outline);
        c.drawRect(cx - w / 2f - 1, cy - h / 2f + 3, cx + w / 2f + 1, cy + h / 2f - 3, p);
        c.drawCircle(cx, cy - h / 2f + 4, w / 2f + 1, p);
        c.drawCircle(cx, cy + h / 2f - 4, w / 2f + 1, p);

        p.setColor(pal.base);
        c.drawRect(cx - w / 2f, cy - h / 2f + 4, cx + w / 2f, cy + h / 2f - 4, p);
        p.setColor(pal.mid);
        c.drawCircle(cx, cy - h / 2f + 4, w / 2f - 1, p);
        c.drawCircle(cx, cy + h / 2f - 4, w / 2f - 1, p);

        p.setColor(pal.dark);
        c.drawRect(cx - w / 2f, cy - 8, cx + w / 2f, cy - 5, p);
        c.drawRect(cx - w / 2f, cy + 6, cx + w / 2f, cy + 9, p);

        p.setColor(pal.light);
        float lx = dir < 2 ? cx - w / 2f + 4 : cx + w / 2f - 6;
        c.drawRect(lx, cy - h / 2f + 7, lx + 3, cy + h / 2f - 7, p);

        p.setColor(pal.hi);
        c.drawRect(lx, cy - h / 2f + 8, lx + 2, cy - h / 2f + 11, p);
    }

    private void fillPoly(Canvas c, Paint p, int color, PointF... pts) {
        Path path = new Path();
        path.moveTo(pts[0].x, pts[0].y);
        for (int i = 1; i < pts.length; i++) path.lineTo(pts[i].x, pts[i].y);
        path.close();
        p.setColor(color);
        p.setStyle(Paint.Style.FILL);
        c.drawPath(path, p);
    }

    private void drawPolyOutline(Canvas c, Paint p, PointF... pts) {
        Path path = new Path();
        path.moveTo(pts[0].x, pts[0].y);
        for (int i = 1; i < pts.length; i++) path.lineTo(pts[i].x, pts[i].y);
        path.close();
        c.drawPath(path, p);
    }

    private void buildStaticPhysics() {
        platforms.add(new Platform(2.1f, -0.5f, 1.5f, 1.5f, 0.72f));
        platforms.add(new Platform(-2.8f, 1.0f, 1.7f, 1.4f, 0.80f));
        ramps.add(new RampSurface(-3.3f, 0.1f, -1.5f, -0.25f, 1.35f, 0f));
    }

    private void createStarterSet() {
        spawnInternal(PropType.WOOD_BALL, MaterialKind.MAHOGANY, -2.4f, -0.85f, 1.9f, 0f, false, nextId++);
        spawnInternal(PropType.CRATE, MaterialKind.MAHOGANY, -3.0f, 1.0f, 1.55f, 0f, false, nextId++);
        spawnInternal(PropType.WHEEL, MaterialKind.MAHOGANY, -1.7f, -2.4f, 0.75f, 0.1f, false, nextId++);
        spawnInternal(PropType.CUBE, MaterialKind.METAL, -0.6f, -2.5f, 0.55f, 0f, false, nextId++);
        spawnInternal(PropType.CUBE, MaterialKind.MAHOGANY, 2.1f, -0.5f, 1.72f, 0f, false, nextId++);
        spawnInternal(PropType.CUBE, MaterialKind.RUBBER, 2.1f, -0.5f, 2.72f, 0f, false, nextId++);
        spawnInternal(PropType.CUBE, MaterialKind.MAHOGANY, 2.1f, -0.5f, 3.72f, 0f, false, nextId++);
        spawnInternal(PropType.BEAM_SHORT, MaterialKind.MAHOGANY, 1.0f, -2.3f, 0.55f, 0.15f, false, nextId++);
        spawnInternal(PropType.WEIGHT, MaterialKind.METAL, 0.7f, 2.6f, 2.2f, 0f, true, nextId++);
        saveWorld();
    }

    private Prop spawnInternal(PropType type, MaterialKind material, float x, float y, float z, float yaw, boolean frozen, int id) {
        if (props.size() >= MAX_PROPS) return null;
        Prop p = new Prop();
        p.id = id > 0 ? id : nextId++;
        nextId = Math.max(nextId, p.id + 1);
        p.type = type;
        p.material = material;
        p.x = x;
        p.y = y;
        p.z = z;
        p.yaw = yaw;
        p.frozen = frozen;
        props.add(p);
        byId.put(p.id, p);
        return p;
    }

    private void removeProp(Prop p) {
        if (p == null) return;
        if (grabbed == p) endGrab();
        props.remove(p);
        byId.remove(p.id);
    }

    private void spawnWithUndo(PropType type) {
        float n = (props.size() % 5) - 2f;
        Prop p = spawnInternal(type, type.defaultMaterial, n * 0.35f, 0.1f, 3.3f, 0f, false, nextId++);
        if (p == null) return;
        final int id = p.id;
        pushUndo(() -> {
            Prop q = byId.get(id);
            if (q != null) removeProp(q);
        });
        feedback();
        saveWorld();
    }

    private void deleteWithUndo(Prop p) {
        if (p == null) return;
        SaveState s = p.snapshot();
        pushUndo(() -> spawnFromState(s));
        removeProp(p);
        feedback();
        saveWorld();
    }

    private void duplicateWithUndo(Prop source) {
        if (source == null) return;
        Prop p = spawnInternal(source.type, source.material, source.x + 0.35f, source.y - 0.35f, source.z + 0.45f,
                source.yaw, source.frozen, nextId++);
        if (p == null) return;
        final int id = p.id;
        pushUndo(() -> {
            Prop q = byId.get(id);
            if (q != null) removeProp(q);
        });
        feedback();
        saveWorld();
    }

    private void setFrozen(Prop p, boolean frozen, boolean record) {
        if (p == null || p.frozen == frozen) return;
        final int id = p.id;
        final boolean prior = p.frozen;
        if (record) {
            pushUndo(() -> {
                Prop q = byId.get(id);
                if (q != null) setFrozen(q, prior, false);
            });
        }
        p.frozen = frozen;
        p.vx = p.vy = p.vz = p.spin = 0f;
        feedback();
        saveWorld();
    }

    private void cycleMaterial(Prop p) {
        if (p == null) return;
        final int id = p.id;
        final MaterialKind prior = p.material;
        MaterialKind next = prior == MaterialKind.MAHOGANY ? MaterialKind.METAL
                : prior == MaterialKind.METAL ? MaterialKind.RUBBER : MaterialKind.MAHOGANY;
        pushUndo(() -> {
            Prop q = byId.get(id);
            if (q != null) q.material = prior;
        });
        p.material = next;
        feedback();
        saveWorld();
    }

    private void pushUndo(UndoAction action) {
        undo.addLast(action);
        while (undo.size() > 32) undo.removeFirst();
    }

    private void doUndo() {
        if (undo.isEmpty()) return;
        undo.removeLast().undo();
        feedback();
        saveWorld();
    }

    private void physicsStep(float dt) {
        for (Prop p : props) {
            if (p.frozen) continue;

            if (p != grabbed) p.vz -= GRAVITY * dt;

            p.x += p.vx * dt;
            p.y += p.vy * dt;
            p.z += p.vz * dt;
            p.yaw += p.spin * dt;
            p.spin *= 0.992f;

            float r = p.radius();

            if (p.x - r < -ROOM) {
                p.x = -ROOM + r;
                p.vx = Math.abs(p.vx) * restitution(p.material);
            }
            if (p.x + r > ROOM) {
                p.x = ROOM - r;
                p.vx = -Math.abs(p.vx) * restitution(p.material);
            }
            if (p.y - r < -ROOM) {
                p.y = -ROOM + r;
                p.vy = Math.abs(p.vy) * restitution(p.material);
            }
            if (p.y + r > ROOM) {
                p.y = ROOM - r;
                p.vy = -Math.abs(p.vy) * restitution(p.material);
            }

            float support = supportHeight(p);
            if (p.bottom() < support) {
                p.z = support + p.type.h * 0.5f;

                if (p.vz < 0f) p.vz = -p.vz * restitution(p.material);
                if (Math.abs(p.vz) < 0.35f) p.vz = 0f;

                float fr = friction(p.material);
                p.vx *= Math.max(0f, 1f - fr * dt * 2.4f);
                p.vy *= Math.max(0f, 1f - fr * dt * 2.4f);

                RampSurface ramp = rampUnder(p);
                if (ramp != null && Math.abs(p.vz) < 0.5f) {
                    p.vx += ramp.downhillSign() * 3.0f * dt;
                }
            }

            if (p.z < -3f) {
                p.x = 0f;
                p.y = 0f;
                p.z = 4f;
                p.vx = p.vy = p.vz = 0f;
            }
        }

        solveHorizontalCollisions();
    }

    private void solveHorizontalCollisions() {
        for (int i = 0; i < props.size(); i++) {
            Prop a = props.get(i);
            for (int j = i + 1; j < props.size(); j++) {
                Prop b = props.get(j);

                if (a.frozen && b.frozen) continue;
                if (a.top() < b.bottom() + 0.04f || b.top() < a.bottom() + 0.04f) continue;

                float dx = b.x - a.x;
                float dy = b.y - a.y;
                float min = (a.radius() + b.radius()) * 0.72f;
                float d2 = dx * dx + dy * dy;

                if (d2 >= min * min) continue;

                float d = (float) Math.sqrt(Math.max(d2, 0.0001f));
                float nx = dx / d;
                float ny = dy / d;
                float penetration = min - d;

                float wa = a.frozen ? 0f : 1f;
                float wb = b.frozen ? 0f : 1f;
                float sum = wa + wb;
                if (sum <= 0f) continue;

                if (!a.frozen) {
                    a.x -= nx * penetration * (wa / sum);
                    a.y -= ny * penetration * (wa / sum);
                }
                if (!b.frozen) {
                    b.x += nx * penetration * (wb / sum);
                    b.y += ny * penetration * (wb / sum);
                }

                float rvx = b.vx - a.vx;
                float rvy = b.vy - a.vy;
                float rel = rvx * nx + rvy * ny;

                if (rel < 0f) {
                    float e = Math.min(restitution(a.material), restitution(b.material));
                    float invA = a.frozen ? 0f : 1f / a.type.mass;
                    float invB = b.frozen ? 0f : 1f / b.type.mass;
                    float impulse = -(1f + e) * rel / Math.max(0.0001f, invA + invB);

                    if (!a.frozen) {
                        a.vx -= impulse * nx * invA;
                        a.vy -= impulse * ny * invA;
                    }
                    if (!b.frozen) {
                        b.vx += impulse * nx * invB;
                        b.vy += impulse * ny * invB;
                    }
                }
            }
        }
    }

    private float supportHeight(Prop p) {
        float best = 0f;
        float margin = p.radius() * 0.35f;

        for (Platform platform : platforms) {
            if (platform.contains(p.x, p.y, margin)) best = Math.max(best, platform.top);
        }

        for (RampSurface ramp : ramps) {
            if (ramp.contains(p.x, p.y, margin)) best = Math.max(best, ramp.heightAt(p.x));
        }

        for (Prop q : props) {
            if (q == p) continue;

            float dx = p.x - q.x;
            float dy = p.y - q.y;
            float rr = (p.radius() + q.radius()) * 0.72f;

            if (dx * dx + dy * dy > rr * rr) continue;

            float top = q.top();
            if (top <= p.z + 0.18f && top > best) best = top;
        }

        return best;
    }

    private RampSurface rampUnder(Prop p) {
        for (RampSurface ramp : ramps) {
            if (ramp.contains(p.x, p.y, p.radius() * 0.2f)) return ramp;
        }
        return null;
    }

    private float restitution(MaterialKind m) {
        return m == MaterialKind.RUBBER ? 0.72f : m == MaterialKind.METAL ? 0.12f : 0.16f;
    }

    private float friction(MaterialKind m) {
        return m == MaterialKind.RUBBER ? 0.85f : m == MaterialKind.METAL ? 0.45f : 0.68f;
    }

    private void updateGrab(float dt) {
        if (mode != Mode.GRAB || grabbed == null || grabbed.frozen) return;

        float mass = Math.max(0.25f, grabbed.type.mass);
        float massScale = (float) Math.pow(mass, 0.34);
        float kp = 55f * massScale;
        float kd = 9f * (float) Math.sqrt(massScale);

        float fx = (targetX - grabbed.x) * kp - grabbed.vx * kd;
        float fy = (targetY - grabbed.y) * kp - grabbed.vy * kd;
        float fz = (targetZ - grabbed.z) * kp - grabbed.vz * kd;

        float cap = 95f * (float) Math.pow(Math.max(1f, mass), 0.58);
        float len = (float) Math.sqrt(fx * fx + fy * fy + fz * fz);

        if (len > cap) {
            float s = cap / len;
            fx *= s;
            fy *= s;
            fz *= s;
        }

        grabbed.vx += fx / mass * dt;
        grabbed.vy += fy / mass * dt;
        grabbed.vz += fz / mass * dt;
    }

    private PointF project(float x, float y, float z) {
        return new PointF(
                ORIGIN_X + (x - y) * ISO_X,
                ORIGIN_Y - (x + y) * ISO_Y - z * Z_PX
        );
    }

    private PointF unprojectAtZ(float sx, float sy, float z) {
        float a = (sx - ORIGIN_X) / ISO_X;
        float b = (ORIGIN_Y - sy - z * Z_PX) / ISO_Y;
        return new PointF((a + b) * 0.5f, (b - a) * 0.5f);
    }

    private int directionIndex(float yaw) {
        int q = Math.round((float) Math.toDegrees(yaw) / 90f) % 4;
        if (q < 0) q += 4;
        return q;
    }

    private float depthKey(Prop p) {
        return ORIGIN_Y - (p.x + p.y) * ISO_Y - p.z * Z_PX;
    }

    private void drawRoomBack(Canvas c) {
        paint.setStyle(Paint.Style.FILL);
        paint.setStrokeWidth(1f);
        paint.setAntiAlias(false);

        c.drawColor(rgb("111822"));

        PointF front = project(-ROOM, -ROOM, 0f);
        PointF right = project(ROOM, -ROOM, 0f);
        PointF back = project(ROOM, ROOM, 0f);
        PointF left = project(-ROOM, ROOM, 0f);

        PointF rightTop = project(ROOM, -ROOM, WALL_H);
        PointF backTop = project(ROOM, ROOM, WALL_H);
        PointF leftTop = project(-ROOM, ROOM, WALL_H);

        fillPoly(c, paint, rgb("56433F"), right, back, backTop, rightTop);
        fillPoly(c, paint, rgb("5D4640"), left, back, backTop, leftTop);

        paint.setColor(rgb("6B5650"));
        paint.setStrokeWidth(1f);
        for (float z = 0.45f; z < WALL_H; z += 0.72f) {
            line(c, project(ROOM, -ROOM, z), project(ROOM, ROOM, z), 1f, rgb("6B5650"));
            line(c, project(-ROOM, ROOM, z), project(ROOM, ROOM, z), 1f, rgb("6B5650"));
        }

        fillPoly(c, paint, rgb("AD7A4E"), front, right, back, left);

        for (float v = -ROOM; v <= ROOM + 0.01f; v += 0.86f) {
            line(c, project(v, -ROOM, 0.01f), project(v, ROOM, 0.01f), 1f, rgb("7D583C"));
            line(c, project(-ROOM, v, 0.01f), project(ROOM, v, 0.01f), 1f, rgb("7D583C"));
        }

        drawWallRailY(c, ROOM, 3.18f, 5f, rgb("6A2F1D"));
        drawWallRailX(c, ROOM, 3.18f, 5f, rgb("6A2F1D"));
        drawWallRailY(c, ROOM, 0.25f, 4f, rgb("4B241A"));
        drawWallRailX(c, ROOM, 0.25f, 4f, rgb("4B241A"));

        drawVerticalPost(c, ROOM, ROOM, 0f, WALL_H, 9f, rgb("4A2118"), rgb("A55730"));
        drawVerticalPost(c, -ROOM, ROOM, 0f, WALL_H, 7f, rgb("6B301D"), rgb("A55730"));
        drawVerticalPost(c, ROOM, -ROOM, 0f, WALL_H, 7f, rgb("6B301D"), rgb("A55730"));

        drawWallRectY(c, ROOM, -3.35f, -1.95f, 1.35f, 2.65f, rgb("3A231C"));
        drawWallRectY(c, ROOM, -3.22f, -2.08f, 1.48f, 2.52f, rgb("F0B85D"));
        drawWallLineY(c, ROOM, -2.65f, 1.48f, -2.65f, 2.52f, 3f, rgb("5A2A1B"));
        drawWallLineY(c, ROOM, -3.22f, 2.0f, -2.08f, 2.0f, 3f, rgb("5A2A1B"));

        drawWallRectY(c, ROOM, -1.55f, 0.45f, 1.12f, 2.35f, rgb("202B3A"));
        drawWallLineY(c, ROOM, -1.35f, 1.34f, 0.18f, 2.05f, 2f, rgb("8A8580"));
        drawWallLineY(c, ROOM, -1.20f, 1.48f, -0.55f, 1.22f, 2f, rgb("8A8580"));

        drawShelfY(c, ROOM, -0.8f, 0.75f, 2.65f);
        drawShelfX(c, ROOM, -1.4f, 0.6f, 2.45f);

        PointF lampTop = project(-1.95f, 3.9f, 3.45f);
        paint.setColor(rgb("2D241F"));
        c.drawRect(lampTop.x - 1, lampTop.y, lampTop.x + 1, lampTop.y + 24, paint);
        paint.setColor(rgb("C06C24"));
        c.drawRect(lampTop.x - 7, lampTop.y + 24, lampTop.x + 7, lampTop.y + 29, paint);
        paint.setColor(rgb("FFD66B"));
        c.drawRect(lampTop.x - 5, lampTop.y + 29, lampTop.x + 5, lampTop.y + 40, paint);
        paint.setColor(rgb("FFF0A8"));
        c.drawRect(lampTop.x - 2, lampTop.y + 31, lampTop.x + 2, lampTop.y + 38, paint);

        PointF spring = project(3.92f, 1.55f, 2.95f);
        paint.setColor(rgb("20252B"));
        c.drawRect(spring.x - 5, spring.y - 2, spring.x + 5, spring.y + 3, paint);
        for (int i = 0; i < 8; i++) {
            paint.setColor(i % 2 == 0 ? rgb("AEB8BF") : rgb("59636B"));
            float yy = spring.y + 7 + i * 6;
            if (i % 2 == 0) c.drawRect(spring.x - 7, yy, spring.x + 7, yy + 3, paint);
            else c.drawRect(spring.x - 4, yy, spring.x + 4, yy + 3, paint);
        }
        paint.setColor(rgb("20252B"));
        c.drawRect(spring.x - 1, spring.y + 55, spring.x + 1, spring.y + 65, paint);

        for (Platform platform : platforms) drawPlatform(c, platform);
        for (RampSurface ramp : ramps) drawStaticRamp(c, ramp);

        PointF ropeTop = project(0.7f, 2.6f, WALL_H);
        PointF ropeBottom = project(0.7f, 2.6f, 2.72f);
        line(c, ropeTop, ropeBottom, 3f, rgb("9B5D2B"));
    }

    private void drawRoomFront(Canvas c) {
        PointF front = project(-ROOM, -ROOM, 0f);
        PointF left = project(-ROOM, ROOM, 0f);
        PointF right = project(ROOM, -ROOM, 0f);

        line(c, left, front, 9f, rgb("5B2A1D"));
        line(c, front, right, 9f, rgb("5B2A1D"));
        line(c, left, front, 2f, rgb("B25F34"));
        line(c, front, right, 2f, rgb("B25F34"));

        drawCornerCap(c, left.x, left.y, false);
        drawCornerCap(c, front.x, front.y - 2, true);
        drawCornerCap(c, right.x, right.y, false);
    }

    private void drawCornerCap(Canvas c, float x, float y, boolean large) {
        float s = large ? 9f : 7f;
        paint.setStyle(Paint.Style.FILL);
        paint.setColor(rgb("343844"));
        c.drawRect(x - s, y - s, x + s, y + s, paint);
        paint.setColor(rgb("707783"));
        c.drawRect(x - s + 2, y - s + 2, x + s - 2, y - s + 4, paint);
        paint.setColor(rgb("151820"));
        c.drawRect(x - s + 3, y - 1, x - s + 5, y + 1, paint);
        c.drawRect(x + s - 5, y - 1, x + s - 3, y + 1, paint);
    }

    private void drawPlatform(Canvas c, Platform p) {
        PointF a = project(p.cx - p.w / 2f, p.cy - p.d / 2f, p.top);
        PointF b = project(p.cx + p.w / 2f, p.cy - p.d / 2f, p.top);
        PointF cc = project(p.cx + p.w / 2f, p.cy + p.d / 2f, p.top);
        PointF d = project(p.cx - p.w / 2f, p.cy + p.d / 2f, p.top);
        fillPoly(c, paint, rgb("7A3C23"), a, b, cc, d);
        line(c, d, cc, 2f, rgb("D17A3C"));

        PointF base = project(p.cx, p.cy, 0f);
        PointF top = project(p.cx, p.cy, p.top);
        paint.setColor(rgb("4A2419"));
        c.drawRect(base.x - 4, top.y, base.x + 4, base.y, paint);
    }

    private void drawStaticRamp(Canvas c, RampSurface r) {
        PointF a = project(r.x0, r.y0, r.z0);
        PointF b = project(r.x1, r.y0, r.z1);
        PointF cc = project(r.x1, r.y1, r.z1);
        PointF d = project(r.x0, r.y1, r.z0);
        fillPoly(c, paint, rgb("B65E2E"), a, b, cc, d);
        line(c, d, cc, 2f, rgb("E58A48"));
        line(c, a, b, 2f, rgb("552619"));

        for (int i = 1; i < 5; i++) {
            float t = i / 5f;
            float x = lerp(r.x0, r.x1, t);
            float z = lerp(r.z0, r.z1, t);
            line(c, project(x, r.y0, z + 0.01f), project(x, r.y1, z + 0.01f), 1f, rgb("6F351F"));
        }
    }

    private void drawPropShadows(Canvas c) {
        paint.setStyle(Paint.Style.FILL);
        paint.setColor(Color.argb(65, 25, 18, 20));

        for (Prop p : props) {
            float support = Math.max(0f, supportHeightIgnoringProp(p));
            PointF s = project(p.x, p.y, support + 0.01f);
            float rx = Math.max(6f, p.radius() * 12f);
            float ry = Math.max(3f, p.radius() * 4.5f);
            c.drawOval(s.x - rx, s.y - ry, s.x + rx, s.y + ry, paint);
        }
    }

    private float supportHeightIgnoringProp(Prop p) {
        float best = 0f;
        float margin = p.radius() * 0.35f;

        for (Platform platform : platforms) {
            if (platform.contains(p.x, p.y, margin)) best = Math.max(best, platform.top);
        }
        for (RampSurface ramp : ramps) {
            if (ramp.contains(p.x, p.y, margin)) best = Math.max(best, ramp.heightAt(p.x));
        }
        return best;
    }

    private void drawProps(Canvas c) {
        drawOrder.clear();
        drawOrder.addAll(props);
        drawOrder.sort(Comparator.comparingDouble(this::depthKey));

        for (Prop p : drawOrder) {
            PointF s = project(p.x, p.y, p.z);
            Bitmap frame = sprites.get(spriteKey(p.type, p.material))[directionIndex(p.yaw)];
            c.drawBitmap(frame, Math.round(s.x - FRAME * 0.5f), Math.round(s.y - FRAME * 0.5f), spritePaint);
        }

        for (Prop p : props) {
            if (p != grabbed && !p.frozen) continue;

            PointF s = project(p.x, p.y, p.z);
            int color = p == grabbed ? rgb("F4D35E") : rgb("70D6FF");

            paint.setStyle(Paint.Style.FILL);
            paint.setColor(rgb("09090B"));
            c.drawRect(s.x - 15, s.y - 15, s.x + 16, s.y - 13, paint);
            c.drawRect(s.x - 15, s.y + 14, s.x + 16, s.y + 16, paint);
            c.drawRect(s.x - 15, s.y - 15, s.x - 13, s.y + 16, paint);
            c.drawRect(s.x + 14, s.y - 15, s.x + 16, s.y + 16, paint);

            paint.setColor(color);
            c.drawRect(s.x - 14, s.y - 14, s.x - 7, s.y - 13, paint);
            c.drawRect(s.x + 7, s.y - 14, s.x + 14, s.y - 13, paint);
            c.drawRect(s.x - 14, s.y + 13, s.x - 7, s.y + 14, paint);
            c.drawRect(s.x + 7, s.y + 13, s.x + 14, s.y + 14, paint);
        }
    }

    private void drawHud(Canvas c) {
        panel(c, 7, H - 35, 38, 28, rgb("25303A"));
        panel(c, W - 66, H - 35, 59, 28, rgb("25303A"));
        panel(c, W - 61, 7, 54, 28, rgb("25303A"));

        if (mode == Mode.SPAWN) {
            panel(c, 0, H - 116, W, 116, rgb("15171D"));
            float cellW = W / 4f;
            float cellsTop = H - 94f;
            for (int i = 0; i < PropType.values().length; i++) {
                int row = i / 4;
                int col = i % 4;
                panel(c, col * cellW + 3, cellsTop + row * 23, cellW - 6, 21, (i & 1) == 0 ? rgb("2C211D") : rgb("22252A"));
            }
        }

        if (mode == Mode.SETTINGS) {
            panel(c, W - 158, 42, 151, 142, rgb("14181E"));
        }

        if (mode == Mode.CONTEXT) {
            float cx = clamp(contextX, 68f, W - 68f);
            float cy = clamp(contextY, 48f, H - 48f);
            panel(c, cx - 66, cy - 41, 132, 82, rgb("17151A"));
            paint.setColor(rgb("5A321F"));
            c.drawRect(cx - 1, cy - 40, cx + 1, cy + 40, paint);
            c.drawRect(cx - 65, cy - 1, cx + 65, cy + 1, paint);
        }

        pixelText(c, "+", 20, H - 16, Color.WHITE);
        pixelText(c, "UNDO", W - 60, H - 17, Color.WHITE);
        pixelText(c, "MENU", W - 55, 25, Color.WHITE);
        pixelText(c, "PIXEL PHYSICS 2.5D", 8, 12, rgb("F0C06A"));
        pixelText(c, "ANDROID CANVAS / ZERO NATIVE LIBS", 8, 24, rgb("8FA3AD"));

        if (grabbed != null) {
            pixelText(c, grabbed.type.label + " / " + (directionIndex(grabbed.yaw) * 90) + " DEG", 8, 37, Color.WHITE);
        }

        if (mode == Mode.SPAWN) drawSpawnDrawer(c);
        if (mode == Mode.SETTINGS) drawSettings(c);
        if (mode == Mode.CONTEXT) drawContext(c);
    }

    private void panel(Canvas c, float x, float y, float w, float h, int fill) {
        paint.setStyle(Paint.Style.FILL);
        paint.setColor(rgb("09090C"));
        c.drawRect(x - 2, y - 2, x + w + 2, y + h + 2, paint);
        paint.setColor(fill);
        c.drawRect(x, y, x + w, y + h, paint);
        paint.setColor(rgb("6A5A50"));
        c.drawRect(x, y, x + w, y + 1, paint);
        c.drawRect(x, y, x + 1, y + h, paint);
        paint.setColor(rgb("101218"));
        c.drawRect(x, y + h - 1, x + w, y + h, paint);
        c.drawRect(x + w - 1, y, x + w, y + h, paint);
    }

    private void drawSpawnDrawer(Canvas c) {
        pixelText(c, "SPAWN / 4 AUTHORED ANGLES", 8, H - 105, rgb("F0C06A"));

        int preview = (int) ((SystemClock.uptimeMillis() / 650L) % 4L);
        float cellW = W / 4f;
        float cellsTop = H - 94f;
        PropType[] values = PropType.values();

        for (int i = 0; i < values.length; i++) {
            int row = i / 4;
            int col = i % 4;
            float x = col * cellW + 5;
            float y = cellsTop + row * 23;

            Bitmap frame = sprites.get(spriteKey(values[i], values[i].defaultMaterial))[preview];
            Rect from = new Rect(0, 0, FRAME, FRAME);
            RectF to = new RectF(x, y - 1, x + 22, y + 21);
            c.drawBitmap(frame, from, to, spritePaint);

            pixelText(c, values[i].label, x + 27, y + 9, Color.WHITE);
            pixelText(c, (preview * 90) + "", x + 27, y + 18, rgb("8C9AA7"));
        }
    }

    private void drawSettings(Canvas c) {
        float x = W - 149;
        float y = 56;
        pixelText(c, "TRUE 2.5D PIXEL", x, y, rgb("F0C06A"));
        pixelText(c, "RESET WORLD", x, y + 32, Color.WHITE);
        pixelText(c, "HAPTICS: " + (haptics ? "ON" : "OFF"), x, y + 62, Color.WHITE);
        pixelText(c, "XYZ CUSTOM PHYSICS", x, y + 92, rgb("70D6FF"));
        pixelText(c, "NO JNI / NO ENGINE", x, y + 108, rgb("70D6FF"));
        pixelText(c, "MENU = CLOSE", x, y + 127, rgb("776F72"));
    }

    private void drawContext(Canvas c) {
        float cx = clamp(contextX, 68f, W - 68f);
        float cy = clamp(contextY, 48f, H - 48f);

        pixelText(c, contextProp != null && contextProp.frozen ? "UNFREEZE" : "FREEZE", cx - 58, cy - 21, rgb("70D6FF"));
        pixelText(c, "DELETE", cx + 10, cy - 21, rgb("FF806C"));
        pixelText(c, "DUPLICATE", cx - 58, cy + 23, Color.WHITE);
        pixelText(c, "MATERIAL", cx + 10, cy + 23, rgb("F0C06A"));
    }

    private void pixelText(Canvas c, String text, float x, float y, int color) {
        paint.setStyle(Paint.Style.FILL);
        paint.setTextSize(7f);
        paint.setTypeface(Typeface.create(Typeface.MONOSPACE, Typeface.BOLD));
        paint.setAntiAlias(false);

        paint.setColor(rgb("08080A"));
        c.drawText(text, x + 1, y + 1, paint);
        paint.setColor(color);
        c.drawText(text, x, y, paint);
    }

    private void drawWallRailY(Canvas c, float y, float z, float width, int color) {
        line(c, project(-ROOM, y, z), project(ROOM, y, z), width, color);
    }

    private void drawWallRailX(Canvas c, float x, float z, float width, int color) {
        line(c, project(x, -ROOM, z), project(x, ROOM, z), width, color);
    }

    private void drawVerticalPost(Canvas c, float x, float y, float z0, float z1, float width, int base, int highlight) {
        PointF a = project(x, y, z0);
        PointF b = project(x, y, z1);
        paint.setColor(base);
        c.drawRect(a.x - width * 0.5f, b.y, a.x + width * 0.5f, a.y, paint);
        paint.setColor(highlight);
        c.drawRect(a.x - width * 0.5f + 2, b.y, a.x - width * 0.5f + 4, a.y, paint);
    }

    private void drawWallRectY(Canvas c, float y, float x0, float x1, float z0, float z1, int color) {
        fillPoly(c, paint, color,
                project(x0, y, z0), project(x1, y, z0), project(x1, y, z1), project(x0, y, z1));
    }

    private void drawWallLineY(Canvas c, float y, float x0, float z0, float x1, float z1, float width, int color) {
        line(c, project(x0, y, z0), project(x1, y, z1), width, color);
    }

    private void drawShelfY(Canvas c, float y, float x0, float x1, float z) {
        line(c, project(x0, y - 0.05f, z), project(x1, y - 0.05f, z), 6f, rgb("B65D2D"));
        line(c, project(x0, y - 0.05f, z + 0.08f), project(x1, y - 0.05f, z + 0.08f), 2f, rgb("E28A49"));
    }

    private void drawShelfX(Canvas c, float x, float y0, float y1, float z) {
        line(c, project(x - 0.05f, y0, z), project(x - 0.05f, y1, z), 6f, rgb("B65D2D"));
        line(c, project(x - 0.05f, y0, z + 0.08f), project(x - 0.05f, y1, z + 0.08f), 2f, rgb("E28A49"));
    }

    private void line(Canvas c, PointF a, PointF b, float width, int color) {
        paint.setStyle(Paint.Style.STROKE);
        paint.setStrokeWidth(width);
        paint.setStrokeCap(Paint.Cap.SQUARE);
        paint.setAntiAlias(false);
        paint.setColor(color);
        c.drawLine(Math.round(a.x), Math.round(a.y), Math.round(b.x), Math.round(b.y), paint);
        paint.setStyle(Paint.Style.FILL);
        paint.setStrokeWidth(1f);
    }

    private Prop pick(float sx, float sy) {
        drawOrder.clear();
        drawOrder.addAll(props);
        drawOrder.sort(Comparator.comparingDouble(this::depthKey));

        for (int i = drawOrder.size() - 1; i >= 0; i--) {
            Prop p = drawOrder.get(i);
            PointF s = project(p.x, p.y, p.z);
            float hw = Math.max(14f, p.type.w * 10f);
            float hh = Math.max(14f, p.type.h * 13f);
            if (Math.abs(sx - s.x) <= hw && Math.abs(sy - s.y) <= hh) return p;
        }
        return null;
    }

    private void startGrab(Prop p, float sx, float sy) {
        if (p == null || p.frozen) return;

        grabbed = p;
        mode = Mode.GRAB;

        PointF world = unprojectAtZ(sx, sy, p.z);
        grabOffsetX = world.x - p.x;
        grabOffsetY = world.y - p.y;

        targetX = p.x;
        targetY = p.y;
        targetZ = p.z;
    }

    private void updateGrabTarget(float sx, float sy) {
        if (grabbed == null) return;

        PointF world = unprojectAtZ(sx, sy, targetZ);
        targetX = world.x - grabOffsetX;
        targetY = world.y - grabOffsetY;

        targetX = clamp(targetX, -ROOM + grabbed.radius(), ROOM - grabbed.radius());
        targetY = clamp(targetY, -ROOM + grabbed.radius(), ROOM - grabbed.radius());
    }

    private void endGrab() {
        grabbed = null;
        secondaryId = -1;
        if (mode == Mode.GRAB) mode = Mode.IDLE;
    }

    private void checkLongPress(long nowNanos) {
        if (pressed == null || primaryId < 0 || contextTriggered) return;
        if (mode != Mode.GRAB && mode != Mode.IDLE) return;

        float moved = distance(downX, downY, primaryX, primaryY);
        if (moved < 6f && nowNanos - downNanos > 550_000_000L) {
            contextTriggered = true;
            contextProp = pressed;
            contextX = downX;
            contextY = downY;

            if (mode == Mode.GRAB) endGrab();
            primaryId = -1;
            mode = Mode.CONTEXT;
            feedback();
        }
    }

    @Override
    public boolean onTouchEvent(MotionEvent event) {
        int action = event.getActionMasked();
        int index = event.getActionIndex();
        int pointerId = event.getPointerId(index);

        if (action == MotionEvent.ACTION_DOWN) {
            float x = toLogicalX(event.getX(index));
            float y = toLogicalY(event.getY(index));
            if (!insidePresentation(event.getX(index), event.getY(index))) return false;

            if (handleUiDown(x, y)) return true;

            primaryId = pointerId;
            secondaryId = -1;
            primaryX = downX = x;
            primaryY = downY = y;
            downNanos = System.nanoTime();
            contextTriggered = false;

            pressed = pick(x, y);
            if (pressed != null && !pressed.frozen) {
                startGrab(pressed, x, y);
            } else {
                mode = Mode.IDLE;
            }

            return true;
        }

        if (action == MotionEvent.ACTION_POINTER_DOWN) {
            if (mode == Mode.GRAB && secondaryId < 0) {
                secondaryId = pointerId;
                int pIndex = event.findPointerIndex(primaryId);
                if (pIndex >= 0) {
                    float x1 = toLogicalX(event.getX(pIndex));
                    float y1 = toLogicalY(event.getY(pIndex));
                    float x2 = toLogicalX(event.getX(index));
                    float y2 = toLogicalY(event.getY(index));
                    lastTwoDistance = distance(x1, y1, x2, y2);
                    lastTwoAngle = angle(x1, y1, x2, y2);
                }
            }
            return true;
        }

        if (action == MotionEvent.ACTION_MOVE) {
            int pIndex = event.findPointerIndex(primaryId);
            if (pIndex >= 0) {
                primaryX = toLogicalX(event.getX(pIndex));
                primaryY = toLogicalY(event.getY(pIndex));

                if (mode == Mode.GRAB && grabbed != null) {
                    updateGrabTarget(primaryX, primaryY);
                }
            }

            if (mode == Mode.GRAB && grabbed != null && secondaryId >= 0) {
                int sIndex = event.findPointerIndex(secondaryId);
                pIndex = event.findPointerIndex(primaryId);

                if (pIndex >= 0 && sIndex >= 0) {
                    float x1 = toLogicalX(event.getX(pIndex));
                    float y1 = toLogicalY(event.getY(pIndex));
                    float x2 = toLogicalX(event.getX(sIndex));
                    float y2 = toLogicalY(event.getY(sIndex));

                    float dist = distance(x1, y1, x2, y2);
                    float dd = dist - lastTwoDistance;
                    targetZ = clamp(targetZ + dd * 0.022f, grabbed.type.h * 0.5f, 5f);

                    float a = angle(x1, y1, x2, y2);
                    float da = wrapAngle(a - lastTwoAngle);
                    grabbed.yaw += da;
                    grabbed.spin += da * 2.5f;

                    lastTwoDistance = dist;
                    lastTwoAngle = a;
                    updateGrabTarget(primaryX, primaryY);
                }
            }
            return true;
        }

        if (action == MotionEvent.ACTION_POINTER_UP) {
            if (pointerId == secondaryId) {
                secondaryId = -1;
                return true;
            }
            if (pointerId == primaryId) {
                endGrab();
                pressed = null;
                primaryId = -1;
                return true;
            }
            return true;
        }

        if (action == MotionEvent.ACTION_UP || action == MotionEvent.ACTION_CANCEL) {
            if (mode == Mode.GRAB) endGrab();
            pressed = null;
            primaryId = -1;
            secondaryId = -1;
            return true;
        }

        return true;
    }

    private boolean handleUiDown(float x, float y) {
        if (mode == Mode.SPAWN) {
            float cellsTop = H - 94f;
            if (y >= cellsTop) {
                int col = clampInt((int) (x / (W / 4f)), 0, 3);
                int row = clampInt((int) ((y - cellsTop) / 23f), 0, 3);
                int idx = row * 4 + col;
                if (idx >= 0 && idx < PropType.values().length) {
                    spawnWithUndo(PropType.values()[idx]);
                }
            }
            mode = Mode.IDLE;
            return true;
        }

        if (mode == Mode.SETTINGS) {
            float left = W - 158f;
            if (x < left || y < 42f || y > 184f) {
                mode = Mode.IDLE;
                return true;
            }

            if (y >= 73f && y < 108f) {
                resetWorld();
                mode = Mode.IDLE;
                return true;
            }

            if (y >= 108f && y < 139f) {
                haptics = !haptics;
                feedback();
                saveWorld();
                return true;
            }

            return true;
        }

        if (mode == Mode.CONTEXT) {
            float cx = clamp(contextX, 68f, W - 68f);
            float cy = clamp(contextY, 48f, H - 48f);
            boolean left = x < cx;
            boolean top = y < cy;
            Prop target = contextProp;

            mode = Mode.IDLE;
            contextProp = null;

            if (target != null) {
                if (top && left) setFrozen(target, !target.frozen, true);
                else if (top) deleteWithUndo(target);
                else if (left) duplicateWithUndo(target);
                else cycleMaterial(target);
            }

            return true;
        }

        if (x < 52f && y > H - 42f) {
            mode = Mode.SPAWN;
            feedback();
            return true;
        }

        if (x > W - 73f && y > H - 42f) {
            doUndo();
            return true;
        }

        if (x > W - 70f && y < 42f) {
            mode = Mode.SETTINGS;
            feedback();
            return true;
        }

        return false;
    }

    private boolean insidePresentation(float sx, float sy) {
        return sx >= presentX && sx <= presentX + presentW && sy >= presentY && sy <= presentY + presentH;
    }

    private float toLogicalX(float sx) {
        return clamp((sx - presentX) / Math.max(0.0001f, presentScale), 0f, W - 1f);
    }

    private float toLogicalY(float sy) {
        return clamp((sy - presentY) / Math.max(0.0001f, presentScale), 0f, H - 1f);
    }

    private void feedback() {
        if (!haptics) return;
        performHapticFeedback(HapticFeedbackConstants.KEYBOARD_TAP);
    }

    private void saveWorld() {
        try {
            JSONArray array = new JSONArray();
            for (Prop p : props) {
                SaveState s = p.snapshot();
                JSONObject o = new JSONObject();
                o.put("id", s.id);
                o.put("type", s.type);
                o.put("material", s.material);
                o.put("x", s.x);
                o.put("y", s.y);
                o.put("z", s.z);
                o.put("yaw", s.yaw);
                o.put("frozen", s.frozen);
                array.put(o);
            }
            prefs.edit()
                    .putString("world", array.toString())
                    .putBoolean("haptics", haptics)
                    .apply();
        } catch (Exception ignored) {
        }
    }

    private boolean restoreWorld() {
        String data = prefs.getString("world", "");
        if (data == null || data.isEmpty()) return false;

        try {
            JSONArray array = new JSONArray(data);
            if (array.length() == 0) return false;

            for (int i = 0; i < array.length(); i++) {
                JSONObject o = array.getJSONObject(i);
                SaveState s = new SaveState();
                s.id = o.getInt("id");
                s.type = o.getString("type");
                s.material = o.getString("material");
                s.x = (float) o.getDouble("x");
                s.y = (float) o.getDouble("y");
                s.z = (float) o.getDouble("z");
                s.yaw = (float) o.getDouble("yaw");
                s.frozen = o.optBoolean("frozen", false);
                spawnFromState(s);
            }

            return !props.isEmpty();
        } catch (Exception ignored) {
            props.clear();
            byId.clear();
            nextId = 1;
            return false;
        }
    }

    private void spawnFromState(SaveState s) {
        spawnInternal(
                PropType.valueOf(s.type),
                MaterialKind.valueOf(s.material),
                s.x, s.y, s.z, s.yaw, s.frozen, s.id
        );
    }

    private void resetWorld() {
        endGrab();
        props.clear();
        byId.clear();
        undo.clear();
        nextId = 1;
        pressed = null;
        contextProp = null;
        mode = Mode.IDLE;
        prefs.edit().remove("world").apply();
        createStarterSet();
        feedback();
        saveWorld();
    }

    private static float distance(float x1, float y1, float x2, float y2) {
        float dx = x2 - x1;
        float dy = y2 - y1;
        return (float) Math.sqrt(dx * dx + dy * dy);
    }

    private static float angle(float x1, float y1, float x2, float y2) {
        return (float) Math.atan2(y2 - y1, x2 - x1);
    }

    private static float wrapAngle(float a) {
        while (a > Math.PI) a -= (float) (Math.PI * 2.0);
        while (a < -Math.PI) a += (float) (Math.PI * 2.0);
        return a;
    }

    private static float clamp(float v, float min, float max) {
        return Math.max(min, Math.min(max, v));
    }

    private static int clampInt(int v, int min, int max) {
        return Math.max(min, Math.min(max, v));
    }

    private static float lerp(float a, float b, float t) {
        return a + (b - a) * t;
    }
}
