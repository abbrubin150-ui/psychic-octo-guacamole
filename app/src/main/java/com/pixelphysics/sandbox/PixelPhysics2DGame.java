package com.pixelphysics.sandbox;

import com.badlogic.gdx.*;
import com.badlogic.gdx.graphics.*;
import com.badlogic.gdx.graphics.g2d.*;
import com.badlogic.gdx.graphics.glutils.*;
import com.badlogic.gdx.math.*;
import com.badlogic.gdx.physics.box2d.*;
import com.badlogic.gdx.utils.*;

import java.util.Locale;

public class PixelPhysics2DGame extends ApplicationAdapter implements InputProcessor {
    private static final int W = 480;
    private static final int H = 270;
    private static final float PPM = 24f;
    private static final float STEP = 1f / 60f;
    private static final float MAX_ACCUM = 0.12f;
    private static final int FRAME_SIZE = 96;
    private static final int MAX_PROPS = 100;

    private SpriteBatch batch;
    private ShapeRenderer shapes;
    private BitmapFont font;
    private FrameBuffer buffer;
    private TextureRegion bufferRegion;
    private final Matrix4 projection = new Matrix4();

    private float presentScale = 1f;
    private float presentX, presentY, presentW = W, presentH = H;

    private World world;
    private Body groundBody;
    private final Array<StaticSurface> surfaces = new Array<>();
    private final Array<Prop> props = new Array<>();
    private final ObjectMap<Integer, Prop> byId = new ObjectMap<>();
    private final ObjectMap<String, Texture[]> spriteSets = new ObjectMap<>();
    private final Array<UndoAction> undo = new Array<>();
    private int nextId = 1;

    private Preferences prefs;
    private Json json;
    private float autosaveClock;
    private float accumulator;
    private boolean haptics = true;

    private enum Mode { IDLE, GRAB, CONTEXT, SPAWN, SETTINGS }
    private Mode mode = Mode.IDLE;

    private final int[] px = new int[10];
    private final int[] py = new int[10];
    private int primaryPointer = -1;
    private int secondPointer = -1;
    private float downX, downY;
    private long downNanos;
    private float lastTwoAngle;
    private float accumulatedTorque;
    private Prop pressed;
    private Prop grabbed;
    private final Vector2 localGrab = new Vector2();
    private final Vector2 grabTarget = new Vector2();
    private boolean contextTriggered;
    private float contextX, contextY;
    private Prop contextProp;

    private final Vector2 tmpA = new Vector2();
    private final Vector2 tmpB = new Vector2();

    private enum MaterialKind { MAHOGANY, METAL, RUBBER }

    private enum PropType {
        CUBE("CUBE", 28, 28, false, 1.0f, MaterialKind.MAHOGANY),
        BEAM_SHORT("BEAM S", 58, 14, false, 1.3f, MaterialKind.MAHOGANY),
        BEAM_LONG("BEAM L", 88, 14, false, 2.1f, MaterialKind.MAHOGANY),
        PLANK("PLANK", 70, 10, false, 1.2f, MaterialKind.MAHOGANY),
        WOOD_BALL("WOOD BALL", 25, 25, true, 0.7f, MaterialKind.MAHOGANY),
        METAL_BALL("METAL BALL", 25, 25, true, 4.0f, MaterialKind.METAL),
        WEIGHT("WEIGHT", 31, 31, false, 9.0f, MaterialKind.METAL),
        WHEEL("WHEEL", 31, 31, true, 1.1f, MaterialKind.MAHOGANY),
        RUBBER_BALL("RUBBER", 27, 27, true, 0.8f, MaterialKind.RUBBER),
        BARREL("BARREL", 28, 40, false, 3.2f, MaterialKind.METAL),
        CRATE("CRATE", 35, 35, false, 1.8f, MaterialKind.MAHOGANY),
        RAMP("RAMP", 60, 30, false, 1.8f, MaterialKind.MAHOGANY);

        final String label;
        final int wPx, hPx;
        final boolean circle;
        final float mass;
        final MaterialKind defaultMaterial;

        PropType(String label, int wPx, int hPx, boolean circle, float mass, MaterialKind defaultMaterial) {
            this.label = label;
            this.wPx = wPx;
            this.hPx = hPx;
            this.circle = circle;
            this.mass = mass;
            this.defaultMaterial = defaultMaterial;
        }
    }

    private static class Palette {
        final Color outline, dark, base, mid, light, hi;
        Palette(String outline, String dark, String base, String mid, String light, String hi) {
            this.outline = Color.valueOf(outline);
            this.dark = Color.valueOf(dark);
            this.base = Color.valueOf(base);
            this.mid = Color.valueOf(mid);
            this.light = Color.valueOf(light);
            this.hi = Color.valueOf(hi);
        }
    }

    private static class StaticSurface {
        float x, y, w, h;
        String kind;
        Body body;
    }

    private class Prop {
        int id;
        PropType type;
        MaterialKind material;
        Body body;
        boolean frozen;

        SaveState save() {
            SaveState s = new SaveState();
            s.id = id;
            s.type = type.name();
            s.material = material.name();
            s.x = body.getPosition().x;
            s.y = body.getPosition().y;
            s.angle = body.getAngle();
            s.frozen = frozen;
            return s;
        }
    }

    private static class SaveState {
        public int id;
        public String type;
        public String material;
        public float x, y, angle;
        public boolean frozen;
        public SaveState() {}
    }

    private interface UndoAction { void undo(); }

    @Override
    public void create() {
        Locale.setDefault(Locale.US);
        prefs = Gdx.app.getPreferences("pixel-physics-2d-world-v1");
        json = new Json();

        batch = new SpriteBatch();
        shapes = new ShapeRenderer();
        font = new BitmapFont();
        font.getData().setScale(0.72f);
        font.getRegion().getTexture().setFilter(Texture.TextureFilter.Nearest, Texture.TextureFilter.Nearest);
        projection.setToOrtho2D(0, 0, W, H);

        world = new World(new Vector2(0, -9.81f), true);
        BodyDef gd = new BodyDef();
        gd.type = BodyDef.BodyType.StaticBody;
        groundBody = world.createBody(gd);

        createAllSpriteSets();
        buildWorkshopPhysics();
        recreateBuffer();
        if (!restoreWorld()) createStarterSet();

        Gdx.input.setInputProcessor(this);
    }

    private Palette palette(MaterialKind k) {
        switch (k) {
            case METAL: return new Palette("11151A","283039","4B5964","71808A","A7B1B7","E0E4E6");
            case RUBBER: return new Palette("101713","1B271E","2D3C30","405442","607761","9AAE94");
            default: return new Palette("24100D","421712","6D2619","963A21","C35D30","E78745");
        }
    }

    private void createAllSpriteSets() {
        for (PropType t : PropType.values()) {
            for (MaterialKind m : MaterialKind.values()) {
                Texture[] frames = new Texture[4];
                for (int d = 0; d < 4; d++) frames[d] = makeSprite(t, m, d);
                spriteSets.put(spriteKey(t, m), frames);
            }
        }
    }

    private String spriteKey(PropType t, MaterialKind m) {
        return t.name() + ":" + m.name();
    }

    private Texture makeSprite(PropType type, MaterialKind material, int dir) {
        Pixmap p = new Pixmap(FRAME_SIZE, FRAME_SIZE, Pixmap.Format.RGBA8888);
        p.setBlending(Pixmap.Blending.None);
        p.setColor(0,0,0,0);
        p.fill();
        p.setBlending(Pixmap.Blending.SourceOver);
        Palette pal = palette(material);

        if (type == PropType.WOOD_BALL || type == PropType.METAL_BALL || type == PropType.RUBBER_BALL) {
            drawBall(p, type, pal, dir);
        } else if (type == PropType.WHEEL) {
            drawWheel(p, pal, dir);
        } else if (type == PropType.CRATE) {
            drawCrate(p, pal, dir);
        } else if (type == PropType.BARREL) {
            drawBarrel(p, pal, dir);
        } else if (type == PropType.RAMP) {
            drawRamp(p, pal, dir);
        } else if (type == PropType.WEIGHT) {
            drawWeight(p, pal, dir);
        } else {
            drawRectProp(p, type, pal, dir);
        }

        Texture tex = new Texture(p);
        p.dispose();
        tex.setFilter(Texture.TextureFilter.Nearest, Texture.TextureFilter.Nearest);
        return tex;
    }

    private void drawRectProp(Pixmap p, PropType type, Palette pal, int dir) {
        int w = (dir % 2 == 0) ? type.wPx : type.hPx;
        int h = (dir % 2 == 0) ? type.hPx : type.wPx;
        int x = (FRAME_SIZE - w) / 2, y = (FRAME_SIZE - h) / 2;

        p.setColor(pal.outline); p.fillRectangle(x-1,y-1,w+2,h+2);
        p.setColor(pal.base); p.fillRectangle(x,y,w,h);

        if (dir == 0) {
            p.setColor(pal.light); p.fillRectangle(x+1,y+1,w-2,2);
            p.setColor(pal.dark); p.fillRectangle(x+1,y+h-3,w-2,2);
        } else if (dir == 1) {
            p.setColor(pal.light); p.fillRectangle(x+w-3,y+1,2,h-2);
            p.setColor(pal.dark); p.fillRectangle(x+1,y+1,2,h-2);
        } else if (dir == 2) {
            p.setColor(pal.light); p.fillRectangle(x+1,y+h-3,w-2,2);
            p.setColor(pal.dark); p.fillRectangle(x+1,y+1,w-2,2);
        } else {
            p.setColor(pal.light); p.fillRectangle(x+1,y+1,2,h-2);
            p.setColor(pal.dark); p.fillRectangle(x+w-3,y+1,2,h-2);
        }

        p.setColor(pal.mid);
        if (w >= h) {
            for (int yy=y+4; yy<y+h-3; yy+=4) {
                int start = x + 4 + ((yy + dir*3) % 7);
                p.drawLine(start,yy,Math.min(x+w-4,start+Math.max(4,w/3)),yy);
            }
        } else {
            for (int xx=x+4; xx<x+w-3; xx+=4) {
                int start = y + 4 + ((xx + dir*5) % 7);
                p.drawLine(xx,start,xx,Math.min(y+h-4,start+Math.max(4,h/3)));
            }
        }

        p.setColor(pal.hi);
        int hx = dir==1 || dir==2 ? x+w-5 : x+3;
        int hy = dir>=2 ? y+h-5 : y+3;
        p.fillRectangle(hx,hy,2,2);
    }

    private void drawBall(Pixmap p, PropType type, Palette pal, int dir) {
        int r = type.wPx/2;
        int cx = FRAME_SIZE/2, cy = FRAME_SIZE/2;
        p.setColor(pal.outline); p.fillCircle(cx,cy,r+1);
        p.setColor(pal.dark); p.fillCircle(cx,cy,r);
        p.setColor(pal.base); p.fillCircle(cx,cy,r-2);
        p.setColor(pal.mid); p.fillCircle(cx,cy,Math.max(3,r-5));

        int[][] offsets={{-4,4},{4,4},{4,-4},{-4,-4}};
        int ox=offsets[dir][0], oy=offsets[dir][1];
        p.setColor(pal.light); p.fillRectangle(cx+ox-2,cy+oy-1,5,3);
        p.setColor(pal.hi); p.fillRectangle(cx+ox-1,cy+oy,2,1);
        p.setColor(pal.outline);
        if (dir%2==0) p.drawLine(cx-r+4,cy,cx+r-4,cy);
        else p.drawLine(cx,cy-r+4,cx,cy+r-4);
    }

    private void drawWheel(Pixmap p, Palette pal, int dir) {
        int cx=FRAME_SIZE/2, cy=FRAME_SIZE/2, r=15;
        p.setColor(pal.outline); p.fillCircle(cx,cy,r+1);
        p.setColor(pal.dark); p.fillCircle(cx,cy,r);
        p.setColor(pal.base); p.fillCircle(cx,cy,r-4);
        p.setColor(pal.outline); p.fillCircle(cx,cy,4);
        p.setColor(pal.hi); p.fillCircle(cx,cy,2);

        p.setColor(pal.light);
        if (dir%2==0) {
            p.fillRectangle(cx-r+5,cy-1,(r-5)*2,3);
            p.fillRectangle(cx-1,cy-r+5,3,(r-5)*2);
        } else {
            for(int i=-8;i<=8;i++){
                p.drawPixel(cx+i,cy+i);
                p.drawPixel(cx+i,cy-i);
                if(i%3==0){ p.drawPixel(cx+i,cy+i+1); p.drawPixel(cx+i,cy-i-1); }
            }
        }
        p.setColor(pal.hi);
        int[][] h={{-8,8},{8,8},{8,-8},{-8,-8}};
        p.fillRectangle(cx+h[dir][0]-1,cy+h[dir][1]-1,3,3);
    }

    private void drawCrate(Pixmap p, Palette pal, int dir) {
        int s=35,x=(FRAME_SIZE-s)/2,y=(FRAME_SIZE-s)/2;
        p.setColor(pal.outline); p.fillRectangle(x-1,y-1,s+2,s+2);
        p.setColor(pal.base); p.fillRectangle(x,y,s,s);
        p.setColor(pal.dark); p.fillRectangle(x+3,y+3,s-6,3); p.fillRectangle(x+3,y+s-6,s-6,3);
        p.fillRectangle(x+3,y+3,3,s-6); p.fillRectangle(x+s-6,y+3,3,s-6);
        p.setColor(pal.mid);
        if(dir==0||dir==2){
            p.drawLine(x+7,y+7,x+s-8,y+s-8);
            p.drawLine(x+7,y+s-8,x+s-8,y+7);
        }else{
            p.fillRectangle(x+s/2-2,y+7,4,s-14);
            p.fillRectangle(x+7,y+s/2-2,s-14,4);
        }
        p.setColor(pal.hi);
        int[][] q={{x+6,y+6},{x+s-8,y+6},{x+s-8,y+s-8},{x+6,y+s-8}};
        p.fillRectangle(q[dir][0],q[dir][1],2,2);
    }

    private void drawBarrel(Pixmap p, Palette pal, int dir) {
        int w=(dir%2==0)?28:40, h=(dir%2==0)?40:28;
        int x=(FRAME_SIZE-w)/2,y=(FRAME_SIZE-h)/2;
        p.setColor(pal.outline); p.fillRectangle(x-1,y+3,w+2,h-6);
        p.fillRectangle(x+2,y-1,w-4,h+2);
        p.setColor(pal.base); p.fillRectangle(x,y+3,w,h-6); p.fillRectangle(x+3,y,w-6,h);
        p.setColor(pal.mid);
        if(dir%2==0){
            p.fillRectangle(x+1,y+8,w-2,3); p.fillRectangle(x+1,y+h-11,w-2,3);
        }else{
            p.fillRectangle(x+8,y+1,3,h-2); p.fillRectangle(x+w-11,y+1,3,h-2);
        }
        p.setColor(pal.light);
        if(dir==0) p.fillRectangle(x+4,y+5,3,h-10);
        else if(dir==1) p.fillRectangle(x+5,y+h-7,w-10,3);
        else if(dir==2) p.fillRectangle(x+w-7,y+5,3,h-10);
        else p.fillRectangle(x+5,y+4,w-10,3);
    }

    private void drawRamp(Pixmap p, Palette pal, int dir) {
        int cx=FRAME_SIZE/2,cy=FRAME_SIZE/2, hw=30,hh=15;
        p.setColor(pal.outline);
        if(dir==0) p.fillTriangle(cx-hw-2,cy-hh-2,cx+hw+2,cy-hh-2,cx+hw+2,cy+hh+2);
        if(dir==1) p.fillTriangle(cx-hh-2,cy-hw-2,cx+hh+2,cy-hw-2,cx-hh-2,cy+hw+2);
        if(dir==2) p.fillTriangle(cx-hw-2,cy-hh-2,cx+hw+2,cy-hh-2,cx-hw-2,cy+hh+2);
        if(dir==3) p.fillTriangle(cx-hh-2,cy+hw+2,cx+hh+2,cy+hw+2,cx+hh+2,cy-hw-2);
        p.setColor(pal.base);
        if(dir==0) p.fillTriangle(cx-hw,cy-hh,cx+hw,cy-hh,cx+hw,cy+hh);
        if(dir==1) p.fillTriangle(cx-hh,cy-hw,cx+hh,cy-hw,cx-hh,cy+hw);
        if(dir==2) p.fillTriangle(cx-hw,cy-hh,cx+hw,cy-hh,cx-hw,cy+hh);
        if(dir==3) p.fillTriangle(cx-hh,cy+hw,cx+hh,cy+hw,cx+hh,cy-hw);
        p.setColor(pal.light);
        if(dir==0) p.drawLine(cx-hw+4,cy-hh+3,cx+hw-4,cy+hh-3);
        if(dir==1) p.drawLine(cx-hh+3,cy+hw-4,cx+hh-3,cy-hw+4);
        if(dir==2) p.drawLine(cx-hw+4,cy+hh-3,cx+hw-4,cy-hh+3);
        if(dir==3) p.drawLine(cx-hh+3,cy-hw+4,cx+hh-3,cy+hw-4);
    }

    private void drawWeight(Pixmap p, Palette pal, int dir) {
        int s=31,x=(FRAME_SIZE-s)/2,y=(FRAME_SIZE-s)/2;
        p.setColor(pal.outline); p.fillRectangle(x-1,y-1,s+2,s+2);
        p.setColor(pal.dark); p.fillRectangle(x,y,s,s);
        p.setColor(pal.base); p.fillRectangle(x+3,y+3,s-6,s-6);
        p.setColor(pal.mid); p.fillRectangle(x+7,y+7,s-14,s-14);
        p.setColor(pal.light);
        if(dir==0) p.fillRectangle(x+5,y+s-8,s-10,3);
        if(dir==1) p.fillRectangle(x+5,y+5,3,s-10);
        if(dir==2) p.fillRectangle(x+5,y+5,s-10,3);
        if(dir==3) p.fillRectangle(x+s-8,y+5,3,s-10);
        p.setColor(pal.hi);
        int[][] q={{x+7,y+s-10},{x+7,y+7},{x+s-10,y+7},{x+s-10,y+s-10}};
        p.fillRectangle(q[dir][0],q[dir][1],3,3);
    }

    private void buildWorkshopPhysics() {
        createStatic(0,0,W,14,"floor");
        createStatic(38,86,150,10,"wood");
        createStatic(294,72,112,10,"wood");
        createStatic(300,164,132,8,"metal");
        createStatic(206,48,74,8,"wood");
    }

    private void createStatic(float x,float y,float w,float h,String kind) {
        BodyDef bd=new BodyDef();
        bd.type=BodyDef.BodyType.StaticBody;
        bd.position.set((x+w/2f)/PPM,(y+h/2f)/PPM);
        Body b=world.createBody(bd);
        PolygonShape ps=new PolygonShape();
        ps.setAsBox(w/(2f*PPM),h/(2f*PPM));
        FixtureDef fd=new FixtureDef();
        fd.shape=ps; fd.friction=0.82f; fd.restitution=0.05f;
        b.createFixture(fd); ps.dispose();

        StaticSurface s=new StaticSurface();
        s.x=x;s.y=y;s.w=w;s.h=h;s.kind=kind;s.body=b;
        surfaces.add(s);
    }

    private void createStarterSet() {
        spawnInternal(PropType.CUBE, MaterialKind.MAHOGANY, 118, 132, 0f, false, nextId++);
        spawnInternal(PropType.BEAM_SHORT, MaterialKind.MAHOGANY, 208, 122, 0.18f, false, nextId++);
        spawnInternal(PropType.RUBBER_BALL, MaterialKind.RUBBER, 270, 155, 0f, false, nextId++);
        spawnInternal(PropType.CRATE, MaterialKind.MAHOGANY, 350, 118, 0f, false, nextId++);
        spawnInternal(PropType.METAL_BALL, MaterialKind.METAL, 320, 210, 0f, false, nextId++);
        saveWorld();
    }

    private Prop spawnInternal(PropType type, MaterialKind material, float pxCenter, float pyCenter,
                               float angle, boolean frozen, int requestedId) {
        if(props.size>=MAX_PROPS) return null;
        Prop p=new Prop();
        p.id=requestedId>0?requestedId:nextId++;
        nextId=Math.max(nextId,p.id+1);
        p.type=type;p.material=material;p.frozen=frozen;

        BodyDef bd=new BodyDef();
        bd.type=frozen?BodyDef.BodyType.StaticBody:BodyDef.BodyType.DynamicBody;
        bd.position.set(pxCenter/PPM,pyCenter/PPM);
        bd.angle=angle;
        bd.angularDamping=0.45f;
        bd.linearDamping=0.08f;
        p.body=world.createBody(bd);
        p.body.setUserData(p.id);

        Shape shape;
        float area;
        if(type.circle){
            CircleShape cs=new CircleShape();
            float r=type.wPx/(2f*PPM);
            cs.setRadius(r);
            shape=cs;
            area=MathUtils.PI*r*r;
        } else {
            PolygonShape ps=new PolygonShape();
            float hw=type.wPx/(2f*PPM),hh=type.hPx/(2f*PPM);
            ps.setAsBox(hw,hh);
            shape=ps;
            area=(type.wPx/PPM)*(type.hPx/PPM);
        }
        FixtureDef fd=new FixtureDef();
        fd.shape=shape;
        fd.density=Math.max(0.15f,type.mass/Math.max(0.05f,area));
        fd.friction=material==MaterialKind.RUBBER?0.9f:(material==MaterialKind.METAL?0.55f:0.78f);
        fd.restitution=material==MaterialKind.RUBBER?0.68f:(material==MaterialKind.METAL?0.08f:0.10f);
        p.body.createFixture(fd);
        shape.dispose();

        props.add(p);byId.put(p.id,p);
        return p;
    }

    private void removeProp(Prop p) {
        if(p==null)return;
        if(grabbed==p) endGrab();
        props.removeValue(p,true);
        byId.remove(p.id);
        world.destroyBody(p.body);
    }

    private void spawnWithUndo(PropType type) {
        float x=235+((props.size*23)%100)-50;
        float y=206+((props.size*11)%28);
        Prop p=spawnInternal(type,type.defaultMaterial,x,y,0f,false,nextId++);
        if(p==null)return;
        final int id=p.id;
        pushUndo(()->{Prop q=byId.get(id);if(q!=null)removeProp(q);});
        feedback();
        saveWorld();
    }

    private void deleteWithUndo(Prop p) {
        if(p==null)return;
        SaveState s=p.save();
        pushUndo(()->spawnFromState(s));
        removeProp(p);
        feedback();
        saveWorld();
    }

    private void duplicateWithUndo(Prop source) {
        if(source==null)return;
        Vector2 pos=source.body.getPosition();
        Prop p=spawnInternal(source.type,source.material,pos.x*PPM+16,pos.y*PPM+16,
                source.body.getAngle(),source.frozen,nextId++);
        if(p==null)return;
        final int id=p.id;
        pushUndo(()->{Prop q=byId.get(id);if(q!=null)removeProp(q);});
        feedback();saveWorld();
    }

    private void setFrozen(Prop p,boolean frozen,boolean record) {
        if(p==null||p.frozen==frozen)return;
        final int id=p.id;
        final boolean prior=p.frozen;
        if(record)pushUndo(()->{Prop q=byId.get(id);if(q!=null)setFrozen(q,prior,false);});
        p.frozen=frozen;
        p.body.setLinearVelocity(0,0);
        p.body.setAngularVelocity(0);
        p.body.setType(frozen?BodyDef.BodyType.StaticBody:BodyDef.BodyType.DynamicBody);
        feedback();saveWorld();
    }

    private void cycleMaterial(Prop p) {
        if(p==null)return;
        final int id=p.id;
        final MaterialKind prior=p.material;
        MaterialKind next=prior==MaterialKind.MAHOGANY?MaterialKind.METAL:
                (prior==MaterialKind.METAL?MaterialKind.RUBBER:MaterialKind.MAHOGANY);
        pushUndo(()->{Prop q=byId.get(id);if(q!=null)applyMaterial(q,prior);});
        applyMaterial(p,next);feedback();saveWorld();
    }

    private void applyMaterial(Prop p,MaterialKind k) {
        if(p==null)return;
        p.material=k;
        for(Fixture f:p.body.getFixtureList()){
            f.setFriction(k==MaterialKind.RUBBER?0.9f:(k==MaterialKind.METAL?0.55f:0.78f));
            f.setRestitution(k==MaterialKind.RUBBER?0.68f:(k==MaterialKind.METAL?0.08f:0.10f));
        }
    }

    private void pushUndo(UndoAction a){
        undo.add(a);
        while(undo.size>32)undo.removeIndex(0);
    }

    private void doUndo(){
        if(undo.size==0)return;
        undo.pop().undo();
        feedback();saveWorld();
    }

    private Prop pick(float wx,float wy) {
        final Prop[] hit={null};
        world.QueryAABB(fixture->{
            Body b=fixture.getBody();
            if(b.getType()!=BodyDef.BodyType.StaticBody || byId.containsKey((Integer)b.getUserData())){
                if(fixture.testPoint(wx,wy)){
                    Object data=b.getUserData();
                    if(data instanceof Integer) hit[0]=byId.get((Integer)data);
                    return false;
                }
            }
            return true;
        },wx-0.02f,wy-0.02f,wx+0.02f,wy+0.02f);
        return hit[0];
    }

    private void startGrab(Prop p,float wx,float wy) {
        if(p==null||p.frozen)return;
        grabbed=p;
        mode=Mode.GRAB;
        localGrab.set(p.body.getLocalPoint(tmpA.set(wx,wy)));
        grabTarget.set(wx,wy);
        p.body.setAwake(true);
    }

    private void endGrab() {
        grabbed=null;
        accumulatedTorque=0;
        secondPointer=-1;
        if(mode==Mode.GRAB)mode=Mode.IDLE;
    }

    private void updateGrabPhysics() {
        if(grabbed==null||mode!=Mode.GRAB||grabbed.frozen)return;
        Vector2 current=grabbed.body.getWorldPoint(localGrab);
        Vector2 velocity=grabbed.body.getLinearVelocityFromWorldPoint(current);
        float mass=Math.max(0.2f,grabbed.body.getMass());
        float massScale=(float)Math.pow(mass,0.34);
        float kp=42f*massScale;
        float kd=7.5f*(float)Math.sqrt(massScale);
        float maxF=75f*(float)Math.pow(Math.max(1f,mass),0.58);

        tmpA.set(grabTarget).sub(current).scl(kp).mulAdd(velocity,-kd);
        if(tmpA.len2()>maxF*maxF)tmpA.nor().scl(maxF);
        grabbed.body.applyForce(tmpA,current,true);

        if(Math.abs(accumulatedTorque)>0.0001f){
            grabbed.body.applyTorque(MathUtils.clamp(accumulatedTorque*4f,-22f,22f),true);
            accumulatedTorque*=0.55f;
        }
    }

    private void recreateBuffer(){
        if(buffer!=null)buffer.dispose();
        buffer=new FrameBuffer(Pixmap.Format.RGBA8888,W,H,false);
        buffer.getColorBufferTexture().setFilter(Texture.TextureFilter.Nearest,Texture.TextureFilter.Nearest);
        bufferRegion=new TextureRegion(buffer.getColorBufferTexture());
        bufferRegion.flip(false,true);
        updatePresentation();
    }

    private void updatePresentation(){
        float raw=Math.min(Gdx.graphics.getWidth()/(float)W,Gdx.graphics.getHeight()/(float)H);
        float integer=(float)Math.floor(raw);
        presentScale=integer>=1f?integer:raw;
        presentW=W*presentScale;presentH=H*presentScale;
        presentX=(Gdx.graphics.getWidth()-presentW)*0.5f;
        presentY=(Gdx.graphics.getHeight()-presentH)*0.5f;
    }

    private boolean insidePresentation(int sx,int sy){
        return sx>=presentX&&sx<=presentX+presentW&&sy>=presentY&&sy<=presentY+presentH;
    }

    private int vx(int sx){return MathUtils.clamp(Math.round((sx-presentX)/presentScale),0,W-1);}
    private int vy(int sy){return MathUtils.clamp(Math.round((sy-presentY)/presentScale),0,H-1);}
    private float worldX(int vx){return vx/PPM;}
    private float worldY(int vy){return (H-vy)/PPM;}

    @Override
    public void render(){
        float dt=Math.min(Gdx.graphics.getDeltaTime(),0.05f);

        if(pressed!=null&&primaryPointer>=0&&!contextTriggered&&
                (mode==Mode.GRAB||mode==Mode.IDLE)){
            float moved=Vector2.dst(downX,downY,px[primaryPointer],py[primaryPointer]);
            if(moved<6f&&(TimeUtils.nanoTime()-downNanos)>550_000_000L){
                contextTriggered=true;
                contextProp=pressed;
                contextX=downX;
                contextY=H-downY;
                if(mode==Mode.GRAB)endGrab();
                primaryPointer=-1;
                mode=Mode.CONTEXT;
            }
        }

        accumulator=Math.min(MAX_ACCUM,accumulator+dt);
        while(accumulator>=STEP){
            updateGrabPhysics();
            world.step(STEP,6,2);
            accumulator-=STEP;
        }

        buffer.begin();
        Gdx.gl.glViewport(0,0,W,H);
        Gdx.gl.glClearColor(0.055f,0.043f,0.050f,1f);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT);
        drawWorkshop();
        drawProps();
        drawOverlay();
        buffer.end();

        updatePresentation();
        Gdx.gl.glViewport(0,0,Gdx.graphics.getWidth(),Gdx.graphics.getHeight());
        Gdx.gl.glClearColor(0.012f,0.010f,0.014f,1f);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT);
        batch.setProjectionMatrix(new Matrix4().setToOrtho2D(0,0,Gdx.graphics.getWidth(),Gdx.graphics.getHeight()));
        batch.begin();
        batch.draw(bufferRegion,presentX,presentY,presentW,presentH);
        batch.end();

        autosaveClock+=dt;
        if(autosaveClock>3f){autosaveClock=0;saveWorld();}
    }

    private void drawWorkshop(){
        shapes.setProjectionMatrix(projection);
        shapes.begin(ShapeRenderer.ShapeType.Filled);

        shapes.setColor(Color.valueOf("171319"));shapes.rect(0,0,W,H);
        shapes.setColor(Color.valueOf("38231E"));shapes.rect(0,14,W,H-14);

        shapes.setColor(Color.valueOf("4A2B22"));
        for(int y=18;y<H;y+=18)shapes.rect(0,y,W,2);
        shapes.setColor(Color.valueOf("2B1A18"));
        for(int x=16;x<W;x+=32)shapes.rect(x,14,1,H-14);

        // Window.
        shapes.setColor(Color.valueOf("11151A"));shapes.rect(332,188,118,67);
        shapes.setColor(Color.valueOf("76A5B7"));shapes.rect(336,192,110,59);
        shapes.setColor(Color.valueOf("A7CDD2"));shapes.rect(336,225,110,26);
        shapes.setColor(Color.valueOf("6F8D58"));
        shapes.rect(336,192,110,19);shapes.rect(350,207,24,11);shapes.rect(405,204,35,16);
        shapes.setColor(Color.valueOf("26303A"));shapes.rect(389,192,4,59);shapes.rect(336,220,110,4);

        // Pegboard and simple tools.
        shapes.setColor(Color.valueOf("72503A"));shapes.rect(18,176,120,70);
        shapes.setColor(Color.valueOf("3A281F"));
        for(int y=182;y<240;y+=8)for(int x=24;x<132;x+=8)shapes.rect(x,y,2,2);
        shapes.setColor(Color.valueOf("87949C"));
        shapes.rect(38,196,5,31);shapes.rect(34,221,13,5);
        shapes.rect(70,192,4,38);shapes.rect(66,190,12,5);
        shapes.setColor(Color.valueOf("B94E2A"));shapes.rect(100,200,8,30);shapes.rect(96,224,16,5);

        // Poster.
        shapes.setColor(Color.valueOf("18212A"));shapes.rect(156,196,104,48);
        shapes.setColor(Color.valueOf("E0B15A"));shapes.rect(162,202,92,4);
        shapes.setColor(Color.valueOf("8C776B"));shapes.rect(162,213,70,3);shapes.rect(162,221,80,3);shapes.rect(162,229,52,3);

        for(StaticSurface s:surfaces)drawStaticSurface(s);

        // Drop zone.
        shapes.setColor(Color.valueOf("101318"));shapes.rect(410,14,60,44);
        for(int x=412;x<468;x+=8){
            shapes.setColor(((x/8)&1)==0?Color.valueOf("E2B33F"):Color.valueOf("25262B"));
            shapes.rect(x,54,8,4);
        }
        shapes.end();

        batch.setProjectionMatrix(projection);
        batch.begin();
        pixelText("PIXEL PHYSICS / 2D",8,264,Color.valueOf("F0C06A"));
        pixelText("4-DIRECTION SPRITES",8,251,Color.valueOf("8FA3AD"));
        pixelText("DROP",426,45,Color.valueOf("E6C79C"));
        pixelText("ZONE",424,34,Color.valueOf("E6C79C"));
        batch.end();
    }

    private void drawStaticSurface(StaticSurface s){
        Color base=s.kind.equals("metal")?Color.valueOf("38444E"):Color.valueOf("6A3823");
        Color dark=s.kind.equals("metal")?Color.valueOf("1E252B"):Color.valueOf("351812");
        Color light=s.kind.equals("metal")?Color.valueOf("788995"):Color.valueOf("B45B31");
        shapes.setColor(dark);shapes.rect(s.x-1,s.y-1,s.w+2,s.h+2);
        shapes.setColor(base);shapes.rect(s.x,s.y,s.w,s.h);
        shapes.setColor(light);shapes.rect(s.x+1,s.y+s.h-2,s.w-2,1);
        if(!s.kind.equals("metal")){
            shapes.setColor(Color.valueOf("4A2318"));
            for(float x=s.x+9;x<s.x+s.w-4;x+=17)shapes.rect(x,s.y+3,8,1);
        }
    }

    private int directionIndex(float radians){
        int q=Math.round(radians*MathUtils.radiansToDegrees/90f);
        q%=4;if(q<0)q+=4;
        return q;
    }

    private void drawProps(){
        batch.setProjectionMatrix(projection);
        batch.begin();
        for(Prop p:props){
            Vector2 pos=p.body.getPosition();
            Texture[] frames=spriteSets.get(spriteKey(p.type,p.material));
            Texture tex=frames[directionIndex(p.body.getAngle())];
            float x=Math.round(pos.x*PPM)-FRAME_SIZE/2f;
            float y=Math.round(pos.y*PPM)-FRAME_SIZE/2f;
            batch.draw(tex,x,y);
        }
        batch.end();
    }

    private void drawOverlay(){
        shapes.setProjectionMatrix(projection);
        shapes.begin(ShapeRenderer.ShapeType.Filled);

        for(Prop p:props){
            if(p!=grabbed&&!p.frozen)continue;
            Vector2 pos=p.body.getPosition();
            float x=Math.round(pos.x*PPM),y=Math.round(pos.y*PPM);
            Color c=p==grabbed?Color.valueOf("F4D35E"):Color.valueOf("70D6FF");
            shapes.setColor(Color.valueOf("09090B"));
            shapes.rect(x-13,y-13,27,2);shapes.rect(x-13,y+12,27,2);
            shapes.rect(x-13,y-13,2,27);shapes.rect(x+12,y-13,2,27);
            shapes.setColor(c);
            shapes.rect(x-12,y-12,7,1);shapes.rect(x+5,y-12,7,1);
            shapes.rect(x-12,y+11,7,1);shapes.rect(x+5,y+11,7,1);
            shapes.rect(x-12,y-12,1,7);shapes.rect(x-12,y+5,1,7);
            shapes.rect(x+11,y-12,1,7);shapes.rect(x+11,y+5,1,7);
        }

        panel(7,7,37,28,"26303A");
        panel(W-65,7,58,28,"26303A");
        panel(W-61,H-35,54,28,"26303A");

        if(mode==Mode.SPAWN){
            panel(0,0,W,116,"171419");
            int cols=4;
            float cell=W/(float)cols;
            for(int i=0;i<PropType.values().length;i++){
                int row=i/cols,col=i%cols;
                panel(col*cell+3,5+(2-row)*35,cell-6,31,(i&1)==0?"2C201C":"252027");
            }
        }

        if(mode==Mode.SETTINGS)panel(W-154,H-148,147,141,"161A20");

        if(mode==Mode.CONTEXT){
            float cx=MathUtils.clamp(contextX,68,W-68),cy=MathUtils.clamp(contextY,48,H-48);
            panel(cx-66,cy-41,132,82,"17151A");
            shapes.setColor(Color.valueOf("5A321F"));shapes.rect(cx-1,cy-40,2,80);shapes.rect(cx-65,cy-1,130,2);
        }
        shapes.end();

        batch.setProjectionMatrix(projection);
        batch.begin();
        pixelText("+",20,28,Color.WHITE);
        pixelText("UNDO",W-59,26,Color.WHITE);
        pixelText("MENU",W-55,H-18,Color.WHITE);

        if(grabbed!=null)pixelText(grabbed.type.label+" / "+(directionIndex(grabbed.body.getAngle())*90)+" DEG",8,H-34,Color.WHITE);

        if(mode==Mode.SPAWN)drawSpawnDrawer();
        if(mode==Mode.SETTINGS)drawSettings();
        if(mode==Mode.CONTEXT)drawContext();
        batch.end();
    }

    private void panel(float x,float y,float w,float h,String fill){
        shapes.setColor(Color.valueOf("09090C"));shapes.rect(x-2,y-2,w+4,h+4);
        shapes.setColor(Color.valueOf(fill));shapes.rect(x,y,w,h);
        shapes.setColor(Color.valueOf("6A5A50"));shapes.rect(x,y+h-1,w,1);shapes.rect(x,y,1,h);
        shapes.setColor(Color.valueOf("101218"));shapes.rect(x,y,w,1);shapes.rect(x+w-1,y,1,h);
    }

    private void drawSpawnDrawer(){
        pixelText("SPAWN / EACH ICON HAS 4 AUTHORED ANGLES",8,109,Color.valueOf("F0C06A"));
        int cols=4;
        float cell=W/(float)cols;
        int preview=(int)(TimeUtils.millis()/600L)%4;
        PropType[] vals=PropType.values();
        for(int i=0;i<vals.length;i++){
            int row=i/cols,col=i%cols;
            float cx=col*cell+19;
            float cy=7+(2-row)*35+15;
            Texture tex=spriteSets.get(spriteKey(vals[i],vals[i].defaultMaterial))[preview];
            batch.draw(tex,cx-16,cy-16,32,32);
            pixelText(vals[i].label,col*cell+39,cy+4,Color.WHITE);
            pixelText((preview*90)+"",col*cell+39,cy-7,Color.valueOf("8C9AA7"));
        }
    }

    private void drawSettings(){
        float x=W-145,y=H-18;
        pixelText("PURE 2D PIXEL",x,y,Color.valueOf("F0C06A"));
        pixelText("RESET WORLD",x,y-32,Color.WHITE);
        pixelText("HAPTICS: "+(haptics?"ON":"OFF"),x,y-62,Color.WHITE);
        pixelText("CANVAS 480x270",x,y-92,Color.valueOf("70D6FF"));
        pixelText("NO 3D / NO MESH",x,y-108,Color.valueOf("70D6FF"));
        pixelText("MENU = CLOSE",x,y-127,Color.valueOf("776F72"));
    }

    private void drawContext(){
        float cx=MathUtils.clamp(contextX,68,W-68),cy=MathUtils.clamp(contextY,48,H-48);
        pixelText(contextProp!=null&&contextProp.frozen?"UNFREEZE":"FREEZE",cx-58,cy+23,Color.valueOf("70D6FF"));
        pixelText("DELETE",cx+10,cy+23,Color.valueOf("FF806C"));
        pixelText("DUPLICATE",cx-58,cy-19,Color.WHITE);
        pixelText("MATERIAL",cx+10,cy-19,Color.valueOf("F0C06A"));
    }

    private void pixelText(String s,float x,float y,Color c){
        font.setColor(Color.valueOf("08080A"));font.draw(batch,s,x+1,y-1);
        font.setColor(c);font.draw(batch,s,x,y);
    }

    private void feedback(){if(haptics)Gdx.input.vibrate(12);}

    private void saveWorld(){
        Array<SaveState> states=new Array<>();
        for(Prop p:props)states.add(p.save());
        prefs.putString("world",json.toJson(states,Array.class,SaveState.class));
        prefs.putBoolean("haptics",haptics);
        prefs.flush();
    }

    @SuppressWarnings("unchecked")
    private boolean restoreWorld(){
        haptics=prefs.getBoolean("haptics",true);
        String data=prefs.getString("world","");
        if(data.isEmpty())return false;
        try{
            Array<SaveState> states=json.fromJson(Array.class,SaveState.class,data);
            if(states==null||states.size==0)return false;
            for(SaveState s:states)spawnFromState(s);
            return true;
        }catch(Exception e){return false;}
    }

    private void spawnFromState(SaveState s){
        spawnInternal(PropType.valueOf(s.type),MaterialKind.valueOf(s.material),
                s.x*PPM,s.y*PPM,s.angle,s.frozen,s.id);
    }

    private void resetWorld(){
        endGrab();
        Array<Prop> copy=new Array<>(props);
        for(Prop p:copy)removeProp(p);
        undo.clear();nextId=1;
        mode=Mode.IDLE;contextProp=null;pressed=null;
        createStarterSet();saveWorld();feedback();
    }

    @Override
    public boolean touchDown(int screenX,int screenY,int pointer,int button){
        if(pointer>=px.length||!insidePresentation(screenX,screenY))return false;
        screenX=vx(screenX);screenY=vy(screenY);
        px[pointer]=screenX;py[pointer]=screenY;

        if(mode==Mode.SPAWN){
            if(screenY>H-116+18){
                int cols=4;
                float cell=W/(float)cols;
                int row=MathUtils.clamp((screenY-(H-116+18))/35,0,2);
                int col=MathUtils.clamp((int)(screenX/cell),0,3);
                int idx=row*cols+col;
                if(idx>=0&&idx<PropType.values().length)spawnWithUndo(PropType.values()[idx]);
            }
            mode=Mode.IDLE;
            return true;
        }

        if(mode==Mode.SETTINGS){
            float ux=screenX,uy=H-screenY;
            if(ux<W-154||uy<H-148){mode=Mode.IDLE;return true;}
            if(uy>H-62&&uy<H-25){resetWorld();mode=Mode.IDLE;return true;}
            if(uy>H-94&&uy<=H-62){haptics=!haptics;feedback();saveWorld();return true;}
            return true;
        }

        if(mode==Mode.CONTEXT){
            float ux=screenX,uy=H-screenY;
            float cx=MathUtils.clamp(contextX,68,W-68),cy=MathUtils.clamp(contextY,48,H-48);
            boolean left=ux<cx,top=uy>cy;
            Prop target=contextProp;
            mode=Mode.IDLE;contextProp=null;
            if(target!=null){
                if(top&&left)setFrozen(target,!target.frozen,true);
                else if(top)deleteWithUndo(target);
                else if(left)duplicateWithUndo(target);
                else cycleMaterial(target);
            }
            return true;
        }

        if(screenX<52&&screenY>H-42){mode=Mode.SPAWN;feedback();return true;}
        if(screenX>W-72&&screenY>H-42){doUndo();return true;}
        if(screenX>W-70&&screenY<42){mode=Mode.SETTINGS;feedback();return true;}

        if(mode==Mode.GRAB&&pointer!=primaryPointer&&secondPointer<0){
            secondPointer=pointer;
            lastTwoAngle=pointerAngle(primaryPointer,secondPointer);
            return true;
        }

        if(pointer==0&&mode==Mode.IDLE){
            primaryPointer=pointer;secondPointer=-1;
            downX=screenX;downY=screenY;
            downNanos=TimeUtils.nanoTime();
            contextTriggered=false;
            float wx=worldX(screenX),wy=worldY(screenY);
            pressed=pick(wx,wy);
            if(pressed!=null&&!pressed.frozen)startGrab(pressed,wx,wy);
            return true;
        }
        return false;
    }

    @Override
    public boolean touchDragged(int screenX,int screenY,int pointer){
        if(pointer>=px.length)return false;
        screenX=vx(screenX);screenY=vy(screenY);
        px[pointer]=screenX;py[pointer]=screenY;

        if(mode==Mode.GRAB&&grabbed!=null){
            if(pointer==primaryPointer)grabTarget.set(worldX(screenX),worldY(screenY));
            if(secondPointer>=0){
                float a=pointerAngle(primaryPointer,secondPointer);
                float da=wrapAngle(a-lastTwoAngle);
                accumulatedTorque+=da;
                lastTwoAngle=a;
            }
            return true;
        }
        return false;
    }

    @Override
    public boolean touchUp(int screenX,int screenY,int pointer,int button){
        screenX=vx(screenX);screenY=vy(screenY);
        if(pointer<px.length){px[pointer]=screenX;py[pointer]=screenY;}
        if(mode==Mode.GRAB){
            if(pointer==secondPointer){secondPointer=-1;return true;}
            if(pointer==primaryPointer){
                endGrab();pressed=null;primaryPointer=-1;return true;
            }
        }
        if(pointer==primaryPointer){
            primaryPointer=-1;pressed=null;return true;
        }
        return false;
    }

    private float pointerAngle(int a,int b){return MathUtils.atan2(py[b]-py[a],px[b]-px[a]);}
    private float wrapAngle(float a){while(a>MathUtils.PI)a-=MathUtils.PI2;while(a<-MathUtils.PI)a+=MathUtils.PI2;return a;}

    @Override public boolean touchCancelled(int x,int y,int pointer,int button){return touchUp(x,y,pointer,button);}
    @Override public boolean keyDown(int keycode){return false;}
    @Override public boolean keyUp(int keycode){return false;}
    @Override public boolean keyTyped(char character){return false;}
    @Override public boolean mouseMoved(int screenX,int screenY){return false;}
    @Override public boolean scrolled(float amountX,float amountY){return false;}

    @Override public void resize(int width,int height){updatePresentation();}
    @Override public void pause(){saveWorld();}

    @Override
    public void dispose(){
        saveWorld();
        for(Texture[] frames:spriteSets.values())for(Texture t:frames)t.dispose();
        if(buffer!=null)buffer.dispose();
        if(batch!=null)batch.dispose();
        if(shapes!=null)shapes.dispose();
        if(font!=null)font.dispose();
        if(world!=null)world.dispose();
    }
}
