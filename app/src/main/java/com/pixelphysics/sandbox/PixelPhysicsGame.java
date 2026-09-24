package com.pixelphysics.sandbox;

import com.badlogic.gdx.*;
import com.badlogic.gdx.files.FileHandle;
import com.badlogic.gdx.graphics.*;
import com.badlogic.gdx.graphics.g2d.*;
import com.badlogic.gdx.graphics.g3d.*;
import com.badlogic.gdx.graphics.g3d.attributes.*;
import com.badlogic.gdx.graphics.g3d.environment.DirectionalLight;
import com.badlogic.gdx.graphics.g3d.utils.ModelBuilder;
import com.badlogic.gdx.graphics.glutils.*;
import com.badlogic.gdx.math.*;
import com.badlogic.gdx.math.collision.Ray;
import com.badlogic.gdx.physics.bullet.Bullet;
import com.badlogic.gdx.physics.bullet.collision.*;
import com.badlogic.gdx.physics.bullet.dynamics.*;
import com.badlogic.gdx.physics.bullet.linearmath.btDefaultMotionState;
import com.badlogic.gdx.utils.*;

import java.util.Locale;
import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;

public class PixelPhysicsGame extends ApplicationAdapter implements InputProcessor {
    private static final float FIXED_DT = 1f / 60f;
    private static final float MAX_ACCUM = 0.12f;
    private static final int MAX_OBJECTS = 100;
    private static final int VIRTUAL_W = 480;
    private static final int VIRTUAL_H = 270;

    private PerspectiveCamera camera;
    private ModelBatch modelBatch;
    private Environment environment;
    private SpriteBatch spriteBatch;
    private ShapeRenderer shapes;
    private BitmapFont font;
    private FrameBuffer frameBuffer;
    private TextureRegion frameRegion;
    private final Matrix4 uiProjection = new Matrix4();
    private float presentationScale = 1f;
    private float presentationX = 0f;
    private float presentationY = 0f;
    private float presentationW = VIRTUAL_W;
    private float presentationH = VIRTUAL_H;

    private Texture mahoganyTex;
    private Texture metalTex;
    private Texture rubberTex;
    private Texture floorTex;
    private Texture workshopWoodTex;
    private Texture darkMetalTex;

    private btDefaultCollisionConfiguration collisionConfig;
    private btCollisionDispatcher dispatcher;
    private btDbvtBroadphase broadphase;
    private btSequentialImpulseConstraintSolver solver;
    private btDiscreteDynamicsWorld world;
    private ContactListener contactListener;

    private final Array<PhysicsObject> objects = new Array<>();
    private final Array<StaticPiece> staticPieces = new Array<>();
    private final ObjectMap<Integer, PhysicsObject> byId = new ObjectMap<>();
    private final Array<UndoAction> undo = new Array<>();
    private int nextId = 1;

    private Preferences prefs;
    private Json json;
    private float autosaveClock = 0f;

    private com.badlogic.gdx.audio.Sound woodSound, metalSound, rubberSound, clickSound;
    private boolean haptics = true;
    private int pendingImpactA = -1, pendingImpactB = -1;

    private final Vector3 pivot = new Vector3(0, 1.2f, 0);
    private float orbitYaw = 35f;
    private float orbitPitch = 23f;
    private float orbitDistance = 10.5f;

    private enum Mode { IDLE, CAMERA, GRAB, CONTEXT, SPAWN, SETTINGS }
    private Mode mode = Mode.IDLE;

    private int primaryPointer = -1;
    private int secondPointer = -1;
    private final int[] px = new int[10];
    private final int[] py = new int[10];
    private float downX, downY;
    private long downNanos;
    private float lastPrimaryX, lastPrimaryY;
    private float lastTwoDistance, lastTwoAngle;
    private float lastTwoMidX, lastTwoMidY;
    private float tapX, tapY;
    private long lastTapNanos = 0;

    private PhysicsObject grabbed;
    private final Vector3 grabLocal = new Vector3();
    private float grabDepth = 4f;
    private float accumulatedTwist = 0f;
    private boolean contextTriggered = false;
    private float contextX, contextY;
    private PhysicsObject contextObject;
    private PhysicsObject pressObject;

    private final Vector3 tmp1 = new Vector3();
    private final Vector3 tmp2 = new Vector3();
    private final Vector3 tmp3 = new Vector3();
    private final Quaternion tmpQ = new Quaternion();
    private final Matrix4 tmpM = new Matrix4();

    private float accumulator = 0f;

    private enum MaterialKind { MAHOGANY, METAL, RUBBER }
    private enum PixelSurface { FLOOR, WOOD, DARK_METAL }

    private enum ObjectType {
        CUBE("Cube", 1.0f, 1.0f, 1.0f, false, MaterialKind.MAHOGANY, 1.0f),
        BEAM_SHORT("Beam S", 2.4f, 0.45f, 0.45f, false, MaterialKind.MAHOGANY, 1.3f),
        BEAM_LONG("Beam L", 4.0f, 0.45f, 0.45f, false, MaterialKind.MAHOGANY, 2.0f),
        PLANK("Plank", 2.8f, 0.22f, 1.2f, false, MaterialKind.MAHOGANY, 1.2f),
        WOOD_BALL("Wood Ball", 0.85f, 0.85f, 0.85f, true, MaterialKind.MAHOGANY, 0.7f),
        METAL_BALL("Metal Ball", 0.85f, 0.85f, 0.85f, true, MaterialKind.METAL, 4.0f),
        WEIGHT("Weight", 0.9f, 0.9f, 0.9f, false, MaterialKind.METAL, 12.0f),
        WHEEL("Wheel", 1.2f, 0.35f, 1.2f, false, MaterialKind.MAHOGANY, 1.1f),
        RUBBER_BALL("Rubber", 0.95f, 0.95f, 0.95f, true, MaterialKind.RUBBER, 0.85f),
        BARREL("Barrel", 0.95f, 1.5f, 0.95f, false, MaterialKind.METAL, 3.8f),
        CRATE("Crate", 1.35f, 1.35f, 1.35f, false, MaterialKind.MAHOGANY, 1.7f),
        RAMP("Ramp", 2.8f, 0.32f, 1.4f, false, MaterialKind.MAHOGANY, 1.8f);

        final String label;
        final float w,h,d;
        final boolean sphere;
        final MaterialKind material;
        final float mass;
        ObjectType(String label, float w, float h, float d, boolean sphere, MaterialKind material, float mass) {
            this.label=label; this.w=w; this.h=h; this.d=d; this.sphere=sphere; this.material=material; this.mass=mass;
        }
    }

    private static class MaterialProfile {
        final float friction, restitution;
        final Color color;
        MaterialProfile(float friction, float restitution, Color color) {
            this.friction=friction; this.restitution=restitution; this.color=color;
        }
    }

    private static class SaveState {
        public int id;
        public String type;
        public String material;
        public float x,y,z,qx,qy,qz,qw;
        public boolean frozen;
        public SaveState() {}
    }

    private class PhysicsObject {
        int id;
        ObjectType type;
        MaterialKind material;
        float mass;
        btCollisionShape shape;
        btDefaultMotionState motionState;
        btRigidBody body;
        Model model;
        ModelInstance instance;
        boolean frozen;
        float lastImpactAt = -100f;

        SaveState snapshot() {
            SaveState s = new SaveState();
            s.id=id; s.type=type.name(); s.material=material.name(); s.frozen=frozen;
            body.getWorldTransform(tmpM);
            tmpM.getTranslation(tmp1);
            tmpM.getRotation(tmpQ, true);
            s.x=tmp1.x; s.y=tmp1.y; s.z=tmp1.z;
            s.qx=tmpQ.x; s.qy=tmpQ.y; s.qz=tmpQ.z; s.qw=tmpQ.w;
            return s;
        }
    }

    private static class StaticPiece {
        btCollisionShape shape;
        btRigidBody body;
        btDefaultMotionState motionState;
        Model model;
        ModelInstance instance;
    }

    private interface UndoAction { void undo(); }

    private static class RayHit {
        PhysicsObject object;
        final Vector3 point = new Vector3();
        boolean hit;
    }

    @Override
    public void create() {
        Locale.setDefault(Locale.US);
        Bullet.init();
        prefs = Gdx.app.getPreferences("pixel-physics-world-v1");
        json = new Json();

        modelBatch = new ModelBatch();
        spriteBatch = new SpriteBatch();
        shapes = new ShapeRenderer();
        font = new BitmapFont();
        font.getData().setScale(0.72f);
        font.getRegion().getTexture().setFilter(Texture.TextureFilter.Nearest, Texture.TextureFilter.Nearest);

        createPixelTextures();

        camera = new PerspectiveCamera(58f, VIRTUAL_W, VIRTUAL_H);
        camera.near = 0.1f;
        camera.far = 100f;
        updateCamera();

        environment = new Environment();
        environment.set(new ColorAttribute(ColorAttribute.AmbientLight, 0.52f, 0.46f, 0.40f, 1f));
        environment.add(new DirectionalLight().set(0.92f, 0.82f, 0.68f, -0.55f, -1f, -0.3f));

        collisionConfig = new btDefaultCollisionConfiguration();
        dispatcher = new btCollisionDispatcher(collisionConfig);
        broadphase = new btDbvtBroadphase();
        solver = new btSequentialImpulseConstraintSolver();
        world = new btDiscreteDynamicsWorld(dispatcher, broadphase, solver, collisionConfig);
        world.setGravity(new Vector3(0, -9.81f, 0));

        contactListener = new ContactListener() {
            @Override
            public void onContactStarted(int userValue0, int userValue1) {
                pendingImpactA = userValue0;
                pendingImpactB = userValue1;
            }
        };
        contactListener.enableOnStarted();

        loadSounds();
        buildWorkshop();
        if (!restoreWorld()) createStarterSet();

        recreateFrameBuffer();
        updatePresentation();
        Gdx.input.setInputProcessor(this);
    }

    private void createPixelTextures() {
        mahoganyTex = makeMahoganyTexture();
        metalTex = makeMetalTexture();
        rubberTex = makeRubberTexture();
        floorTex = makeFloorTexture();
        workshopWoodTex = makeWorkshopWoodTexture();
        darkMetalTex = makeDarkMetalTexture();
    }

    private Texture finishPixelTexture(Pixmap p) {
        Texture t = new Texture(p);
        p.dispose();
        t.setFilter(Texture.TextureFilter.Nearest, Texture.TextureFilter.Nearest);
        t.setWrap(Texture.TextureWrap.Repeat, Texture.TextureWrap.Repeat);
        return t;
    }

    private Texture makeMahoganyTexture() {
        Pixmap p = new Pixmap(16,16, Pixmap.Format.RGBA8888);
        Color base=Color.valueOf("6E2418"), dark=Color.valueOf("3C1512"), mid=Color.valueOf("963A20"),
                light=Color.valueOf("C65C2D"), hi=Color.valueOf("E6813E");
        p.setColor(base); p.fill();
        for(int y=1;y<16;y+=4){ p.setColor(dark); p.drawLine(0,y,15,y); }
        for(int y=2;y<16;y+=4){ p.setColor(mid); p.drawLine((y*3)%7,y,Math.min(15,(y*3)%7+7),y); }
        for(int i=0;i<18;i++){
            int x=(i*7+3)%16, y=(i*11+5)%16;
            p.setColor((i%3)==0?hi:light); p.drawPixel(x,y);
            if((i&1)==0 && x<15) p.drawPixel(x+1,y);
        }
        p.setColor(dark); p.drawRectangle(0,0,16,16);
        return finishPixelTexture(p);
    }

    private Texture makeMetalTexture() {
        Pixmap p = new Pixmap(16,16, Pixmap.Format.RGBA8888);
        Color dark=Color.valueOf("252B31"), base=Color.valueOf("4A5660"), mid=Color.valueOf("71808A"),
                light=Color.valueOf("A7B1B7"), hi=Color.valueOf("D6DADD");
        p.setColor(base); p.fill();
        for(int y=0;y<16;y++){
            p.setColor((y%5)==0?dark:((y%5)==1?mid:base));
            p.drawLine(0,y,15,y);
        }
        p.setColor(light); p.drawLine(2,3,12,3); p.drawLine(4,9,14,9);
        p.setColor(hi); p.drawPixel(3,3); p.drawPixel(10,9); p.drawPixel(12,9);
        p.setColor(dark);
        p.drawPixel(1,1); p.drawPixel(14,1); p.drawPixel(1,14); p.drawPixel(14,14);
        p.drawRectangle(0,0,16,16);
        return finishPixelTexture(p);
    }

    private Texture makeRubberTexture() {
        Pixmap p = new Pixmap(16,16, Pixmap.Format.RGBA8888);
        Color dark=Color.valueOf("17201B"), base=Color.valueOf("28352B"), mid=Color.valueOf("3C4D3E"),
                light=Color.valueOf("566D57");
        p.setColor(base); p.fill();
        for(int y=0;y<16;y++) for(int x=0;x<16;x++) {
            int k=(x*5+y*7)%17;
            if(k==0){ p.setColor(light); p.drawPixel(x,y); }
            else if(k==1){ p.setColor(dark); p.drawPixel(x,y); }
        }
        p.setColor(mid); p.drawLine(1,4,14,4); p.drawLine(1,11,14,11);
        p.setColor(dark); p.drawRectangle(0,0,16,16);
        return finishPixelTexture(p);
    }

    private Texture makeFloorTexture() {
        Pixmap p = new Pixmap(32,32, Pixmap.Format.RGBA8888);
        Color dark=Color.valueOf("241512"), base=Color.valueOf("503022"), mid=Color.valueOf("6B4027"),
                light=Color.valueOf("8A5530"), hi=Color.valueOf("A96C39");
        p.setColor(base); p.fill();
        for(int y=0;y<32;y+=8){ p.setColor(dark); p.drawLine(0,y,31,y); }
        for(int band=0;band<4;band++){
            int y=band*8;
            int seam=(band%2==0)?12:22;
            p.setColor(dark); p.drawLine(seam,y,seam,y+7);
            p.setColor(light); p.drawLine(1,y+2,9,y+2); p.drawLine(17,y+5,29,y+5);
            p.setColor(hi); p.drawPixel((band*9+5)%31,y+3);
            p.setColor(mid); p.drawPixel((band*7+19)%31,y+6);
        }
        return finishPixelTexture(p);
    }

    private Texture makeWorkshopWoodTexture() {
        Pixmap p = new Pixmap(16,16, Pixmap.Format.RGBA8888);
        Color dark=Color.valueOf("2C1912"), base=Color.valueOf("5A321F"), mid=Color.valueOf("784526"),
                light=Color.valueOf("9B6033");
        p.setColor(base); p.fill();
        for(int y=2;y<16;y+=5){ p.setColor(dark); p.drawLine(0,y,15,y); }
        p.setColor(mid); p.drawLine(2,5,10,5); p.drawLine(6,12,15,12);
        p.setColor(light); p.drawPixel(3,5); p.drawPixel(11,12); p.drawPixel(13,7);
        p.setColor(dark); p.drawRectangle(0,0,16,16);
        return finishPixelTexture(p);
    }

    private Texture makeDarkMetalTexture() {
        Pixmap p = new Pixmap(16,16, Pixmap.Format.RGBA8888);
        Color dark=Color.valueOf("141820"), base=Color.valueOf("252C36"), mid=Color.valueOf("394554"),
                light=Color.valueOf("596879");
        p.setColor(base); p.fill();
        for(int x=0;x<16;x+=4){ p.setColor(mid); p.drawLine(x,0,x,15); }
        p.setColor(light); p.drawPixel(2,2); p.drawPixel(13,2); p.drawPixel(2,13); p.drawPixel(13,13);
        p.setColor(dark); p.drawRectangle(0,0,16,16);
        return finishPixelTexture(p);
    }

    private Texture textureFor(MaterialKind kind) {
        switch(kind){
            case METAL: return metalTex;
            case RUBBER: return rubberTex;
            default: return mahoganyTex;
        }
    }

    private Texture textureFor(PixelSurface surface) {
        switch(surface){
            case FLOOR: return floorTex;
            case DARK_METAL: return darkMetalTex;
            default: return workshopWoodTex;
        }
    }

    private MaterialProfile profile(MaterialKind kind) {
        switch (kind) {
            case METAL: return new MaterialProfile(0.55f, 0.08f, new Color(0.46f,0.49f,0.52f,1f));
            case RUBBER: return new MaterialProfile(0.9f, 0.76f, new Color(0.22f,0.25f,0.20f,1f));
            default: return new MaterialProfile(0.78f, 0.10f, new Color(0.43f,0.15f,0.075f,1f));
        }
    }

    private void loadSounds() {
        try { woodSound = makeTone("wood",180f,0.12f,0.45f,18f); } catch (Exception ignored) {}
        try { metalSound = makeTone("metal",540f,0.18f,0.32f,12f); } catch (Exception ignored) {}
        try { rubberSound = makeTone("rubber",110f,0.10f,0.42f,22f); } catch (Exception ignored) {}
        try { clickSound = makeTone("click",760f,0.045f,0.24f,35f); } catch (Exception ignored) {}
    }

    private com.badlogic.gdx.audio.Sound makeTone(String name,float freq,float duration,float amp,float decay) throws IOException {
        final int sampleRate=22050;
        int count=Math.max(64,(int)(sampleRate*duration));
        ByteArrayOutputStream pcm=new ByteArrayOutputStream(count*2);
        for(int i=0;i<count;i++){
            float t=i/(float)sampleRate;
            float env=(float)Math.exp(-decay*t);
            short v=(short)(MathUtils.clamp((float)Math.sin(MathUtils.PI2*freq*t)*env*amp,-1f,1f)*32767);
            pcm.write(v & 255); pcm.write((v>>>8)&255);
        }
        byte[] audio=pcm.toByteArray();
        ByteArrayOutputStream wav=new ByteArrayOutputStream(audio.length+44);
        DataOutputStream d=new DataOutputStream(wav);
        writeAscii(d,"RIFF"); writeLE32(d,36+audio.length); writeAscii(d,"WAVEfmt ");
        writeLE32(d,16); writeLE16(d,1); writeLE16(d,1); writeLE32(d,sampleRate);
        writeLE32(d,sampleRate*2); writeLE16(d,2); writeLE16(d,16);
        writeAscii(d,"data"); writeLE32(d,audio.length); d.write(audio); d.flush();
        FileHandle f=Gdx.files.local("pps-audio/"+name+".wav");
        f.writeBytes(wav.toByteArray(),false);
        return Gdx.audio.newSound(f);
    }

    private void writeAscii(DataOutputStream d,String s) throws IOException { for(int i=0;i<s.length();i++) d.writeByte((byte)s.charAt(i)); }
    private void writeLE16(DataOutputStream d,int v) throws IOException { d.writeByte(v&255); d.writeByte((v>>>8)&255); }
    private void writeLE32(DataOutputStream d,int v) throws IOException { d.writeByte(v&255); d.writeByte((v>>>8)&255); d.writeByte((v>>>16)&255); d.writeByte((v>>>24)&255); }

    private void buildWorkshop() {
        createStaticBox(12f, 0.4f, 12f, 0, -0.2f, 0, 0, PixelSurface.FLOOR);
        createStaticBox(4.2f, 0.35f, 2.5f, -2.2f, 1.15f, -2.0f, 0, PixelSurface.WOOD);
        createStaticBox(1.8f, 0.3f, 2.2f, 3.0f, 0.7f, -1.2f, 0, PixelSurface.WOOD);
        createStaticBox(2.5f, 0.28f, 1.7f, 2.3f, 1.45f, 2.7f, 0, PixelSurface.WOOD);
        createStaticBox(2.7f, 0.3f, 1.3f, -2.7f, 0.45f, 2.8f, -18f, PixelSurface.WOOD);
        createStaticBox(0.35f, 2.0f, 4.0f, 5.3f, 1.0f, 0f, 0, PixelSurface.DARK_METAL);
        createStaticBox(3.8f, 0.18f, 0.75f, 0f, 2.65f, -4.4f, 0, PixelSurface.DARK_METAL);
    }

    private void createStaticBox(float w,float h,float d,float x,float y,float z,float rotZ, PixelSurface surface) {
        StaticPiece p = new StaticPiece();
        p.shape = new btBoxShape(new Vector3(w/2f,h/2f,d/2f));
        Matrix4 tr = new Matrix4().idt().translate(x,y,z).rotate(Vector3.Z, rotZ);
        p.motionState = new btDefaultMotionState(tr);
        btRigidBody.btRigidBodyConstructionInfo ci = new btRigidBody.btRigidBodyConstructionInfo(0f,p.motionState,p.shape,new Vector3());
        p.body = new btRigidBody(ci);
        ci.dispose();
        p.body.setUserValue(0);
        p.body.setFriction(0.82f);
        world.addRigidBody(p.body);
        ModelBuilder mb = new ModelBuilder();
        Material mat = new Material(
                ColorAttribute.createDiffuse(Color.WHITE),
                TextureAttribute.createDiffuse(textureFor(surface)));
        long attrs = VertexAttributes.Usage.Position|VertexAttributes.Usage.Normal|VertexAttributes.Usage.TextureCoordinates;
        p.model = mb.createBox(w,h,d,mat,attrs);
        p.instance = new ModelInstance(p.model, tr);
        staticPieces.add(p);
    }

    private void createStarterSet() {
        spawnInternal(ObjectType.CUBE, MaterialKind.MAHOGANY, new Vector3(-0.8f,1.1f,0.5f), new Quaternion(), false, nextId++);
        spawnInternal(ObjectType.BEAM_SHORT, MaterialKind.MAHOGANY, new Vector3(0.9f,1.1f,0.2f), new Quaternion(Vector3.Y,25f), false, nextId++);
        spawnInternal(ObjectType.RUBBER_BALL, MaterialKind.RUBBER, new Vector3(0f,1.4f,-1.1f), new Quaternion(), false, nextId++);
        spawnInternal(ObjectType.CRATE, MaterialKind.MAHOGANY, new Vector3(1.5f,0.9f,1.8f), new Quaternion(), false, nextId++);
        saveWorld();
    }

    private PhysicsObject spawnInternal(ObjectType type, MaterialKind matKind, Vector3 pos, Quaternion rot, boolean frozen, int requestedId) {
        if (objects.size >= MAX_OBJECTS) return null;
        PhysicsObject o = new PhysicsObject();
        o.id = requestedId > 0 ? requestedId : nextId++;
        nextId = Math.max(nextId, o.id + 1);
        o.type=type; o.material=matKind; o.mass=type.mass;
        MaterialProfile mp = profile(matKind);

        o.shape = type.sphere ? new btSphereShape(type.w/2f) : new btBoxShape(new Vector3(type.w/2f,type.h/2f,type.d/2f));
        Vector3 inertia = new Vector3();
        o.shape.calculateLocalInertia(o.mass, inertia);
        Matrix4 tr = new Matrix4().idt().set(pos, rot);
        o.motionState = new btDefaultMotionState(tr);
        btRigidBody.btRigidBodyConstructionInfo ci = new btRigidBody.btRigidBodyConstructionInfo(o.mass,o.motionState,o.shape,inertia);
        o.body = new btRigidBody(ci);
        ci.dispose();
        o.body.setUserValue(o.id);
        o.body.setFriction(mp.friction);
        o.body.setRestitution(mp.restitution);
        o.body.setDamping(0.06f, 0.12f);
        world.addRigidBody(o.body);

        ModelBuilder mb = new ModelBuilder();
        Material material = new Material(
                ColorAttribute.createDiffuse(Color.WHITE),
                TextureAttribute.createDiffuse(textureFor(matKind)),
                FloatAttribute.createShininess(1f));
        long attrs = VertexAttributes.Usage.Position | VertexAttributes.Usage.Normal | VertexAttributes.Usage.TextureCoordinates;
        if (type.sphere) o.model = mb.createSphere(type.w,type.h,type.d,10,7,material,attrs);
        else o.model = mb.createBox(type.w,type.h,type.d,material,attrs);
        o.instance = new ModelInstance(o.model, tr);
        objects.add(o); byId.put(o.id,o);
        if (frozen) setFrozen(o,true,false);
        return o;
    }

    private void removeObject(PhysicsObject o) {
        if (o == null) return;
        if (grabbed == o) endGrab();
        world.removeRigidBody(o.body);
        byId.remove(o.id);
        objects.removeValue(o,true);
        o.body.dispose(); o.motionState.dispose(); o.shape.dispose(); o.model.dispose();
    }

    private void deleteWithUndo(final PhysicsObject o) {
        if (o == null) return;
        final SaveState s = o.snapshot();
        pushUndo(() -> spawnFromState(s));
        removeObject(o);
        feedbackClick(); saveWorld();
    }

    private void spawnWithUndo(ObjectType type) {
        Vector3 pos = chooseSpawnPosition(type);
        final PhysicsObject o = spawnInternal(type,type.material,pos,new Quaternion(),false,nextId++);
        if (o == null) return;
        final int id = o.id;
        pushUndo(() -> { PhysicsObject target=byId.get(id); if (target!=null) removeObject(target); });
        feedbackClick(); saveWorld();
    }

    private void duplicateWithUndo(PhysicsObject source) {
        if (source==null) return;
        SaveState s = source.snapshot();
        Vector3 pos = new Vector3(s.x+0.35f,s.y+0.35f,s.z+0.35f);
        PhysicsObject o = spawnInternal(ObjectType.valueOf(s.type), MaterialKind.valueOf(s.material), pos,
                new Quaternion(s.qx,s.qy,s.qz,s.qw), s.frozen, nextId++);
        if (o==null) return;
        final int id=o.id;
        pushUndo(() -> { PhysicsObject target=byId.get(id); if(target!=null) removeObject(target); });
        feedbackClick(); saveWorld();
    }

    private void spawnFromState(SaveState s) {
        spawnInternal(ObjectType.valueOf(s.type), MaterialKind.valueOf(s.material), new Vector3(s.x,s.y,s.z),
                new Quaternion(s.qx,s.qy,s.qz,s.qw), s.frozen, s.id);
    }

    private Vector3 chooseSpawnPosition(ObjectType type) {
        Ray ray = camera.getPickRay(VIRTUAL_W/2f, VIRTUAL_H/2f);
        RayHit hit = raycast(ray, 40f);
        Vector3 p = hit.hit ? hit.point.cpy() : ray.origin.cpy().mulAdd(ray.direction, 5f);
        p.y += type.h*0.5f + 0.35f;
        for (int i=0;i<objects.size;i++) {
            Vector3 q = objects.get(i).instance.transform.getTranslation(tmp1);
            float min = Math.max(type.w, type.d)*0.55f + Math.max(objects.get(i).type.w, objects.get(i).type.d)*0.55f;
            if (q.dst2(p) < min*min) p.y += type.h + 0.45f;
        }
        return p;
    }

    private void pushUndo(UndoAction a) {
        undo.add(a);
        while (undo.size > 32) undo.removeIndex(0);
    }

    private void doUndo() {
        if (undo.size==0) return;
        UndoAction a=undo.pop();
        a.undo();
        feedbackClick();
        saveWorld();
    }

    private void setFrozen(final PhysicsObject o, boolean frozen, boolean recordUndo) {
        if (o==null || o.frozen==frozen) return;
        final boolean prior=o.frozen;
        if (recordUndo) {
            final int id=o.id;
            pushUndo(() -> {
                PhysicsObject target=byId.get(id);
                if (target!=null) setFrozen(target,prior,false);
            });
        }
        o.frozen=frozen;
        o.body.setLinearVelocity(Vector3.Zero);
        o.body.setAngularVelocity(Vector3.Zero);
        int flags=o.body.getCollisionFlags();
        if (frozen) {
            o.body.setCollisionFlags(flags | btCollisionObject.CollisionFlags.CF_KINEMATIC_OBJECT);
            o.body.setActivationState(Collision.DISABLE_DEACTIVATION);
        } else {
            o.body.setCollisionFlags(flags & ~btCollisionObject.CollisionFlags.CF_KINEMATIC_OBJECT);
            o.body.activate();
        }
        feedbackClick(); saveWorld();
    }

    private void cycleMaterial(final PhysicsObject o) {
        if (o==null) return;
        final MaterialKind prior=o.material;
        MaterialKind next = prior==MaterialKind.MAHOGANY?MaterialKind.METAL:(prior==MaterialKind.METAL?MaterialKind.RUBBER:MaterialKind.MAHOGANY);
        final int id=o.id;
        pushUndo(() -> {
            PhysicsObject target=byId.get(id);
            if (target!=null) applyMaterial(target,prior);
        });
        applyMaterial(o,next);
        feedbackClick();
        saveWorld();
    }

    private void applyMaterial(PhysicsObject o, MaterialKind kind) {
        if (o==null || !byId.containsKey(o.id)) return;
        o.material=kind;
        MaterialProfile mp=profile(kind);
        o.body.setFriction(mp.friction);
        o.body.setRestitution(mp.restitution);
        if (o.instance.materials.size>0) {
            o.instance.materials.get(0).set(ColorAttribute.createDiffuse(Color.WHITE));
            o.instance.materials.get(0).set(TextureAttribute.createDiffuse(textureFor(kind)));
        }
    }

    private void updateGrabPhysics() {
        if (mode!=Mode.GRAB || grabbed==null || grabbed.frozen) return;
        Ray ray = camera.getPickRay(px[primaryPointer], py[primaryPointer]);
        Vector3 planePoint = tmp1.set(camera.position).mulAdd(camera.direction,grabDepth);
        float denom = camera.direction.dot(ray.direction);
        if (Math.abs(denom) < 0.0001f) return;
        float t = camera.direction.dot(tmp2.set(planePoint).sub(ray.origin)) / denom;
        Vector3 target = tmp3.set(ray.origin).mulAdd(ray.direction,t);

        grabbed.body.getWorldTransform(tmpM);
        Vector3 current = new Vector3(grabLocal).mul(tmpM);
        Vector3 center = tmpM.getTranslation(new Vector3());
        Vector3 rel = current.cpy().sub(center);
        Vector3 vel = new Vector3();
        vel.set(grabbed.body.getLinearVelocity());

        float massScale = (float)Math.pow(Math.max(0.25f,grabbed.mass),0.34);
        float kp = 62f*massScale;
        float kd = 10.5f*(float)Math.sqrt(massScale);
        float maxF = 115f*(float)Math.pow(Math.max(1f,grabbed.mass),0.58);
        Vector3 force = target.cpy().sub(current).scl(kp).mulAdd(vel,-kd);
        if (force.len2()>maxF*maxF) force.nor().scl(maxF);
        grabbed.body.applyForce(force,rel);
        grabbed.body.activate();

        if (Math.abs(accumulatedTwist)>0.0001f) {
            float torque = MathUtils.clamp(accumulatedTwist*4.0f,-18f,18f) * (float)Math.pow(grabbed.mass,0.30);
            grabbed.body.applyTorque(camera.direction.cpy().scl(torque));
            accumulatedTwist *= 0.55f;
        }
    }

    private void startGrab(PhysicsObject o, Vector3 hitPoint, int pointer) {
        if (o==null || o.frozen) return;
        mode=Mode.GRAB;
        grabbed=o;
        primaryPointer=pointer;
        secondPointer=-1;
        contextTriggered=false;
        o.body.getWorldTransform(tmpM);
        grabLocal.set(hitPoint).mul(tmpM.cpy().inv());
        grabDepth = MathUtils.clamp(camera.position.dst(hitPoint),1.0f,18f);
        o.body.activate();
    }

    private void endGrab() {
        grabbed=null;
        primaryPointer=-1;
        secondPointer=-1;
        accumulatedTwist=0f;
        if (mode==Mode.GRAB) mode=Mode.IDLE;
    }

    private RayHit raycast(Ray ray,float maxDist) {
        RayHit out=new RayHit();
        Vector3 from=ray.origin.cpy();
        Vector3 to=ray.origin.cpy().mulAdd(ray.direction,maxDist);
        ClosestRayResultCallback cb=new ClosestRayResultCallback(from,to);
        world.rayTest(from,to,cb);
        if (cb.hasHit()) {
            float f=cb.getClosestHitFraction();
            out.point.set(from).lerp(to,f);
            out.hit=true;
            int id=cb.getCollisionObject().getUserValue();
            if (id>0) out.object=byId.get(id);
        }
        cb.dispose();
        return out;
    }

    private void updateCamera() {
        float yaw=orbitYaw*MathUtils.degreesToRadians;
        float pitch=orbitPitch*MathUtils.degreesToRadians;
        float cp=MathUtils.cos(pitch);
        camera.position.set(
                pivot.x + orbitDistance*cp*MathUtils.sin(yaw),
                pivot.y + orbitDistance*MathUtils.sin(pitch),
                pivot.z + orbitDistance*cp*MathUtils.cos(yaw));
        camera.up.set(Vector3.Y);
        camera.lookAt(pivot);
        camera.update();
    }

    private void recreateFrameBuffer() {
        if (frameBuffer!=null) frameBuffer.dispose();
        frameBuffer=new FrameBuffer(Pixmap.Format.RGBA8888,VIRTUAL_W,VIRTUAL_H,true);
        frameBuffer.getColorBufferTexture().setFilter(Texture.TextureFilter.Nearest,Texture.TextureFilter.Nearest);
        frameRegion=new TextureRegion(frameBuffer.getColorBufferTexture());
        frameRegion.flip(false,true);
        uiProjection.setToOrtho2D(0,0,VIRTUAL_W,VIRTUAL_H);
        updatePresentation();
    }

    private void updatePresentation() {
        float raw=Math.min(Gdx.graphics.getWidth()/(float)VIRTUAL_W, Gdx.graphics.getHeight()/(float)VIRTUAL_H);
        float integer=(float)Math.floor(raw);
        presentationScale=integer>=1f?integer:raw;
        presentationW=VIRTUAL_W*presentationScale;
        presentationH=VIRTUAL_H*presentationScale;
        presentationX=(Gdx.graphics.getWidth()-presentationW)*0.5f;
        presentationY=(Gdx.graphics.getHeight()-presentationH)*0.5f;
    }

    private boolean insidePresentation(int sx,int sy) {
        return sx>=presentationX && sx<=presentationX+presentationW &&
               sy>=presentationY && sy<=presentationY+presentationH;
    }

    private int virtualX(int sx) {
        return MathUtils.clamp(Math.round((sx-presentationX)/presentationScale),0,VIRTUAL_W-1);
    }

    private int virtualY(int sy) {
        return MathUtils.clamp(Math.round((sy-presentationY)/presentationScale),0,VIRTUAL_H-1);
    }

    @Override
    public void render() {
        float delta=Math.min(Gdx.graphics.getDeltaTime(),0.05f);
        if ((mode==Mode.GRAB || mode==Mode.CAMERA) && !contextTriggered && pressObject!=null && primaryPointer>=0) {
            float moved=Vector2.dst(downX,downY,px[primaryPointer],py[primaryPointer]);
            if (moved<7f && (TimeUtils.nanoTime()-downNanos)>550_000_000L) {
                contextTriggered=true;
                contextObject=pressObject;
                contextX=downX;
                contextY=downY;
                if (mode==Mode.GRAB) {
                    endGrab();
                } else {
                    primaryPointer=-1;
                    secondPointer=-1;
                }
                mode=Mode.CONTEXT;
            }
        }

        accumulator=Math.min(MAX_ACCUM,accumulator+delta);
        while (accumulator>=FIXED_DT) {
            updateGrabPhysics();
            world.stepSimulation(FIXED_DT,0,FIXED_DT);
            accumulator-=FIXED_DT;
        }
        syncVisuals();
        processImpact();

        frameBuffer.begin();
        Gdx.gl.glViewport(0,0,VIRTUAL_W,VIRTUAL_H);
        Gdx.gl.glEnable(GL20.GL_DEPTH_TEST);
        Gdx.gl.glClearColor(0.055f,0.043f,0.046f,1f);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT|GL20.GL_DEPTH_BUFFER_BIT);
        modelBatch.begin(camera);
        for (StaticPiece p:staticPieces) modelBatch.render(p.instance,environment);
        for (PhysicsObject o:objects) modelBatch.render(o.instance,environment);
        modelBatch.end();

        Gdx.gl.glDisable(GL20.GL_DEPTH_TEST);
        drawHud();
        frameBuffer.end();

        updatePresentation();
        Gdx.gl.glViewport(0,0,Gdx.graphics.getWidth(),Gdx.graphics.getHeight());
        Gdx.gl.glClearColor(0.015f,0.012f,0.015f,1f);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT);
        spriteBatch.setProjectionMatrix(new Matrix4().setToOrtho2D(0,0,Gdx.graphics.getWidth(),Gdx.graphics.getHeight()));
        spriteBatch.begin();
        spriteBatch.draw(frameRegion,presentationX,presentationY,presentationW,presentationH);
        spriteBatch.end();

        autosaveClock+=delta;
        if (autosaveClock>3.0f) {
            autosaveClock=0f;
            saveWorld();
        }
    }

    private void syncVisuals() {
        for (PhysicsObject o:objects) {
            o.body.getWorldTransform(o.instance.transform);
            if (o.instance.transform.getTranslation(tmp1).y < -8f) {
                o.body.setWorldTransform(new Matrix4().idt().setToTranslation(0,4f,0));
                o.body.setLinearVelocity(Vector3.Zero);
                o.body.setAngularVelocity(Vector3.Zero);
            }
        }
    }

    private void processImpact() {
        if (pendingImpactA<0 && pendingImpactB<0) return;
        PhysicsObject a=byId.get(pendingImpactA), b=byId.get(pendingImpactB);
        PhysicsObject src=a!=null?a:b;
        pendingImpactA=pendingImpactB=-1;
        if (src==null) return;
        Vector3 v=new Vector3();
        v.set(src.body.getLinearVelocity());
        float speed=v.len();
        float now=(float)TimeUtils.nanoTime()/1_000_000_000f;
        if (speed<1.25f || now-src.lastImpactAt<0.09f) return;
        src.lastImpactAt=now;
        float vol=MathUtils.clamp((speed-1f)/7f,0.08f,0.75f);
        com.badlogic.gdx.audio.Sound s=src.material==MaterialKind.METAL?metalSound:(src.material==MaterialKind.RUBBER?rubberSound:woodSound);
        if (s!=null) s.play(vol, MathUtils.random(0.94f,1.06f),0f);
        if (haptics && speed>3.3f) Gdx.input.vibrate(MathUtils.clamp((int)(speed*3f),8,35));
    }

    private void feedbackClick() {
        if (clickSound!=null) clickSound.play(0.28f);
        if (haptics) Gdx.input.vibrate(12);
    }

    private void drawHud() {
        final int w=VIRTUAL_W, h=VIRTUAL_H;
        shapes.setProjectionMatrix(uiProjection);
        shapes.begin(ShapeRenderer.ShapeType.Filled);

        drawPixelPanel(7,7,36,28, Color.valueOf("26303A"));
        drawPixelPanel(w-60,7,53,28, Color.valueOf("26303A"));
        drawPixelPanel(w-58,h-35,51,28, Color.valueOf("26303A"));

        if (mode==Mode.SPAWN) {
            drawPixelPanel(0,0,w,122, Color.valueOf("1A161B"));
            ObjectType[] vals=ObjectType.values();
            int cols=4;
            float cellW=w/(float)cols;
            for(int i=0;i<vals.length;i++){
                int row=i/cols,col=i%cols;
                float x=col*cellW+4, y=8+(2-row)*30;
                drawPixelPanel(x,y,cellW-8,25, Color.valueOf(i%2==0?"33241D":"2A2020"));
                shapes.setColor(i%3==0?Color.valueOf("B94E2A"):(i%3==1?Color.valueOf("788692"):Color.valueOf("4E654D")));
                shapes.rect(x+5,y+6,9,9);
                shapes.setColor(Color.valueOf("E6C79C")); shapes.rect(x+6,y+7,3,3);
            }
        }

        if (mode==Mode.SETTINGS) {
            drawPixelPanel(w-154,h-142,147,135, Color.valueOf("171A20"));
        }

        if (mode==Mode.CONTEXT) {
            float cx=MathUtils.clamp(contextX,68,w-68);
            float cy=MathUtils.clamp(h-contextY,50,h-50);
            drawPixelPanel(cx-66,cy-43,132,86, Color.valueOf("17151A"));
            shapes.setColor(Color.valueOf("5A321F")); shapes.rect(cx-1,cy-42,2,84);
            shapes.rect(cx-65,cy-1,130,2);
        }

        drawObjectMarkers();
        shapes.end();

        spriteBatch.setProjectionMatrix(uiProjection);
        spriteBatch.begin();
        pixelText("+",20,27,Color.WHITE);
        pixelText("UNDO",w-53,25,Color.WHITE);
        pixelText("MENU",w-52,h-17,Color.WHITE);
        pixelText("PIXEL PHYSICS",8,h-8,Color.valueOf("F0C06A"));
        pixelText("480x270 // TRUE PIXEL CANVAS",8,h-21,Color.valueOf("8C776B"));

        if (grabbed!=null) {
            pixelText(grabbed.type.label.toUpperCase()+"  "+String.format(Locale.US,"%.1f KG",grabbed.mass),8,h-34,Color.WHITE);
        }

        if (mode==Mode.SPAWN) drawSpawnText(spriteBatch,w,h);
        if (mode==Mode.SETTINGS) drawSettingsText(spriteBatch,w,h);
        if (mode==Mode.CONTEXT) drawContextText(spriteBatch,w,h);
        spriteBatch.end();
    }

    private void drawPixelPanel(float x,float y,float w,float h,Color fill) {
        shapes.setColor(Color.valueOf("09090C")); shapes.rect(x-2,y-2,w+4,h+4);
        shapes.setColor(fill); shapes.rect(x,y,w,h);
        shapes.setColor(Color.valueOf("6A5A50")); shapes.rect(x,y+h-1,w,1); shapes.rect(x,y,1,h);
        shapes.setColor(Color.valueOf("101218")); shapes.rect(x,y,w,1); shapes.rect(x+w-1,y,1,h);
    }

    private void pixelText(String s,float x,float y,Color color) {
        font.setColor(Color.valueOf("08080A"));
        font.draw(spriteBatch,s,x+1,y-1);
        font.setColor(color);
        font.draw(spriteBatch,s,x,y);
    }

    private void drawObjectMarkers() {
        for(PhysicsObject o:objects){
            if(!o.frozen && o!=grabbed) continue;
            Vector3 p=o.instance.transform.getTranslation(new Vector3());
            camera.project(p,0,0,VIRTUAL_W,VIRTUAL_H);
            float x=Math.round(p.x), y=Math.round(p.y);
            Color col=o==grabbed?Color.valueOf("F4D35E"):Color.valueOf("70D6FF");
            shapes.setColor(Color.valueOf("08080A"));
            shapes.rect(x-7,y-7,15,2); shapes.rect(x-7,y+6,15,2);
            shapes.rect(x-7,y-7,2,15); shapes.rect(x+6,y-7,2,15);
            shapes.setColor(col);
            shapes.rect(x-6,y-6,5,1); shapes.rect(x+2,y-6,5,1);
            shapes.rect(x-6,y+5,5,1); shapes.rect(x+2,y+5,5,1);
            shapes.rect(x-6,y-6,1,5); shapes.rect(x-6,y+2,1,5);
            shapes.rect(x+5,y-6,1,5); shapes.rect(x+5,y+2,1,5);
        }
    }

    private void drawSpawnText(SpriteBatch b,int w,int h) {
        pixelText("SPAWN OBJECT",8,115,Color.valueOf("F0C06A"));
        ObjectType[] vals=ObjectType.values();
        int cols=4;
        float cellW=w/(float)cols;
        for(int i=0;i<vals.length;i++) {
            int row=i/cols,col=i%cols;
            float x=col*cellW+20;
            float y=8+(2-row)*30+17;
            pixelText(vals[i].label.toUpperCase(),x,y,Color.WHITE);
        }
    }

    private void drawSettingsText(SpriteBatch b,int w,int h) {
        float x=w-145,y=h-20;
        pixelText("WORKSHOP",x,y,Color.valueOf("F0C06A"));
        pixelText("RESET WORLD",x,y-31,Color.WHITE);
        pixelText("HAPTICS: "+(haptics?"ON":"OFF"),x,y-59,Color.WHITE);
        pixelText("PIXEL GRID: LOCKED",x,y-87,Color.valueOf("70D6FF"));
        pixelText("480 x 270",x,y-103,Color.valueOf("8C9AA7"));
        pixelText("MENU = CLOSE",x,y-120,Color.valueOf("776F72"));
    }

    private void drawContextText(SpriteBatch b,int w,int h) {
        float cx=MathUtils.clamp(contextX,68,w-68);
        float cy=MathUtils.clamp(h-contextY,50,h-50);
        pixelText(contextObject!=null&&contextObject.frozen?"UNFREEZE":"FREEZE",cx-58,cy+24,Color.valueOf("70D6FF"));
        pixelText("DELETE",cx+9,cy+24,Color.valueOf("FF806C"));
        pixelText("DUPLICATE",cx-58,cy-20,Color.WHITE);
        pixelText("MATERIAL",cx+9,cy-20,Color.valueOf("F0C06A"));
    }

    private void saveWorld() {
        Array<SaveState> states=new Array<>();
        for(PhysicsObject o:objects) states.add(o.snapshot());
        prefs.putString("world",json.toJson(states,Array.class,SaveState.class));
        prefs.putBoolean("haptics",haptics);
        prefs.putInteger("pixelCanvasWidth",VIRTUAL_W);
        prefs.flush();
    }

    @SuppressWarnings("unchecked")
    private boolean restoreWorld() {
        haptics=prefs.getBoolean("haptics",true);
        // Pixel canvas is deliberately fixed; old render scale preferences are ignored.
        String data=prefs.getString("world","");
        if (data.isEmpty()) return false;
        try {
            Array<SaveState> states=json.fromJson(Array.class,SaveState.class,data);
            if (states==null || states.size==0) return false;
            for(SaveState s:states) spawnFromState(s);
            return true;
        } catch(Exception e) {
            return false;
        }
    }

    private void resetWorld() {
        endGrab();
        contextObject=null;
        mode=Mode.IDLE;
        Array<PhysicsObject> copy=new Array<>(objects);
        for(PhysicsObject o:copy) removeObject(o);
        undo.clear();
        nextId=1;
        createStarterSet();
        pivot.set(0,1.2f,0);
        orbitYaw=35;
        orbitPitch=23;
        orbitDistance=10.5f;
        updateCamera();
        feedbackClick();
        saveWorld();
    }

    @Override
    public boolean touchDown(int screenX,int screenY,int pointer,int button) {
        if (pointer>=px.length) return false;
        if (!insidePresentation(screenX,screenY)) return false;
        screenX=virtualX(screenX);
        screenY=virtualY(screenY);
        px[pointer]=screenX;
        py[pointer]=screenY;
        int w=VIRTUAL_W,h=VIRTUAL_H;

        if (mode==Mode.SPAWN) {
            int drawerTop=h-122;
            if (screenY>=drawerTop+22) {
                int cols=4;
                float cellW=w/(float)cols;
                int row=MathUtils.clamp((int)((screenY-(drawerTop+22))/30f),0,2);
                int col=MathUtils.clamp((int)(screenX/cellW),0,3);
                int idx=row*cols+col;
                ObjectType[] vals=ObjectType.values();
                if(idx>=0&&idx<vals.length) spawnWithUndo(vals[idx]);
            }
            mode=Mode.IDLE;
            return true;
        }

        if (mode==Mode.SETTINGS) {
            float ux=screenX, uy=h-screenY;
            if (ux<w-154 || uy<h-142) {
                mode=Mode.IDLE;
                return true;
            }
            if (uy>h-60 && uy<h-28) {
                resetWorld();
                mode=Mode.IDLE;
                return true;
            }
            if (uy>h-90 && uy<=h-60) {
                haptics=!haptics;
                feedbackClick();
                saveWorld();
                return true;
            }
            return true;
        }

        if (mode==Mode.CONTEXT) {
            float cy=h-contextY;
            float ux=screenX, uy=h-screenY;
            boolean left=ux<contextX, top=uy>cy;
            PhysicsObject target=contextObject;
            mode=Mode.IDLE;
            contextObject=null;
            if(target!=null) {
                if(top&&left) setFrozen(target,!target.frozen,true);
                else if(top) deleteWithUndo(target);
                else if(left) duplicateWithUndo(target);
                else cycleMaterial(target);
            }
            return true;
        }

        if (screenX<52 && screenY>h-42) {
            mode=Mode.SPAWN;
            feedbackClick();
            return true;
        }
        if (screenX>w-68 && screenY>h-42) {
            doUndo();
            return true;
        }
        if (screenX>w-68 && screenY<42) {
            mode=Mode.SETTINGS;
            feedbackClick();
            return true;
        }

        if (mode==Mode.GRAB && pointer!=primaryPointer && secondPointer<0) {
            secondPointer=pointer;
            lastTwoDistance=distancePointers(primaryPointer,secondPointer);
            lastTwoAngle=anglePointers(primaryPointer,secondPointer);
            return true;
        }

        if (mode==Mode.CAMERA && pointer!=primaryPointer && secondPointer<0) {
            secondPointer=pointer;
            lastTwoDistance=distancePointers(primaryPointer,secondPointer);
            lastTwoMidX=(px[primaryPointer]+px[secondPointer])*0.5f;
            lastTwoMidY=(py[primaryPointer]+py[secondPointer])*0.5f;
            return true;
        }

        if (pointer==0 && (mode==Mode.IDLE || mode==Mode.CAMERA)) {
            primaryPointer=pointer;
            secondPointer=-1;
            downX=screenX;
            downY=screenY;
            lastPrimaryX=screenX;
            lastPrimaryY=screenY;
            downNanos=TimeUtils.nanoTime();
            RayHit hit=raycast(camera.getPickRay(screenX,screenY),50f);
            pressObject=hit.object;
            if (hit.object!=null && !hit.object.frozen) {
                startGrab(hit.object,hit.point,pointer);
            } else {
                mode=Mode.CAMERA;
                long now=TimeUtils.nanoTime();
                if (now-lastTapNanos<330_000_000L && Vector2.dst(tapX,tapY,screenX,screenY)<18f && hit.hit) {
                    pivot.set(hit.point);
                    updateCamera();
                }
                lastTapNanos=now;
                tapX=screenX;
                tapY=screenY;
            }
            return true;
        }
        return false;
    }

    @Override
    public boolean touchDragged(int screenX,int screenY,int pointer) {
        if(pointer>=px.length) return false;
        screenX=virtualX(screenX);
        screenY=virtualY(screenY);
        px[pointer]=screenX;
        py[pointer]=screenY;
        if(mode==Mode.GRAB && grabbed!=null) {
            if(secondPointer>=0) {
                float dist=distancePointers(primaryPointer,secondPointer);
                float dd=dist-lastTwoDistance;
                grabDepth=MathUtils.clamp(grabDepth+dd*0.026f,0.8f,20f);
                float a=anglePointers(primaryPointer,secondPointer);
                float da=wrapAngle(a-lastTwoAngle);
                accumulatedTwist+=da;
                lastTwoDistance=dist;
                lastTwoAngle=a;
            }
            return true;
        }
        if(mode==Mode.CAMERA && primaryPointer>=0) {
            if(secondPointer>=0) {
                float dist=distancePointers(primaryPointer,secondPointer);
                if(lastTwoDistance>2f) orbitDistance=MathUtils.clamp(orbitDistance*(lastTwoDistance/Math.max(2f,dist)),3.0f,22f);
                float mx=(px[primaryPointer]+px[secondPointer])*0.5f;
                float my=(py[primaryPointer]+py[secondPointer])*0.5f;
                float dx=mx-lastTwoMidX, dy=my-lastTwoMidY;
                Vector3 right=tmp1.set(camera.direction).crs(camera.up).nor();
                float panScale=orbitDistance*0.0045f;
                pivot.mulAdd(right,-dx*panScale).mulAdd(camera.up,dy*panScale);
                lastTwoDistance=dist;
                lastTwoMidX=mx;
                lastTwoMidY=my;
            } else if(pointer==primaryPointer) {
                float dx=screenX-lastPrimaryX,dy=screenY-lastPrimaryY;
                orbitYaw-=dx*0.78f;
                orbitPitch=MathUtils.clamp(orbitPitch+dy*0.62f,-12f,72f);
                lastPrimaryX=screenX;
                lastPrimaryY=screenY;
            }
            updateCamera();
            return true;
        }
        return false;
    }

    @Override
    public boolean touchUp(int screenX,int screenY,int pointer,int button) {
        screenX=virtualX(screenX);
        screenY=virtualY(screenY);
        if(pointer<px.length){
            px[pointer]=screenX;
            py[pointer]=screenY;
        }
        if(mode==Mode.GRAB) {
            if(pointer==secondPointer) {
                secondPointer=-1;
                return true;
            }
            if(pointer==primaryPointer) {
                endGrab();
                pressObject=null;
                return true;
            }
        }
        if(mode==Mode.CAMERA) {
            if(pointer==secondPointer) {
                secondPointer=-1;
                return true;
            }
            if(pointer==primaryPointer) {
                primaryPointer=-1;
                secondPointer=-1;
                pressObject=null;
                mode=Mode.IDLE;
                return true;
            }
        }
        return false;
    }

    private float distancePointers(int a,int b){return Vector2.dst(px[a],py[a],px[b],py[b]);}
    private float anglePointers(int a,int b){return MathUtils.atan2(py[b]-py[a],px[b]-px[a]);}
    private float wrapAngle(float a){while(a>MathUtils.PI)a-=MathUtils.PI2;while(a<-MathUtils.PI)a+=MathUtils.PI2;return a;}

    @Override public boolean touchCancelled(int x,int y,int pointer,int button){return touchUp(x,y,pointer,button);}
    @Override public boolean keyDown(int keycode){return false;}
    @Override public boolean keyUp(int keycode){return false;}
    @Override public boolean keyTyped(char character){return false;}
    @Override public boolean mouseMoved(int screenX,int screenY){return false;}
    @Override public boolean scrolled(float amountX,float amountY){
        orbitDistance=MathUtils.clamp(orbitDistance+amountY,3f,22f);
        updateCamera();
        return true;
    }

    @Override
    public void resize(int width,int height) {
        camera.viewportWidth=VIRTUAL_W;
        camera.viewportHeight=VIRTUAL_H;
        updateCamera();
        updatePresentation();
    }

    @Override
    public void pause() {
        saveWorld();
    }

    @Override
    public void dispose() {
        saveWorld();
        if(contactListener!=null) contactListener.dispose();
        for(PhysicsObject o:new Array<>(objects)) removeObject(o);
        for(StaticPiece p:staticPieces){
            world.removeRigidBody(p.body);
            p.body.dispose();
            p.motionState.dispose();
            p.shape.dispose();
            p.model.dispose();
        }
        world.dispose();
        solver.dispose();
        broadphase.dispose();
        dispatcher.dispose();
        collisionConfig.dispose();
        if(frameBuffer!=null) frameBuffer.dispose();
        if(modelBatch!=null) modelBatch.dispose();
        if(spriteBatch!=null) spriteBatch.dispose();
        if(shapes!=null) shapes.dispose();
        if(font!=null) font.dispose();
        if(mahoganyTex!=null) mahoganyTex.dispose();
        if(metalTex!=null) metalTex.dispose();
        if(rubberTex!=null) rubberTex.dispose();
        if(floorTex!=null) floorTex.dispose();
        if(workshopWoodTex!=null) workshopWoodTex.dispose();
        if(darkMetalTex!=null) darkMetalTex.dispose();
        if(woodSound!=null)woodSound.dispose();
        if(metalSound!=null)metalSound.dispose();
        if(rubberSound!=null)rubberSound.dispose();
        if(clickSound!=null)clickSound.dispose();
    }
}
