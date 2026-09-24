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

    private PerspectiveCamera camera;
    private ModelBatch modelBatch;
    private Environment environment;
    private SpriteBatch spriteBatch;
    private ShapeRenderer shapes;
    private BitmapFont font;
    private FrameBuffer frameBuffer;
    private TextureRegion frameRegion;
    private int renderDivisor = 2;

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
        font.getData().setScale(1.05f);

        camera = new PerspectiveCamera(58f, Gdx.graphics.getWidth(), Gdx.graphics.getHeight());
        camera.near = 0.1f;
        camera.far = 100f;
        updateCamera();

        environment = new Environment();
        environment.set(new ColorAttribute(ColorAttribute.AmbientLight, 0.48f, 0.44f, 0.42f, 1f));
        environment.add(new DirectionalLight().set(0.95f, 0.88f, 0.78f, -0.55f, -1f, -0.3f));

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
        Gdx.input.setInputProcessor(this);
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
        createStaticBox(12f, 0.4f, 12f, 0, -0.2f, 0, 0, Color.valueOf("4B4540"));
        createStaticBox(4.2f, 0.35f, 2.5f, -2.2f, 1.15f, -2.0f, 0, Color.valueOf("5B4031"));
        createStaticBox(1.8f, 0.3f, 2.2f, 3.0f, 0.7f, -1.2f, 0, Color.valueOf("544B44"));
        createStaticBox(2.5f, 0.28f, 1.7f, 2.3f, 1.45f, 2.7f, 0, Color.valueOf("544B44"));
        createStaticBox(2.7f, 0.3f, 1.3f, -2.7f, 0.45f, 2.8f, -18f, Color.valueOf("5B4031"));
        createStaticBox(0.35f, 2.0f, 4.0f, 5.3f, 1.0f, 0f, 0, Color.valueOf("39383B"));
    }

    private void createStaticBox(float w,float h,float d,float x,float y,float z,float rotZ, Color color) {
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
        p.model = mb.createBox(w,h,d,new Material(ColorAttribute.createDiffuse(color)), VertexAttributes.Usage.Position|VertexAttributes.Usage.Normal);
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
        Material material = new Material(ColorAttribute.createDiffuse(mp.color), FloatAttribute.createShininess(matKind==MaterialKind.METAL ? 22f : 5f));
        long attrs = VertexAttributes.Usage.Position | VertexAttributes.Usage.Normal;
        if (type.sphere) o.model = mb.createSphere(type.w,type.h,type.d,14,10,material,attrs);
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
        Ray ray = camera.getPickRay(Gdx.graphics.getWidth()/2f, Gdx.graphics.getHeight()/2f);
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
        if (o.instance.materials.size>0) o.instance.materials.get(0).set(ColorAttribute.createDiffuse(mp.color));
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
        int rw=Math.max(320,Gdx.graphics.getWidth()/Math.max(1,renderDivisor));
        int rh=Math.max(180,Gdx.graphics.getHeight()/Math.max(1,renderDivisor));
        frameBuffer=new FrameBuffer(Pixmap.Format.RGBA8888,rw,rh,true);
        frameBuffer.getColorBufferTexture().setFilter(Texture.TextureFilter.Nearest,Texture.TextureFilter.Nearest);
        frameRegion=new TextureRegion(frameBuffer.getColorBufferTexture());
        frameRegion.flip(false,true);
    }

    @Override
    public void render() {
        float delta=Math.min(Gdx.graphics.getDeltaTime(),0.05f);
        if ((mode==Mode.GRAB || mode==Mode.CAMERA) && !contextTriggered && pressObject!=null && primaryPointer>=0) {
            float moved=Vector2.dst(downX,downY,px[primaryPointer],py[primaryPointer]);
            if (moved<18f && (TimeUtils.nanoTime()-downNanos)>550_000_000L) {
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
        Gdx.gl.glViewport(0,0,frameBuffer.getWidth(),frameBuffer.getHeight());
        Gdx.gl.glEnable(GL20.GL_DEPTH_TEST);
        Gdx.gl.glClearColor(0.095f,0.085f,0.09f,1f);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT|GL20.GL_DEPTH_BUFFER_BIT);
        modelBatch.begin(camera);
        for (StaticPiece p:staticPieces) modelBatch.render(p.instance,environment);
        for (PhysicsObject o:objects) modelBatch.render(o.instance,environment);
        modelBatch.end();
        frameBuffer.end();

        Gdx.gl.glViewport(0,0,Gdx.graphics.getWidth(),Gdx.graphics.getHeight());
        Gdx.gl.glDisable(GL20.GL_DEPTH_TEST);
        spriteBatch.begin();
        spriteBatch.draw(frameRegion,0,0,Gdx.graphics.getWidth(),Gdx.graphics.getHeight());
        spriteBatch.end();
        drawHud();

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
        int w=Gdx.graphics.getWidth(), h=Gdx.graphics.getHeight();
        shapes.begin(ShapeRenderer.ShapeType.Filled);
        shapes.setColor(0.05f,0.045f,0.055f,0.80f);
        shapes.rect(18,18,78,56);
        shapes.rect(w-112,18,94,56);
        shapes.rect(w-112,h-74,94,56);
        if (mode==Mode.SPAWN) {
            shapes.setColor(0.035f,0.032f,0.04f,0.95f);
            shapes.rect(0,0,w,h*0.46f);
        }
        if (mode==Mode.SETTINGS) {
            shapes.setColor(0.035f,0.032f,0.04f,0.95f);
            shapes.rect(w-330,h-330,312,240);
        }
        if (mode==Mode.CONTEXT) {
            float cy=h-contextY;
            shapes.setColor(0.04f,0.035f,0.045f,0.94f);
            shapes.rect(contextX-120,cy-95,240,190);
        }
        shapes.end();

        spriteBatch.begin();
        font.setColor(Color.WHITE);
        font.draw(spriteBatch,"+",50,57);
        font.draw(spriteBatch,"UNDO",w-100,53);
        font.draw(spriteBatch,"MENU",w-100,h-38);
        font.setColor(0.83f,0.72f,0.62f,1f);
        font.draw(spriteBatch,"PIXEL PHYSICS",18,h-22);
        if (grabbed!=null) {
            font.setColor(Color.WHITE);
            font.draw(spriteBatch, grabbed.type.label+"  "+String.format(Locale.US,"%.1f kg",grabbed.mass),18,h-46);
        }
        if (mode==Mode.SPAWN) drawSpawnText(spriteBatch,w,h);
        if (mode==Mode.SETTINGS) drawSettingsText(spriteBatch,w,h);
        if (mode==Mode.CONTEXT) drawContextText(spriteBatch,w,h);
        spriteBatch.end();
    }

    private void drawSpawnText(SpriteBatch b,int w,int h) {
        font.setColor(Color.WHITE);
        font.draw(b,"SPAWN",22,h*0.46f-18);
        ObjectType[] vals=ObjectType.values();
        int cols=4;
        float cellW=w/(float)cols;
        float base=h*0.46f-52;
        for(int i=0;i<vals.length;i++) {
            int row=i/cols,col=i%cols;
            font.setColor(i%2==0?Color.LIGHT_GRAY:Color.WHITE);
            font.draw(b,vals[i].label,col*cellW+20,base-row*54);
        }
        font.setColor(Color.GRAY);
        font.draw(b,"tap outside drawer to close",22,20);
    }

    private void drawSettingsText(SpriteBatch b,int w,int h) {
        float x=w-310,y=h-118;
        font.setColor(Color.WHITE);
        font.draw(b,"WORKSHOP",x,y);
        font.draw(b,"RESET WORLD",x,y-54);
        font.draw(b,"HAPTICS: "+(haptics?"ON":"OFF"),x,y-104);
        font.draw(b,"PIXEL SCALE: 1/"+renderDivisor,x,y-154);
        font.setColor(Color.GRAY);
        font.draw(b,"tap MENU to close",x,y-200);
    }

    private void drawContextText(SpriteBatch b,int w,int h) {
        float cy=h-contextY;
        font.setColor(Color.WHITE);
        font.draw(b,contextObject!=null&&contextObject.frozen?"UNFREEZE":"FREEZE",contextX-100,cy+55);
        font.draw(b,"DELETE",contextX+20,cy+55);
        font.draw(b,"DUPLICATE",contextX-100,cy-35);
        font.draw(b,"MATERIAL",contextX+20,cy-35);
    }

    private void saveWorld() {
        Array<SaveState> states=new Array<>();
        for(PhysicsObject o:objects) states.add(o.snapshot());
        prefs.putString("world",json.toJson(states,Array.class,SaveState.class));
        prefs.putBoolean("haptics",haptics);
        prefs.putInteger("renderDivisor",renderDivisor);
        prefs.flush();
    }

    @SuppressWarnings("unchecked")
    private boolean restoreWorld() {
        haptics=prefs.getBoolean("haptics",true);
        renderDivisor=MathUtils.clamp(prefs.getInteger("renderDivisor",2),1,3);
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
        px[pointer]=screenX;
        py[pointer]=screenY;
        int w=Gdx.graphics.getWidth(),h=Gdx.graphics.getHeight();

        if (mode==Mode.SPAWN) {
            if (screenY > h*0.54f) {
                int cols=4;
                float cellW=w/(float)cols;
                float localY=screenY-h*0.54f+52;
                int row=(int)(localY/54f);
                int col=MathUtils.clamp((int)(screenX/cellW),0,3);
                int idx=row*cols+col;
                ObjectType[] vals=ObjectType.values();
                if(idx>=0&&idx<vals.length) spawnWithUndo(vals[idx]);
                mode=Mode.IDLE;
                return true;
            }
            mode=Mode.IDLE;
            return true;
        }

        if (mode==Mode.SETTINGS) {
            float left=w-330;
            if (screenX<left) {
                mode=Mode.IDLE;
                return true;
            }
            float yFromTop=screenY;
            if (yFromTop>130 && yFromTop<195) {
                resetWorld();
                mode=Mode.IDLE;
                return true;
            }
            if (yFromTop>=195 && yFromTop<250) {
                haptics=!haptics;
                feedbackClick();
                saveWorld();
                return true;
            }
            if (yFromTop>=250 && yFromTop<315) {
                renderDivisor=renderDivisor==1?2:(renderDivisor==2?3:1);
                recreateFrameBuffer();
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

        if (screenX<108 && screenY>h-86) {
            mode=Mode.SPAWN;
            feedbackClick();
            return true;
        }
        if (screenX>w-125 && screenY>h-86) {
            doUndo();
            return true;
        }
        if (screenX>w-125 && screenY<86) {
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
                if (now-lastTapNanos<330_000_000L && Vector2.dst(tapX,tapY,screenX,screenY)<42f && hit.hit) {
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
        px[pointer]=screenX;
        py[pointer]=screenY;
        if(mode==Mode.GRAB && grabbed!=null) {
            if(secondPointer>=0) {
                float dist=distancePointers(primaryPointer,secondPointer);
                float dd=dist-lastTwoDistance;
                grabDepth=MathUtils.clamp(grabDepth+dd*0.009f,0.8f,20f);
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
                float panScale=orbitDistance*0.0016f;
                pivot.mulAdd(right,-dx*panScale).mulAdd(camera.up,dy*panScale);
                lastTwoDistance=dist;
                lastTwoMidX=mx;
                lastTwoMidY=my;
            } else if(pointer==primaryPointer) {
                float dx=screenX-lastPrimaryX,dy=screenY-lastPrimaryY;
                orbitYaw-=dx*0.24f;
                orbitPitch=MathUtils.clamp(orbitPitch+dy*0.20f,-12f,72f);
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
        camera.viewportWidth=width;
        camera.viewportHeight=height;
        updateCamera();
        recreateFrameBuffer();
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
        if(woodSound!=null)woodSound.dispose();
        if(metalSound!=null)metalSound.dispose();
        if(rubberSound!=null)rubberSound.dispose();
        if(clickSound!=null)clickSound.dispose();
    }
}
