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
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

public final class PixelPhysicsVoxelView extends View {
    /*
     * Raster contract
     * ---------------
     * 0.01 m in simulation space = 1 voxel unit.
     * The fixed camera is orthographic, yawed 45 degrees and pitched down 45 degrees.
     *
     * In an ordinary orthographic 45-degree projection one world-axis step lands on
     * sqrt(1/2) of a screen unit. The original technique compensates the Y/Z terms
     * by sqrt(2). After that compensation, and after integer pixel snap, the raster
     * basis is exactly:
     *
     *      px = X - Y
     *      py = (X + Y) / 2 - Z
     *
     * where X/Y/Z are voxel integers. Each exposed top/front face is emitted as one
     * viewport pixel. This class performs that compensated projection directly in
     * software instead of using a GPU matrix, which keeps the APK pure Java/Canvas.
     */
    private static final float SQRT2 = 1.41421356237f;
    private static final float CAMERA_45_COMPONENT = 1f / SQRT2;
    private static final float YZ_PROJECTION_COMPENSATION = SQRT2;
    private static final float VOXEL_METERS = 0.01f;

    private static final int W = 480;
    private static final int H = 270;
    private static final int SPRITE_SIZE = 150;

    private static final float ROOM = 0.90f;
    private static final float WALL_H = 0.65f;
    private static final float ORIGIN_X = 240f;
    private static final float ORIGIN_Y = 160f;

    private static final float FIXED_DT = 1f / 60f;
    private static final float MAX_ACCUM = 0.12f;
    private static final float GRAVITY = 9.81f;
    private static final int MAX_PROPS = 100;

    // Physics v7: sequential impulses + adaptive substeps.
    private static final int VELOCITY_ITERATIONS = 10;
    private static final int POSITION_ITERATIONS = 4;
    private static final int MAX_SUBSTEPS = 6;
    private static final float CONTACT_SLOP = 0.0015f;
    private static final float POSITION_BETA = 0.62f;
    private static final float RESTITUTION_VELOCITY_THRESHOLD = 0.85f;
    private static final float SLEEP_LINEAR = 0.035f;
    private static final float SLEEP_VERTICAL = 0.045f;
    private static final float SLEEP_ANGULAR = 0.08f;
    private static final float SLEEP_TIME = 0.70f;

    private final Paint paint = new Paint();
    private final Paint pixelPaint = new Paint();
    private final Bitmap logicalBitmap = Bitmap.createBitmap(W, H, Bitmap.Config.ARGB_8888);
    private final Canvas logical = new Canvas(logicalBitmap);
    private final Rect srcRect = new Rect(0, 0, W, H);
    private final RectF dstRect = new RectF();

    private float presentationScale = 1f;
    private float presentationX;
    private float presentationY;
    private float presentationW = W;
    private float presentationH = H;

    private final ArrayList<Prop> props = new ArrayList<>();
    private final ArrayList<Prop> drawOrder = new ArrayList<>();
    private final HashMap<Integer, Prop> byId = new HashMap<>();
    private final ArrayDeque<UndoAction> undo = new ArrayDeque<>();
    private final ArrayList<Contact> contacts = new ArrayList<>();

    private final HashMap<PropType, VoxelModel> models = new HashMap<>();
    private final HashMap<String, Bitmap[]> voxelSpriteCache = new HashMap<>();

    private int nextId = 1;
    private final SharedPreferences prefs;

    private long lastFrameNanos;
    private long lastSaveMs;
    private float accumulator;
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

    private float grabLocalX;
    private float grabLocalY;
    private float grabLocalZ;
    private float targetX;
    private float targetY;
    private float targetZ;
    private float targetYaw;
    private float lastTwoDistance;
    private float lastTwoAngle;

    private enum MaterialKind { MAHOGANY, METAL, RUBBER }

    private enum RenderKind { VOXEL, BALL }

    private enum PropType {
        CUBE("CUBE", 18,18,18, 1.0f, MaterialKind.MAHOGANY, RenderKind.VOXEL),
        BEAM_SHORT("BEAM S", 40,10,10, 1.3f, MaterialKind.MAHOGANY, RenderKind.VOXEL),
        BEAM_LONG("BEAM L", 62,10,10, 2.0f, MaterialKind.MAHOGANY, RenderKind.VOXEL),
        PLANK("PLANK", 50,18,6, 1.2f, MaterialKind.MAHOGANY, RenderKind.VOXEL),
        WOOD_BALL("WOOD BALL", 18,18,18, 0.7f, MaterialKind.MAHOGANY, RenderKind.BALL),
        METAL_BALL("METAL BALL", 18,18,18, 4.0f, MaterialKind.METAL, RenderKind.BALL),
        WEIGHT("WEIGHT", 20,20,20, 10.0f, MaterialKind.METAL, RenderKind.VOXEL),
        WHEEL("WHEEL", 26,8,26, 1.1f, MaterialKind.MAHOGANY, RenderKind.VOXEL),
        RUBBER_BALL("RUBBER", 20,20,20, 0.85f, MaterialKind.RUBBER, RenderKind.BALL),
        BARREL("BARREL", 22,22,30, 3.6f, MaterialKind.METAL, RenderKind.VOXEL),
        CRATE("CRATE", 24,24,24, 1.8f, MaterialKind.MAHOGANY, RenderKind.VOXEL),
        STAIRS("STAIRS", 48,26,26, 2.3f, MaterialKind.MAHOGANY, RenderKind.VOXEL),
        SPRING("SPRING", 14,14,34, 1.4f, MaterialKind.METAL, RenderKind.VOXEL);

        final String label;
        final int vx, vy, vz;
        final float mass;
        final MaterialKind defaultMaterial;
        final RenderKind renderKind;

        PropType(String label, int vx, int vy, int vz, float mass,
                 MaterialKind defaultMaterial, RenderKind renderKind) {
            this.label = label;
            this.vx = vx;
            this.vy = vy;
            this.vz = vz;
            this.mass = mass;
            this.defaultMaterial = defaultMaterial;
            this.renderKind = renderKind;
        }

        float wMeters() { return vx * VOXEL_METERS; }
        float dMeters() { return vy * VOXEL_METERS; }
        float hMeters() { return vz * VOXEL_METERS; }
    }

    private static final class Palette {
        final int outline, dark, base, mid, light, hi;
        Palette(String outline, String dark, String base, String mid, String light, String hi) {
            this.outline = rgb(outline);
            this.dark = rgb(dark);
            this.base = rgb(base);
            this.mid = rgb(mid);
            this.light = rgb(light);
            this.hi = rgb(hi);
        }
    }

    private static final class Voxel {
        final short x, y, z;
        final byte tag;
        Voxel(int x, int y, int z, int tag) {
            this.x=(short)x; this.y=(short)y; this.z=(short)z; this.tag=(byte)tag;
        }
    }

    private static final class VoxelModel {
        final int sx, sy, sz;
        final ArrayList<Voxel> voxels = new ArrayList<>();
        VoxelModel(int sx, int sy, int sz) {
            this.sx=sx; this.sy=sy; this.sz=sz;
        }
        VoxelModel add(int x,int y,int z,int tag) {
            voxels.add(new Voxel(x,y,z,tag));
            return this;
        }
    }

    private static final class RotVoxel {
        int x,y,z,tag;
    }

    private static final class FacePixel {
        int x,y,color;
        float depth;
        FacePixel(int x,int y,int color,float depth) {
            this.x=x;this.y=y;this.color=color;this.depth=depth;
        }
    }

    private final class Prop {
        int id;
        PropType type;
        MaterialKind material;

        float x,y,z;
        float vx,vy,vz;
        float yaw;
        float spin;
        boolean frozen;
        boolean sleeping;
        float sleepTimer;

        float radius() {
            return Math.max(type.wMeters(), type.dMeters()) * 0.5f;
        }
        float halfW() { return type.wMeters()*0.5f; }
        float halfD() { return type.dMeters()*0.5f; }
        float bottom() { return z-type.hMeters()*0.5f; }
        float top() { return z+type.hMeters()*0.5f; }
        boolean circleFootprint() {
            return type.renderKind==RenderKind.BALL || type==PropType.BARREL;
        }
        float invMass() {
            return frozen ? 0f : 1f/Math.max(0.001f,type.mass);
        }
        float inertia() {
            float m=Math.max(0.001f,type.mass);
            if(circleFootprint()) {
                float r=Math.max(type.wMeters(),type.dMeters())*0.5f;
                return 0.5f*m*r*r;
            }
            float w=type.wMeters(),d=type.dMeters();
            return m*(w*w+d*d)/12f;
        }
        float invInertia() {
            return frozen ? 0f : 1f/Math.max(0.00001f,inertia());
        }
        void wake() {
            if(frozen)return;
            sleeping=false;
            sleepTimer=0f;
        }

        SaveState snapshot() {
            SaveState s=new SaveState();
            s.id=id;s.type=type.name();s.material=material.name();
            s.x=x;s.y=y;s.z=z;s.yaw=yaw;s.frozen=frozen;
            return s;
        }
    }

    private static final class SaveState {
        int id;
        String type,material;
        float x,y,z,yaw;
        boolean frozen;
    }

    private static final class Contact {
        Prop a;
        Prop b;
        float nx,ny,nz;
        float px,py,pz;
        float penetration;
        float friction;
        float restitution;
    }

    private static final class Hit2 {
        boolean hit;
        float nx,ny,penetration;
        float px,py;
    }

    private interface UndoAction { void undo(); }

    public PixelPhysicsVoxelView(Context context) {
        super(context);
        setFocusable(true);
        setFocusableInTouchMode(true);
        setKeepScreenOn(true);

        paint.setAntiAlias(false);
        paint.setDither(false);
        paint.setFilterBitmap(false);
        paint.setStrokeCap(Paint.Cap.SQUARE);
        paint.setStrokeJoin(Paint.Join.MITER);
        paint.setTypeface(Typeface.create(Typeface.MONOSPACE,Typeface.BOLD));

        pixelPaint.setAntiAlias(false);
        pixelPaint.setDither(false);
        pixelPaint.setFilterBitmap(false);

        prefs=context.getSharedPreferences("pixel_physics_voxel_ortho_v6",Context.MODE_PRIVATE);
        haptics=prefs.getBoolean("haptics",true);

        buildModels();
        if(!restoreWorld()) createStarterSet();

        lastFrameNanos=System.nanoTime();
        lastSaveMs=SystemClock.uptimeMillis();
    }

    private Palette palette(MaterialKind m) {
        switch(m) {
            case METAL:
                return new Palette("11151B","27313A","4D5D68","768894","AEBAC1","E4E8EA");
            case RUBBER:
                return new Palette("101713","1A281F","2B3F31","426047","69806A","AFC3A7");
            default:
                return new Palette("25110D","451813","72291B","9B3D22","CA6231","EE914C");
        }
    }

    private void buildModels() {
        models.put(PropType.CUBE, solidBox(18,18,18,false));
        models.put(PropType.BEAM_SHORT, woodBeam(40,10,10));
        models.put(PropType.BEAM_LONG, woodBeam(62,10,10));
        models.put(PropType.PLANK, woodBeam(50,18,6));
        models.put(PropType.WEIGHT, weightModel());
        models.put(PropType.WHEEL, wheelModel());
        models.put(PropType.BARREL, barrelModel());
        models.put(PropType.CRATE, crateModel());
        models.put(PropType.STAIRS, staircaseModel());
        models.put(PropType.SPRING, springModel());
    }

    private VoxelModel solidBox(int sx,int sy,int sz,boolean accentEdges) {
        VoxelModel m=new VoxelModel(sx,sy,sz);
        for(int z=0;z<sz;z++)for(int y=0;y<sy;y++)for(int x=0;x<sx;x++){
            int boundaries=0;
            if(x==0||x==sx-1)boundaries++;
            if(y==0||y==sy-1)boundaries++;
            if(z==0||z==sz-1)boundaries++;
            int tag=accentEdges&&boundaries>=2?1:0;
            // Keep the full volume in the model. Internal voxels never emit faces,
            // but they are required by the occupancy test to suppress internal faces.
            m.add(x,y,z,tag);
        }
        return m;
    }

    private VoxelModel woodBeam(int sx,int sy,int sz) {
        VoxelModel m=solidBox(sx,sy,sz,false);
        for(Voxel v:m.voxels) {
            // Geometry stays one model; grain comes from deterministic coordinate shading.
        }
        return m;
    }

    private VoxelModel weightModel() {
        VoxelModel m=new VoxelModel(20,20,20);
        for(int z=0;z<20;z++)for(int y=0;y<20;y++)for(int x=0;x<20;x++){
            int inset=(z<4||z>15)?3:0;
            if(x<inset||x>=20-inset||y<inset||y>=20-inset)continue;
            boolean surface=x==inset||x==19-inset||y==inset||y==19-inset||z==0||z==19;
            if(surface)m.add(x,y,z,(x<3||y<3||x>16||y>16)?1:0);
        }
        return m;
    }

    private VoxelModel crateModel() {
        int s=24;
        VoxelModel m=new VoxelModel(s,s,s);
        for(int z=0;z<s;z++)for(int y=0;y<s;y++)for(int x=0;x<s;x++){
            boolean shell=x==0||y==0||z==0||x==s-1||y==s-1||z==s-1;
            if(!shell)continue;
            int edgeBands=0;
            if(x<3||x>=s-3)edgeBands++;
            if(y<3||y>=s-3)edgeBands++;
            if(z<3||z>=s-3)edgeBands++;
            boolean edge=edgeBands>=2;
            boolean brace=(x==y||x+y==s-1) && (z==0||z==s-1);
            m.add(x,y,z,edge||brace?1:0);
        }
        return m;
    }

    private VoxelModel staircaseModel() {
        int sx=48,sy=26,sz=26,steps=8;
        VoxelModel m=new VoxelModel(sx,sy,sz);
        int stepW=sx/steps;
        for(int x=0;x<sx;x++){
            int step=x/stepW;
            int maxZ=Math.min(sz-1,(step+1)*3);
            for(int y=0;y<sy;y++)for(int z=0;z<=maxZ;z++){
                boolean surface=z==maxZ||y==0||y==sy-1||x==0||x==sx-1;
                if(surface)m.add(x,y,z,z==maxZ?2:0);
            }
        }
        return m;
    }

    private VoxelModel wheelModel() {
        int sx=26,sy=8,sz=26;
        VoxelModel m=new VoxelModel(sx,sy,sz);
        float cx=(sx-1)/2f,cz=(sz-1)/2f;
        for(int y=0;y<sy;y++)for(int z=0;z<sz;z++)for(int x=0;x<sx;x++){
            float dx=x-cx,dz=z-cz,r=(float)Math.sqrt(dx*dx+dz*dz);
            boolean rim=r>=9.5f&&r<=12.5f;
            boolean hub=r<=3f;
            boolean spoke=(Math.abs(dx)<1.2f||Math.abs(dz)<1.2f||Math.abs(Math.abs(dx)-Math.abs(dz))<1.1f) && r<10f;
            if(!(rim||hub||spoke))continue;
            boolean surface=y==0||y==sy-1||rim;
            if(surface)m.add(x,y,z,rim?1:(hub?2:0));
        }
        return m;
    }

    private VoxelModel barrelModel() {
        int sx=22,sy=22,sz=30;
        VoxelModel m=new VoxelModel(sx,sy,sz);
        float cx=(sx-1)/2f,cy=(sy-1)/2f;
        for(int z=0;z<sz;z++){
            float bulge=1.5f-(Math.abs(z-(sz-1)/2f)/((sz-1)/2f))*1.5f;
            float radius=8.5f+bulge;
            for(int y=0;y<sy;y++)for(int x=0;x<sx;x++){
                float r=(float)Math.sqrt((x-cx)*(x-cx)+(y-cy)*(y-cy));
                if(r<radius-1f||r>radius+0.2f)continue;
                int tag=(z==4||z==5||z==14||z==15||z==24||z==25)?1:0;
                m.add(x,y,z,tag);
            }
        }
        return m;
    }

    private VoxelModel springModel() {
        int sx=14,sy=14,sz=34;
        VoxelModel m=new VoxelModel(sx,sy,sz);
        float cx=6.5f,cy=6.5f;
        for(int z=0;z<sz;z++){
            float a=z*0.78f;
            int x=Math.round(cx+(float)Math.cos(a)*5f);
            int y=Math.round(cy+(float)Math.sin(a)*5f);
            for(int oy=-1;oy<=1;oy++)for(int ox=-1;ox<=1;ox++){
                int xx=x+ox,yy=y+oy;
                if(xx>=0&&xx<sx&&yy>=0&&yy<sy)m.add(xx,yy,z,1);
            }
        }
        for(int y=2;y<12;y++)for(int x=2;x<12;x++){
            m.add(x,y,0,2);m.add(x,y,1,2);m.add(x,y,32,2);m.add(x,y,33,2);
        }
        return m;
    }

    private static long key(int x,int y,int z) {
        return ((long)(x&0x3ff)<<20)|((long)(y&0x3ff)<<10)|(long)(z&0x3ff);
    }

    private Bitmap[] framesFor(PropType type,MaterialKind material) {
        if(type.renderKind!=RenderKind.VOXEL)return null;
        String cacheKey=type.name()+":"+material.name();
        Bitmap[] cached=voxelSpriteCache.get(cacheKey);
        if(cached!=null)return cached;

        Bitmap[] out=new Bitmap[4];
        VoxelModel model=models.get(type);
        for(int dir=0;dir<4;dir++)out[dir]=renderVoxelModel(model,material,dir);
        voxelSpriteCache.put(cacheKey,out);
        return out;
    }

    private Bitmap renderVoxelModel(VoxelModel model,MaterialKind material,int dir) {
        Bitmap bmp=Bitmap.createBitmap(SPRITE_SIZE,SPRITE_SIZE,Bitmap.Config.ARGB_8888);
        int[] pixels=new int[SPRITE_SIZE*SPRITE_SIZE];
        float[] depths=new float[pixels.length];
        Arrays.fill(depths,Float.POSITIVE_INFINITY);

        int rsx=(dir%2==0)?model.sx:model.sy;
        int rsy=(dir%2==0)?model.sy:model.sx;

        ArrayList<RotVoxel> rv=new ArrayList<>(model.voxels.size());
        HashSet<Long> occupancy=new HashSet<>(model.voxels.size()*2);

        for(Voxel v:model.voxels){
            RotVoxel r=new RotVoxel();
            switch(dir){
                case 1:
                    r.x=v.y;r.y=model.sx-1-v.x;break;
                case 2:
                    r.x=model.sx-1-v.x;r.y=model.sy-1-v.y;break;
                case 3:
                    r.x=model.sy-1-v.y;r.y=v.x;break;
                default:
                    r.x=v.x;r.y=v.y;
            }
            r.z=v.z;r.tag=v.tag;
            rv.add(r);occupancy.add(key(r.x,r.y,r.z));
        }

        Palette pal=palette(material);
        float cx=(rsx-1)*0.5f,cy=(rsy-1)*0.5f,cz=(model.sz-1)*0.5f;
        int ox=SPRITE_SIZE/2,oy=SPRITE_SIZE/2;

        for(RotVoxel v:rv){
            float lx=v.x-cx,ly=v.y-cy,lz=v.z-cz;
            int px=ox+Math.round(lx-ly);
            int py=oy+Math.round((lx+ly)*0.5f-lz);

            // Camera is on the -X/-Y/+Z side: top, -X and -Y are visible.
            if(!occupancy.contains(key(v.x,v.y,v.z+1))){
                int color=faceColor(pal,v.tag,0,v.x,v.y,v.z);
                writeFacePixel(pixels,depths,px,py-1,color,cameraDepth(v.x,v.y,v.z,0));
            }
            if(!occupancy.contains(key(v.x-1,v.y,v.z))){
                int color=faceColor(pal,v.tag,1,v.x,v.y,v.z);
                writeFacePixel(pixels,depths,px-1,py,color,cameraDepth(v.x,v.y,v.z,1));
            }
            if(!occupancy.contains(key(v.x,v.y-1,v.z))){
                int color=faceColor(pal,v.tag,2,v.x,v.y,v.z);
                writeFacePixel(pixels,depths,px+1,py,color,cameraDepth(v.x,v.y,v.z,2));
            }
        }

        bmp.setPixels(pixels,0,SPRITE_SIZE,0,0,SPRITE_SIZE,SPRITE_SIZE);
        return bmp;
    }

    private float cameraDepth(int x,int y,int z,int face) {
        // Same 45-degree camera used by the compensated raster basis.
        float horizontal=(x+y)*CAMERA_45_COMPONENT;
        float vertical=z*CAMERA_45_COMPONENT;
        return horizontal-vertical+(face*0.0001f);
    }

    private int faceColor(Palette p,int tag,int face,int x,int y,int z) {
        int c;
        if(face==0)c=(tag==2?p.hi:(tag==1?p.mid:p.light));
        else if(face==1)c=(tag==2?p.light:(tag==1?p.dark:p.base));
        else c=(tag==2?p.mid:(tag==1?p.outline:p.dark));

        int hash=(x*17+y*31+z*13)&15;
        if(hash==0&&tag==0)c=face==0?p.hi:(face==1?p.mid:p.base);
        return c;
    }

    private void writeFacePixel(int[] pixels,float[] depths,int x,int y,int color,float depth) {
        if(x<0||x>=SPRITE_SIZE||y<0||y>=SPRITE_SIZE)return;
        int idx=y*SPRITE_SIZE+x;
        if(depth<depths[idx]){
            depths[idx]=depth;
            pixels[idx]=color;
        }
    }

    private void createStarterSet() {
        spawnInternal(PropType.CRATE,MaterialKind.MAHOGANY,-0.55f,0.30f,0.13f,0f,false,nextId++);
        spawnInternal(PropType.CUBE,MaterialKind.MAHOGANY,0.42f,-0.12f,0.10f,0f,false,nextId++);
        spawnInternal(PropType.CUBE,MaterialKind.RUBBER,0.42f,-0.12f,0.29f,0f,false,nextId++);
        spawnInternal(PropType.CUBE,MaterialKind.MAHOGANY,0.42f,-0.12f,0.48f,0f,false,nextId++);
        spawnInternal(PropType.WHEEL,MaterialKind.MAHOGANY,-0.28f,-0.56f,0.14f,0f,false,nextId++);
        spawnInternal(PropType.WOOD_BALL,MaterialKind.MAHOGANY,-0.48f,-0.08f,0.30f,0f,false,nextId++);
        spawnInternal(PropType.BEAM_SHORT,MaterialKind.MAHOGANY,0.05f,-0.58f,0.10f,0f,false,nextId++);
        spawnInternal(PropType.STAIRS,MaterialKind.MAHOGANY,0.58f,0.44f,0.14f,(float)Math.PI,false,nextId++);
        saveWorld();
    }

    private Prop spawnInternal(PropType type,MaterialKind material,float x,float y,float z,float yaw,boolean frozen,int id) {
        if(props.size()>=MAX_PROPS)return null;
        Prop p=new Prop();
        p.id=id>0?id:nextId++;
        nextId=Math.max(nextId,p.id+1);
        p.type=type;p.material=material;
        p.x=x;p.y=y;p.z=z;p.yaw=yaw;p.frozen=frozen;
        props.add(p);byId.put(p.id,p);
        return p;
    }

    private void removeProp(Prop p) {
        if(p==null)return;
        if(grabbed==p)endGrab();
        props.remove(p);byId.remove(p.id);
    }

    private void spawnWithUndo(PropType type) {
        float n=(props.size()%5)-2f;
        Prop p=spawnInternal(type,type.defaultMaterial,n*0.09f,0f,0.55f,0f,false,nextId++);
        if(p==null)return;
        final int id=p.id;
        pushUndo(()->{
            Prop q=byId.get(id);
            if(q!=null)removeProp(q);
        });
        feedback();saveWorld();
    }

    private void deleteWithUndo(Prop p) {
        if(p==null)return;
        SaveState s=p.snapshot();
        pushUndo(()->spawnFromState(s));
        removeProp(p);feedback();saveWorld();
    }

    private void duplicateWithUndo(Prop source) {
        if(source==null)return;
        Prop p=spawnInternal(source.type,source.material,source.x+0.08f,source.y-0.08f,source.z+0.08f,
                source.yaw,source.frozen,nextId++);
        if(p==null)return;
        final int id=p.id;
        pushUndo(()->{
            Prop q=byId.get(id);
            if(q!=null)removeProp(q);
        });
        feedback();saveWorld();
    }

    private void setFrozen(Prop p,boolean frozen,boolean record) {
        if(p==null||p.frozen==frozen)return;
        final int id=p.id;
        final boolean prior=p.frozen;
        if(record)pushUndo(()->{
            Prop q=byId.get(id);
            if(q!=null)setFrozen(q,prior,false);
        });
        p.frozen=frozen;
        p.vx=p.vy=p.vz=p.spin=0f;
        feedback();saveWorld();
    }

    private void cycleMaterial(Prop p) {
        if(p==null)return;
        final int id=p.id;
        final MaterialKind prior=p.material;
        MaterialKind next=prior==MaterialKind.MAHOGANY?MaterialKind.METAL:
                (prior==MaterialKind.METAL?MaterialKind.RUBBER:MaterialKind.MAHOGANY);
        pushUndo(()->{
            Prop q=byId.get(id);
            if(q!=null)q.material=prior;
        });
        p.material=next;feedback();saveWorld();
    }

    private void pushUndo(UndoAction action) {
        undo.addLast(action);
        while(undo.size()>32)undo.removeFirst();
    }

    private void doUndo() {
        if(undo.isEmpty())return;
        undo.removeLast().undo();feedback();saveWorld();
    }

    private void physicsStep(float dt) {
        for(int i=0;i<props.size();i++){
            Prop p=props.get(i);
            if(p.frozen)continue;

            if(p!=grabbed)p.vz-=GRAVITY*dt;
            p.x+=p.vx*dt;p.y+=p.vy*dt;p.z+=p.vz*dt;
            p.yaw+=p.spin*dt;p.spin*=0.992f;

            float r=p.radius();
            if(p.x-r<-ROOM){p.x=-ROOM+r;p.vx=Math.abs(p.vx)*restitution(p.material);}
            if(p.x+r> ROOM){p.x= ROOM-r;p.vx=-Math.abs(p.vx)*restitution(p.material);}
            if(p.y-r<-ROOM){p.y=-ROOM+r;p.vy=Math.abs(p.vy)*restitution(p.material);}
            if(p.y+r> ROOM){p.y= ROOM-r;p.vy=-Math.abs(p.vy)*restitution(p.material);}

            float support=supportHeight(p);
            if(p.bottom()<support){
                p.z=support+p.type.hMeters()*0.5f;
                if(p.vz<0)p.vz=-p.vz*restitution(p.material);
                if(Math.abs(p.vz)<0.24f)p.vz=0f;
                float fr=friction(p.material);
                p.vx*=Math.max(0f,1f-fr*dt*3f);
                p.vy*=Math.max(0f,1f-fr*dt*3f);
            }

            if(p.z<-1f){
                p.x=0;p.y=0;p.z=0.55f;p.vx=p.vy=p.vz=0;
            }
        }

        solveHorizontalCollisions();
    }

    private float supportHeight(Prop p) {
        float best=0f;
        for(int i=0;i<props.size();i++){
            Prop q=props.get(i);
            if(q==p)continue;
            float dx=p.x-q.x,dy=p.y-q.y;
            float rr=(p.radius()+q.radius())*0.78f;
            if(dx*dx+dy*dy>rr*rr)continue;
            float top=q.top();
            if(top<=p.z+0.05f&&top>best)best=top;
        }
        return best;
    }

    private void solveHorizontalCollisions() {
        for(int i=0;i<props.size();i++){
            Prop a=props.get(i);
            for(int j=i+1;j<props.size();j++){
                Prop b=props.get(j);
                if(a.frozen&&b.frozen)continue;
                if(a.top()<b.bottom()+0.01f||b.top()<a.bottom()+0.01f)continue;

                float dx=b.x-a.x,dy=b.y-a.y;
                float min=(a.radius()+b.radius())*0.78f;
                float d2=dx*dx+dy*dy;
                if(d2>=min*min)continue;

                float d=(float)Math.sqrt(Math.max(d2,0.000001f));
                float nx=dx/d,ny=dy/d,penetration=min-d;
                float wa=a.frozen?0f:1f,wb=b.frozen?0f:1f,sum=wa+wb;
                if(sum<=0)continue;

                if(!a.frozen){a.x-=nx*penetration*(wa/sum);a.y-=ny*penetration*(wa/sum);}
                if(!b.frozen){b.x+=nx*penetration*(wb/sum);b.y+=ny*penetration*(wb/sum);}

                float rel=(b.vx-a.vx)*nx+(b.vy-a.vy)*ny;
                if(rel<0){
                    float e=Math.min(restitution(a.material),restitution(b.material));
                    float ia=a.frozen?0f:1f/a.type.mass,ib=b.frozen?0f:1f/b.type.mass;
                    float impulse=-(1f+e)*rel/Math.max(0.0001f,ia+ib);
                    if(!a.frozen){a.vx-=impulse*nx*ia;a.vy-=impulse*ny*ia;}
                    if(!b.frozen){b.vx+=impulse*nx*ib;b.vy+=impulse*ny*ib;}
                }
            }
        }
    }

    private float restitution(MaterialKind m) {
        return m==MaterialKind.RUBBER?0.70f:(m==MaterialKind.METAL?0.10f:0.15f);
    }

    private float friction(MaterialKind m) {
        return m==MaterialKind.RUBBER?0.88f:(m==MaterialKind.METAL?0.42f:0.68f);
    }

    private void updateGrab(float dt) {
        if(mode!=Mode.GRAB||grabbed==null||grabbed.frozen)return;
        float mass=Math.max(0.25f,grabbed.type.mass);
        float ms=(float)Math.pow(mass,0.34);
        float kp=58f*ms,kd=9.5f*(float)Math.sqrt(ms);
        float fx=(targetX-grabbed.x)*kp-grabbed.vx*kd;
        float fy=(targetY-grabbed.y)*kp-grabbed.vy*kd;
        float fz=(targetZ-grabbed.z)*kp-grabbed.vz*kd;
        float cap=100f*(float)Math.pow(Math.max(1f,mass),0.58);
        float len=(float)Math.sqrt(fx*fx+fy*fy+fz*fz);
        if(len>cap){
            float s=cap/len;fx*=s;fy*=s;fz*=s;
        }
        grabbed.vx+=fx/mass*dt;grabbed.vy+=fy/mass*dt;grabbed.vz+=fz/mass*dt;
    }

    private PointF project(float xMeters,float yMeters,float zMeters) {
        float x=xMeters/VOXEL_METERS;
        float y=yMeters/VOXEL_METERS;
        float z=zMeters/VOXEL_METERS;

        // 45-degree orthographic camera + sqrt(2) Y/Z compensation -> exact raster basis.
        float sx=(x-y);
        float zProjected=z*CAMERA_45_COMPONENT*YZ_PROJECTION_COMPENSATION; // exactly z
        float sy=(x+y)*0.5f-zProjected;

        return new PointF(
                ORIGIN_X+Math.round(sx),
                ORIGIN_Y+Math.round(sy)
        );
    }

    private PointF unprojectAtZ(float sx,float sy,float zMeters) {
        float z=zMeters/VOXEL_METERS;
        float a=sx-ORIGIN_X;
        float compensatedZ=z*CAMERA_45_COMPONENT*YZ_PROJECTION_COMPENSATION;
        float sum=2f*(sy-ORIGIN_Y+compensatedZ);
        float x=(sum+a)*0.5f;
        float y=(sum-a)*0.5f;
        return new PointF(x*VOXEL_METERS,y*VOXEL_METERS);
    }

    private int directionIndex(float yaw) {
        int q=Math.round((float)Math.toDegrees(yaw)/90f)%4;
        if(q<0)q+=4;
        return q;
    }

    private float depthKey(Prop p) {
        PointF q=project(p.x,p.y,p.z);
        return q.y;
    }

    @Override
    protected void onSizeChanged(int w,int h,int oldw,int oldh) {
        if(w<=0||h<=0)return;
        float raw=Math.min(w/(float)W,h/(float)H);
        float integer=(float)Math.floor(raw);
        presentationScale=integer>=1?integer:raw;
        presentationW=W*presentationScale;
        presentationH=H*presentationScale;
        presentationX=(w-presentationW)*0.5f;
        presentationY=(h-presentationH)*0.5f;
        dstRect.set(presentationX,presentationY,presentationX+presentationW,presentationY+presentationH);
    }

    @Override
    protected void onDraw(Canvas screen) {
        super.onDraw(screen);

        long now=System.nanoTime();
        float delta=Math.min(0.05f,(now-lastFrameNanos)/1_000_000_000f);
        lastFrameNanos=now;

        checkLongPress(now);

        accumulator=Math.min(MAX_ACCUM,accumulator+delta);
        while(accumulator>=FIXED_DT){
            updateGrab(FIXED_DT);
            physicsStep(FIXED_DT);
            accumulator-=FIXED_DT;
        }

        renderLogical();

        screen.drawColor(rgb("0B0E14"));
        pixelPaint.setFilterBitmap(false);
        screen.drawBitmap(logicalBitmap,srcRect,dstRect,pixelPaint);

        long ms=SystemClock.uptimeMillis();
        if(ms-lastSaveMs>3000){
            lastSaveMs=ms;saveWorld();
        }
        postInvalidateOnAnimation();
    }

    private void renderLogical() {
        logical.drawColor(rgb("111821"));
        drawRoom(logical);
        drawShadows(logical);
        drawProps(logical);
        drawHud(logical);
    }

    private void drawRoom(Canvas c) {
        PointF front=project(-ROOM,-ROOM,0);
        PointF right=project( ROOM,-ROOM,0);
        PointF back=project( ROOM, ROOM,0);
        PointF left=project(-ROOM, ROOM,0);

        PointF rightTop=project( ROOM,-ROOM,WALL_H);
        PointF backTop=project( ROOM, ROOM,WALL_H);
        PointF leftTop=project(-ROOM, ROOM,WALL_H);

        fillPoly(c,rgb("4F3C39"),right,back,backTop,rightTop);
        fillPoly(c,rgb("58413A"),left,back,backTop,leftTop);
        fillPoly(c,rgb("A8754C"),front,right,back,left);

        // Every tenth voxel is a floor grid line; the line endpoints use the exact voxel projection.
        for(int cm=-90;cm<=90;cm+=10){
            float m=cm*VOXEL_METERS;
            line(c,project(m,-ROOM,0),project(m,ROOM,0),1,rgb("765039"));
            line(c,project(-ROOM,m,0),project(ROOM,m,0),1,rgb("765039"));
        }

        // Wall block seams are also locked to voxel rows.
        for(int zcm=10;zcm<65;zcm+=10){
            float z=zcm*VOXEL_METERS;
            line(c,project(ROOM,-ROOM,z),project(ROOM,ROOM,z),1,rgb("65504B"));
            line(c,project(-ROOM,ROOM,z),project(ROOM,ROOM,z),1,rgb("695049"));
        }

        // Structural rails.
        line(c,project(-ROOM,ROOM,0.58f),project(ROOM,ROOM,0.58f),6,rgb("632B1A"));
        line(c,project(ROOM,-ROOM,0.58f),project(ROOM,ROOM,0.58f),6,rgb("632B1A"));
        line(c,project(-ROOM,ROOM,0.05f),project(ROOM,ROOM,0.05f),5,rgb("472018"));
        line(c,project(ROOM,-ROOM,0.05f),project(ROOM,ROOM,0.05f),5,rgb("472018"));

        // Window.
        fillPoly(c,rgb("3B251D"),
                project(-0.72f,ROOM,0.22f),project(-0.34f,ROOM,0.22f),
                project(-0.34f,ROOM,0.50f),project(-0.72f,ROOM,0.50f));
        fillPoly(c,rgb("EAB75F"),
                project(-0.68f,ROOM,0.25f),project(-0.38f,ROOM,0.25f),
                project(-0.38f,ROOM,0.47f),project(-0.68f,ROOM,0.47f));
        line(c,project(-0.53f,ROOM,0.25f),project(-0.53f,ROOM,0.47f),3,rgb("572819"));
        line(c,project(-0.68f,ROOM,0.36f),project(-0.38f,ROOM,0.36f),3,rgb("572819"));

        // Blackboard.
        fillPoly(c,rgb("202B39"),
                project(-0.25f,ROOM,0.18f),project(0.28f,ROOM,0.18f),
                project(0.28f,ROOM,0.45f),project(-0.25f,ROOM,0.45f));
        line(c,project(-0.18f,ROOM,0.25f),project(0.18f,ROOM,0.39f),2,rgb("8D8882"));
        line(c,project(-0.14f,ROOM,0.33f),project(0.05f,ROOM,0.25f),2,rgb("8D8882"));

        // Shelves.
        line(c,project(0.35f,ROOM-0.01f,0.43f),project(0.75f,ROOM-0.01f,0.43f),7,rgb("B65D2D"));
        line(c,project(ROOM-0.01f,-0.25f,0.38f),project(ROOM-0.01f,0.25f,0.38f),7,rgb("B65D2D"));

        // Front frame.
        line(c,left,front,9,rgb("552619"));
        line(c,front,right,9,rgb("552619"));
        line(c,left,front,2,rgb("B16035"));
        line(c,front,right,2,rgb("B16035"));

        drawCornerCap(c,left);
        drawCornerCap(c,front);
        drawCornerCap(c,right);

        pixelText(c,"1 VOXEL = 0.01m = 1px",8,14,rgb("F0C06A"));
        pixelText(c,"ORTHO 45° / sqrt(2) YZ COMPENSATION",8,26,rgb("8FA3AD"));
    }

    private void drawCornerCap(Canvas c,PointF p) {
        paint.setStyle(Paint.Style.FILL);
        paint.setColor(rgb("333843"));
        c.drawRect(p.x-7,p.y-7,p.x+7,p.y+7,paint);
        paint.setColor(rgb("757C86"));
        c.drawRect(p.x-5,p.y-5,p.x+5,p.y-3,paint);
        paint.setColor(rgb("151820"));
        c.drawRect(p.x-4,p.y-1,p.x-2,p.y+1,paint);
        c.drawRect(p.x+2,p.y-1,p.x+4,p.y+1,paint);
    }

    private void drawShadows(Canvas c) {
        paint.setStyle(Paint.Style.FILL);
        paint.setColor(Color.argb(55,18,14,18));
        for(Prop p:props){
            PointF q=project(p.x,p.y,0.001f);
            float r=Math.max(3,p.radius()/VOXEL_METERS*0.65f);
            c.drawOval(q.x-r,q.y-r*0.28f,q.x+r,q.y+r*0.28f,paint);
        }
    }

    private void drawProps(Canvas c) {
        drawOrder.clear();drawOrder.addAll(props);
        drawOrder.sort(Comparator.comparingDouble(this::depthKey));

        for(Prop p:drawOrder){
            PointF q=project(p.x,p.y,p.z);
            if(p.type.renderKind==RenderKind.BALL) {
                drawBall(c,p,q);
            } else {
                Bitmap frame=framesFor(p.type,p.material)[directionIndex(p.yaw)];
                c.drawBitmap(frame,Math.round(q.x-SPRITE_SIZE*0.5f),Math.round(q.y-SPRITE_SIZE*0.5f),pixelPaint);
            }
        }

        for(Prop p:props){
            if(p!=grabbed&&!p.frozen)continue;
            PointF q=project(p.x,p.y,p.z);
            int color=p==grabbed?rgb("F4D35E"):rgb("70D6FF");
            paint.setColor(color);paint.setStyle(Paint.Style.STROKE);paint.setStrokeWidth(1);
            c.drawRect(q.x-14,q.y-14,q.x+14,q.y+14,paint);
            paint.setStyle(Paint.Style.FILL);
        }
    }

    private void drawBall(Canvas c,Prop p,PointF q) {
        Palette pal=palette(p.material);
        float r=p.type.vx*0.5f;
        paint.setStyle(Paint.Style.FILL);
        paint.setAntiAlias(false);
        paint.setColor(pal.outline);c.drawCircle(q.x,q.y,r+1,paint);
        paint.setColor(pal.dark);c.drawCircle(q.x,q.y,r,paint);
        paint.setColor(pal.base);c.drawCircle(q.x-1,q.y-1,r-2,paint);
        paint.setColor(pal.mid);c.drawCircle(q.x-2,q.y-2,Math.max(2,r-5),paint);
        int d=directionIndex(p.yaw);
        int ox=(d==1||d==2)?4:-4;
        int oy=(d>=2)?4:-4;
        paint.setColor(pal.light);c.drawRect(q.x+ox-2,q.y+oy-2,q.x+ox+2,q.y+oy+1,paint);
        paint.setColor(pal.hi);c.drawRect(q.x+ox-1,q.y+oy-1,q.x+ox+1,q.y+oy,paint);
    }

    private void drawHud(Canvas c) {
        panel(c,7,H-35,38,28,rgb("25303A"));
        panel(c,W-66,H-35,59,28,rgb("25303A"));
        panel(c,W-61,7,54,28,rgb("25303A"));

        if(mode==Mode.SPAWN){
            panel(c,0,H-118,W,118,rgb("15171D"));
            float cw=W/4f;
            for(int i=0;i<PropType.values().length;i++){
                int row=i/4,col=i%4;
                panel(c,col*cw+3,H-94+row*23,cw-6,21,(i&1)==0?rgb("2C211D"):rgb("22252A"));
            }
        }

        if(mode==Mode.SETTINGS)panel(c,W-160,42,153,145,rgb("14181E"));

        if(mode==Mode.CONTEXT){
            float cx=clamp(contextX,68,W-68),cy=clamp(contextY,48,H-48);
            panel(c,cx-66,cy-41,132,82,rgb("17151A"));
            paint.setColor(rgb("5A321F"));
            c.drawRect(cx-1,cy-40,cx+1,cy+40,paint);
            c.drawRect(cx-65,cy-1,cx+65,cy+1,paint);
        }

        pixelText(c,"+",20,H-16,Color.WHITE);
        pixelText(c,"UNDO",W-60,H-17,Color.WHITE);
        pixelText(c,"MENU",W-55,25,Color.WHITE);

        if(grabbed!=null) {
            pixelText(c,grabbed.type.label+" / CARDINAL "+(directionIndex(grabbed.yaw)*90)+"°",8,39,Color.WHITE);
        }
        if(mode==Mode.SPAWN)drawSpawnDrawer(c);
        if(mode==Mode.SETTINGS)drawSettings(c);
        if(mode==Mode.CONTEXT)drawContext(c);
    }

    private void drawSpawnDrawer(Canvas c) {
        pixelText(c,"SPAWN / ONE VOXEL MODEL, ROTATED CARDINALLY",8,H-105,rgb("F0C06A"));
        int preview=(int)((SystemClock.uptimeMillis()/650L)%4L);
        float cw=W/4f;
        PropType[] types=PropType.values();
        for(int i=0;i<types.length;i++){
            int row=i/4,col=i%4;
            float x=col*cw+5,y=H-94+row*23;
            PropType t=types[i];
            if(t.renderKind==RenderKind.VOXEL){
                Bitmap frame=framesFor(t,t.defaultMaterial)[preview];
                Rect from=new Rect(0,0,SPRITE_SIZE,SPRITE_SIZE);
                RectF to=new RectF(x,y-1,x+22,y+21);
                c.drawBitmap(frame,from,to,pixelPaint);
            } else {
                Palette p=palette(t.defaultMaterial);
                paint.setColor(p.base);c.drawCircle(x+11,y+10,8,paint);
                paint.setColor(p.light);c.drawRect(x+7,y+5,x+10,y+8,paint);
            }
            pixelText(c,t.label,x+27,y+9,Color.WHITE);
            pixelText(c,(preview*90)+"°",x+27,y+18,rgb("8C9AA7"));
        }
    }

    private void drawSettings(Canvas c) {
        float x=W-151,y=56;
        pixelText(c,"VOXEL ORTHO v6",x,y,rgb("F0C06A"));
        pixelText(c,"RESET WORLD",x,y+32,Color.WHITE);
        pixelText(c,"HAPTICS: "+(haptics?"ON":"OFF"),x,y+62,Color.WHITE);
        pixelText(c,"0.01m = 1px",x,y+92,rgb("70D6FF"));
        pixelText(c,"sqrt(2) Y/Z FIX",x,y+108,rgb("70D6FF"));
        pixelText(c,"ZERO NATIVE LIBS",x,y+124,rgb("70D6FF"));
    }

    private void drawContext(Canvas c) {
        float cx=clamp(contextX,68,W-68),cy=clamp(contextY,48,H-48);
        pixelText(c,contextProp!=null&&contextProp.frozen?"UNFREEZE":"FREEZE",cx-58,cy-21,rgb("70D6FF"));
        pixelText(c,"DELETE",cx+10,cy-21,rgb("FF806C"));
        pixelText(c,"DUPLICATE",cx-58,cy+23,Color.WHITE);
        pixelText(c,"MATERIAL",cx+10,cy+23,rgb("F0C06A"));
    }

    private void panel(Canvas c,float x,float y,float w,float h,int fill) {
        paint.setStyle(Paint.Style.FILL);
        paint.setColor(rgb("09090C"));c.drawRect(x-2,y-2,x+w+2,y+h+2,paint);
        paint.setColor(fill);c.drawRect(x,y,x+w,y+h,paint);
        paint.setColor(rgb("6A5A50"));c.drawRect(x,y,x+w,y+1,paint);c.drawRect(x,y,x+1,y+h,paint);
        paint.setColor(rgb("101218"));c.drawRect(x,y+h-1,x+w,y+h,paint);c.drawRect(x+w-1,y,x+w,y+h,paint);
    }

    private void pixelText(Canvas c,String s,float x,float y,int color) {
        paint.setAntiAlias(false);paint.setStyle(Paint.Style.FILL);
        paint.setTypeface(Typeface.create(Typeface.MONOSPACE,Typeface.BOLD));paint.setTextSize(7);
        paint.setColor(rgb("08080A"));c.drawText(s,x+1,y+1,paint);
        paint.setColor(color);c.drawText(s,x,y,paint);
    }

    private void fillPoly(Canvas c,int color,PointF... pts) {
        Path p=new Path();p.moveTo(pts[0].x,pts[0].y);
        for(int i=1;i<pts.length;i++)p.lineTo(pts[i].x,pts[i].y);
        p.close();paint.setStyle(Paint.Style.FILL);paint.setColor(color);c.drawPath(p,paint);
    }

    private void line(Canvas c,PointF a,PointF b,float width,int color) {
        paint.setStyle(Paint.Style.STROKE);paint.setStrokeWidth(width);paint.setStrokeCap(Paint.Cap.SQUARE);
        paint.setAntiAlias(false);paint.setColor(color);
        c.drawLine(Math.round(a.x),Math.round(a.y),Math.round(b.x),Math.round(b.y),paint);
        paint.setStyle(Paint.Style.FILL);paint.setStrokeWidth(1);
    }

    private Prop pick(float sx,float sy) {
        drawOrder.clear();drawOrder.addAll(props);
        drawOrder.sort(Comparator.comparingDouble(this::depthKey));
        for(int i=drawOrder.size()-1;i>=0;i--){
            Prop p=drawOrder.get(i);
            PointF q=project(p.x,p.y,p.z);
            float hw=Math.max(10,p.type.vx*0.8f);
            float hh=Math.max(10,p.type.vz*0.8f);
            if(Math.abs(sx-q.x)<=hw&&Math.abs(sy-q.y)<=hh)return p;
        }
        return null;
    }

    private void startGrab(Prop p,float sx,float sy) {
        if(p==null||p.frozen)return;
        grabbed=p;mode=Mode.GRAB;
        PointF w=unprojectAtZ(sx,sy,p.z);
        grabOffsetX=w.x-p.x;grabOffsetY=w.y-p.y;
        targetX=p.x;targetY=p.y;targetZ=p.z;
    }

    private void updateGrabTarget(float sx,float sy) {
        if(grabbed==null)return;
        PointF w=unprojectAtZ(sx,sy,targetZ);
        targetX=clamp(w.x-grabOffsetX,-ROOM+grabbed.radius(),ROOM-grabbed.radius());
        targetY=clamp(w.y-grabOffsetY,-ROOM+grabbed.radius(),ROOM-grabbed.radius());
    }

    private void endGrab() {
        if(grabbed!=null){
            // Visual rotation is always one of four exact model rotations.
            grabbed.yaw=directionIndex(grabbed.yaw)*(float)(Math.PI/2.0);
        }
        grabbed=null;secondaryId=-1;
        if(mode==Mode.GRAB)mode=Mode.IDLE;
    }

    private void checkLongPress(long now) {
        if(pressed==null||primaryId<0||contextTriggered)return;
        if(mode!=Mode.GRAB&&mode!=Mode.IDLE)return;
        if(distance(downX,downY,primaryX,primaryY)<6f&&now-downNanos>550_000_000L){
            contextTriggered=true;contextProp=pressed;contextX=downX;contextY=downY;
            if(mode==Mode.GRAB)endGrab();
            primaryId=-1;mode=Mode.CONTEXT;feedback();
        }
    }

    @Override
    public boolean onTouchEvent(MotionEvent e) {
        int action=e.getActionMasked(),index=e.getActionIndex(),id=e.getPointerId(index);

        if(action==MotionEvent.ACTION_DOWN){
            if(!inside(e.getX(index),e.getY(index)))return false;
            float x=lx(e.getX(index)),y=ly(e.getY(index));
            if(handleUiDown(x,y))return true;

            primaryId=id;secondaryId=-1;primaryX=downX=x;primaryY=downY=y;
            downNanos=System.nanoTime();contextTriggered=false;
            pressed=pick(x,y);
            if(pressed!=null&&!pressed.frozen)startGrab(pressed,x,y);
            else mode=Mode.IDLE;
            return true;
        }

        if(action==MotionEvent.ACTION_POINTER_DOWN){
            if(mode==Mode.GRAB&&secondaryId<0){
                secondaryId=id;
                int pi=e.findPointerIndex(primaryId);
                if(pi>=0){
                    float x1=lx(e.getX(pi)),y1=ly(e.getY(pi));
                    float x2=lx(e.getX(index)),y2=ly(e.getY(index));
                    lastTwoDistance=distance(x1,y1,x2,y2);
                    lastTwoAngle=angle(x1,y1,x2,y2);
                }
            }
            return true;
        }

        if(action==MotionEvent.ACTION_MOVE){
            int pi=e.findPointerIndex(primaryId);
            if(pi>=0){
                primaryX=lx(e.getX(pi));primaryY=ly(e.getY(pi));
                if(mode==Mode.GRAB&&grabbed!=null)updateGrabTarget(primaryX,primaryY);
            }

            if(mode==Mode.GRAB&&grabbed!=null&&secondaryId>=0){
                int si=e.findPointerIndex(secondaryId);
                pi=e.findPointerIndex(primaryId);
                if(pi>=0&&si>=0){
                    float x1=lx(e.getX(pi)),y1=ly(e.getY(pi));
                    float x2=lx(e.getX(si)),y2=ly(e.getY(si));
                    float dist=distance(x1,y1,x2,y2);
                    targetZ=clamp(targetZ+(dist-lastTwoDistance)*0.0022f,grabbed.type.hMeters()*0.5f,0.85f);

                    float a=angle(x1,y1,x2,y2);
                    float da=wrapAngle(a-lastTwoAngle);
                    grabbed.yaw+=da;grabbed.spin+=da*2.0f;
                    lastTwoDistance=dist;lastTwoAngle=a;
                    updateGrabTarget(primaryX,primaryY);
                }
            }
            return true;
        }

        if(action==MotionEvent.ACTION_POINTER_UP){
            if(id==secondaryId){secondaryId=-1;return true;}
            if(id==primaryId){endGrab();pressed=null;primaryId=-1;return true;}
            return true;
        }

        if(action==MotionEvent.ACTION_UP||action==MotionEvent.ACTION_CANCEL){
            if(mode==Mode.GRAB)endGrab();
            pressed=null;primaryId=-1;secondaryId=-1;
            return true;
        }
        return true;
    }

    private boolean handleUiDown(float x,float y) {
        if(mode==Mode.SPAWN){
            float top=H-94f;
            if(y>=top){
                int col=clampInt((int)(x/(W/4f)),0,3);
                int row=clampInt((int)((y-top)/23f),0,3);
                int idx=row*4+col;
                if(idx>=0&&idx<PropType.values().length)spawnWithUndo(PropType.values()[idx]);
            }
            mode=Mode.IDLE;return true;
        }

        if(mode==Mode.SETTINGS){
            if(x<W-160||y<42||y>187){mode=Mode.IDLE;return true;}
            if(y>=73&&y<108){resetWorld();mode=Mode.IDLE;return true;}
            if(y>=108&&y<139){haptics=!haptics;feedback();saveWorld();return true;}
            return true;
        }

        if(mode==Mode.CONTEXT){
            float cx=clamp(contextX,68,W-68),cy=clamp(contextY,48,H-48);
            boolean left=x<cx,top=y<cy;
            Prop target=contextProp;mode=Mode.IDLE;contextProp=null;
            if(target!=null){
                if(top&&left)setFrozen(target,!target.frozen,true);
                else if(top)deleteWithUndo(target);
                else if(left)duplicateWithUndo(target);
                else cycleMaterial(target);
            }
            return true;
        }

        if(x<52&&y>H-42){mode=Mode.SPAWN;feedback();return true;}
        if(x>W-73&&y>H-42){doUndo();return true;}
        if(x>W-70&&y<42){mode=Mode.SETTINGS;feedback();return true;}
        return false;
    }

    private boolean inside(float sx,float sy) {
        return sx>=presentationX&&sx<=presentationX+presentationW&&sy>=presentationY&&sy<=presentationY+presentationH;
    }

    private float lx(float sx) {
        return clamp((sx-presentationX)/Math.max(0.0001f,presentationScale),0,W-1);
    }

    private float ly(float sy) {
        return clamp((sy-presentationY)/Math.max(0.0001f,presentationScale),0,H-1);
    }

    private void feedback() {
        if(haptics)performHapticFeedback(HapticFeedbackConstants.KEYBOARD_TAP);
    }

    private void saveWorld() {
        try{
            JSONArray arr=new JSONArray();
            for(Prop p:props){
                SaveState s=p.snapshot();
                JSONObject o=new JSONObject();
                o.put("id",s.id);o.put("type",s.type);o.put("material",s.material);
                o.put("x",s.x);o.put("y",s.y);o.put("z",s.z);o.put("yaw",s.yaw);o.put("frozen",s.frozen);
                arr.put(o);
            }
            prefs.edit().putString("world",arr.toString()).putBoolean("haptics",haptics).apply();
        }catch(Exception ignored){}
    }

    private boolean restoreWorld() {
        String data=prefs.getString("world","");
        if(data==null||data.isEmpty())return false;
        try{
            JSONArray arr=new JSONArray(data);
            if(arr.length()==0)return false;
            for(int i=0;i<arr.length();i++){
                JSONObject o=arr.getJSONObject(i);
                SaveState s=new SaveState();
                s.id=o.getInt("id");s.type=o.getString("type");s.material=o.getString("material");
                s.x=(float)o.getDouble("x");s.y=(float)o.getDouble("y");s.z=(float)o.getDouble("z");
                s.yaw=(float)o.getDouble("yaw");s.frozen=o.optBoolean("frozen",false);
                spawnFromState(s);
            }
            return !props.isEmpty();
        }catch(Exception e){
            props.clear();byId.clear();nextId=1;return false;
        }
    }

    private void spawnFromState(SaveState s) {
        spawnInternal(PropType.valueOf(s.type),MaterialKind.valueOf(s.material),s.x,s.y,s.z,s.yaw,s.frozen,s.id);
    }

    private void resetWorld() {
        endGrab();props.clear();byId.clear();undo.clear();nextId=1;
        pressed=null;contextProp=null;mode=Mode.IDLE;
        prefs.edit().remove("world").apply();
        createStarterSet();feedback();saveWorld();
    }

    private static int rgb(String hex) { return Color.parseColor("#"+hex); }
    private static float clamp(float v,float min,float max) { return Math.max(min,Math.min(max,v)); }
    private static int clampInt(int v,int min,int max) { return Math.max(min,Math.min(max,v)); }
    private static float distance(float x1,float y1,float x2,float y2) {
        float dx=x2-x1,dy=y2-y1;return (float)Math.sqrt(dx*dx+dy*dy);
    }
    private static float angle(float x1,float y1,float x2,float y2) {
        return (float)Math.atan2(y2-y1,x2-x1);
    }
    private static float wrapAngle(float a) {
        while(a>Math.PI)a-=(float)(Math.PI*2);
        while(a<-Math.PI)a+=(float)(Math.PI*2);
        return a;
    }
}
