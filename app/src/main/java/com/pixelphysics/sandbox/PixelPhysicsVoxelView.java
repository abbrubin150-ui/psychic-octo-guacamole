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
        float broadRadius() {
            if(circleFootprint())return radius();
            float hw=halfW(),hd=halfD();
            return (float)Math.sqrt(hw*hw+hd*hd);
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
        p.sleeping=frozen;
        p.sleepTimer=0f;
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
        p.sleeping=frozen;
        p.sleepTimer=0f;
        if(!frozen)p.wake();
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
        int substeps=adaptiveSubsteps(dt);
        float h=dt/substeps;

        for(int step=0;step<substeps;step++){
            integrateForces(h);
            integratePositions(h);

            // Positional constraint pass. Contacts are regenerated because the
            // previous correction changes the manifold for stacked bodies.
            for(int iter=0;iter<POSITION_ITERATIONS;iter++){
                collectContacts();
                if(contacts.isEmpty())break;
                sortContactsForSolver();
                for(int i=0;i<contacts.size();i++)solvePosition(contacts.get(i));
            }

            collectContacts();
            sortContactsForSolver();
            for(int iter=0;iter<VELOCITY_ITERATIONS;iter++){
                for(int i=0;i<contacts.size();i++)solveVelocity(contacts.get(i));
            }
        }

        updateSleeping(dt);

        for(int i=0;i<props.size();i++){
            Prop p=props.get(i);
            if(p.z<-1f){
                p.x=0f;p.y=0f;p.z=0.55f;
                p.vx=p.vy=p.vz=p.spin=0f;
                p.wake();
            }
        }
    }

    private int adaptiveSubsteps(float dt) {
        float worst=0f;
        float smallest=Float.POSITIVE_INFINITY;
        for(int i=0;i<props.size();i++){
            Prop p=props.get(i);
            if(p.frozen||p.sleeping)continue;
            float speed=(float)Math.sqrt(p.vx*p.vx+p.vy*p.vy+p.vz*p.vz);
            worst=Math.max(worst,speed);
            float feature=Math.min(p.type.wMeters(),Math.min(p.type.dMeters(),p.type.hMeters()));
            smallest=Math.min(smallest,Math.max(0.025f,feature));
        }
        if(smallest==Float.POSITIVE_INFINITY)return 1;
        float travel=worst*dt;
        int n=(int)Math.ceil(travel/Math.max(0.008f,smallest*0.35f));
        return clampInt(n,1,MAX_SUBSTEPS);
    }

    private void integrateForces(float dt) {
        for(int i=0;i<props.size();i++){
            Prop p=props.get(i);
            if(p.frozen||p.sleeping)continue;
            if(p!=grabbed)p.vz-=GRAVITY*dt;

            // Air drag is deliberately weak; contact friction does the real work.
            float linearDrag=(float)Math.exp(-0.08f*dt);
            float angularDrag=(float)Math.exp(-0.22f*dt);
            p.vx*=linearDrag;p.vy*=linearDrag;p.vz*=linearDrag;
            p.spin*=angularDrag;
        }
    }

    private void integratePositions(float dt) {
        for(int i=0;i<props.size();i++){
            Prop p=props.get(i);
            if(p.frozen||p.sleeping)continue;
            p.x+=p.vx*dt;
            p.y+=p.vy*dt;
            p.z+=p.vz*dt;
            p.yaw=wrapAngle(p.yaw+p.spin*dt);
        }
    }

    private void collectContacts() {
        contacts.clear();

        for(int i=0;i<props.size();i++){
            Prop p=props.get(i);
            addEnvironmentContacts(p);
        }

        for(int i=0;i<props.size();i++){
            Prop a=props.get(i);
            for(int j=i+1;j<props.size();j++){
                Prop b=props.get(j);
                if(a.frozen&&b.frozen)continue;
                addPairContact(a,b);
            }
        }
    }

    private void sortContactsForSolver() {
        contacts.sort((c1,c2)->{
            boolean vertical1=Math.abs(c1.nz)>0.5f;
            boolean vertical2=Math.abs(c2.nz)>0.5f;
            if(vertical1!=vertical2)return vertical1?-1:1;
            if(vertical1){
                int z=Float.compare(c1.pz,c2.pz);
                if(z!=0)return z;
            }
            int aId=c1.a==null?Integer.MAX_VALUE:c1.a.id;
            int bId=c2.a==null?Integer.MAX_VALUE:c2.a.id;
            return Integer.compare(aId,bId);
        });
    }

    private void addEnvironmentContacts(Prop p) {
        float envFriction=(float)Math.sqrt(friction(p.material)*0.72f);
        // The workshop floor/walls are treated as effectively rigid; the body's
        // restitution dominates the pair coefficient.
        float envRest=restitution(p.material);

        float bottom=p.bottom();
        if(bottom<=CONTACT_SLOP){
            Contact c=new Contact();
            c.a=p;c.b=null;
            c.nx=0;c.ny=0;c.nz=-1f; // normal from body toward floor
            c.px=p.x;c.py=p.y;c.pz=0f;
            c.penetration=Math.max(0f,-bottom);
            c.friction=p.type.renderKind==RenderKind.BALL
                    ? rollingResistance(p.material)
                    : envFriction;
            c.restitution=envRest;
            contacts.add(c);
        }

        float ex=footprintExtent(p,1f,0f);
        float ey=footprintExtent(p,0f,1f);

        float left=p.x-ex;
        if(left<=-ROOM+CONTACT_SLOP){
            float[] cp=supportPoint(p,-1f,0f);
            Contact c=environmentContact(p,-1f,0f,0f,
                    cp[0],cp[1],p.z,Math.max(0f,-ROOM-left),envFriction,envRest);
            contacts.add(c);
        }
        float right=p.x+ex;
        if(right>=ROOM-CONTACT_SLOP){
            float[] cp=supportPoint(p,1f,0f);
            Contact c=environmentContact(p,1f,0f,0f,
                    cp[0],cp[1],p.z,Math.max(0f,right-ROOM),envFriction,envRest);
            contacts.add(c);
        }
        float near=p.y-ey;
        if(near<=-ROOM+CONTACT_SLOP){
            float[] cp=supportPoint(p,0f,-1f);
            Contact c=environmentContact(p,0f,-1f,0f,
                    cp[0],cp[1],p.z,Math.max(0f,-ROOM-near),envFriction,envRest);
            contacts.add(c);
        }
        float far=p.y+ey;
        if(far>=ROOM-CONTACT_SLOP){
            float[] cp=supportPoint(p,0f,1f);
            Contact c=environmentContact(p,0f,1f,0f,
                    cp[0],cp[1],p.z,Math.max(0f,far-ROOM),envFriction,envRest);
            contacts.add(c);
        }
    }

    private Contact environmentContact(Prop p,float nx,float ny,float nz,
                                       float px,float py,float pz,float penetration,
                                       float friction,float restitution) {
        Contact c=new Contact();
        c.a=p;c.b=null;c.nx=nx;c.ny=ny;c.nz=nz;
        c.px=px;c.py=py;c.pz=pz;c.penetration=penetration;
        c.friction=friction;c.restitution=restitution;
        return c;
    }

    private void addPairContact(Prop a,Prop b) {
        float zReach=(a.type.hMeters()+b.type.hMeters())*0.5f+CONTACT_SLOP*2f;
        if(Math.abs(a.z-b.z)>zReach)return;

        float dxBroad=b.x-a.x,dyBroad=b.y-a.y;
        float rr=a.broadRadius()+b.broadRadius()+CONTACT_SLOP*2f;
        if(dxBroad*dxBroad+dyBroad*dyBroad>rr*rr)return;

        // Stair treads override the staircase's enclosing box when an object is
        // approaching from above.
        if(a.type==PropType.STAIRS || b.type==PropType.STAIRS){
            Hit2 stairFootprint=horizontalHit(a,b);
            if(stairFootprint.hit){
                Contact stair=stairContact(a,b,stairFootprint);
                if(stair!=null){
                    contacts.add(stair);
                    return;
                }
            }
        }

        // Hybrid non-voxel balls use true 3D sphere contacts.
        if(a.type.renderKind==RenderKind.BALL && b.type.renderKind==RenderKind.BALL){
            Contact sphere=sphereSphereContact(a,b);
            if(sphere!=null)contacts.add(sphere);
            return;
        }

        if(a.type.renderKind==RenderKind.BALL && !b.circleFootprint()){
            Contact sphere=sphereBoxContact(a,b);
            if(sphere!=null)contacts.add(sphere);
            return;
        }
        if(b.type.renderKind==RenderKind.BALL && !a.circleFootprint()){
            Contact sphere=sphereBoxContact(b,a);
            if(sphere!=null)contacts.add(sphere);
            return;
        }

        Hit2 h=horizontalHit(a,b);
        if(!h.hit)return;

        float zPen=Math.min(a.top(),b.top())-Math.max(a.bottom(),b.bottom());
        if(zPen<-CONTACT_SLOP)return;

        boolean verticalAxis=zPen<=h.penetration+CONTACT_SLOP;
        Contact contact=new Contact();
        contact.a=a;contact.b=b;
        contact.friction=(float)Math.sqrt(friction(a.material)*friction(b.material));
        contact.restitution=(float)Math.sqrt(restitution(a.material)*restitution(b.material));

        if(verticalAxis){
            boolean bAbove=b.z>=a.z;
            contact.nx=0f;contact.ny=0f;contact.nz=bAbove?1f:-1f;
            contact.px=(a.x+b.x)*0.5f;
            contact.py=(a.y+b.y)*0.5f;
            contact.pz=bAbove?Math.min(a.top(),b.bottom()):Math.max(a.bottom(),b.top());
            contact.penetration=Math.max(0f,zPen);
        }else{
            contact.nx=h.nx;contact.ny=h.ny;contact.nz=0f;
            contact.px=h.px;contact.py=h.py;
            contact.pz=(Math.max(a.bottom(),b.bottom())+Math.min(a.top(),b.top()))*0.5f;
            contact.penetration=Math.max(0f,h.penetration);
        }
        contacts.add(contact);
    }

    private Contact sphereSphereContact(Prop a,Prop b) {
        float dx=b.x-a.x,dy=b.y-a.y,dz=b.z-a.z;
        float ra=a.type.wMeters()*0.5f,rb=b.type.wMeters()*0.5f;
        float sum=ra+rb;
        float d2=dx*dx+dy*dy+dz*dz;
        float limit=sum+CONTACT_SLOP;
        if(d2>limit*limit)return null;

        float d=(float)Math.sqrt(Math.max(d2,1e-12f));
        float nx,ny,nz;
        if(d<0.000001f){nx=1f;ny=0f;nz=0f;}
        else{nx=dx/d;ny=dy/d;nz=dz/d;}

        Contact contact=new Contact();
        contact.a=a;contact.b=b;
        contact.nx=nx;contact.ny=ny;contact.nz=nz;
        contact.penetration=Math.max(0f,sum-d);
        float along=ra-contact.penetration*0.5f;
        contact.px=a.x+nx*along;
        contact.py=a.y+ny*along;
        contact.pz=a.z+nz*along;
        contact.friction=(float)Math.sqrt(friction(a.material)*friction(b.material));
        contact.restitution=(float)Math.sqrt(restitution(a.material)*restitution(b.material));
        return contact;
    }

    private Contact sphereBoxContact(Prop sphere,Prop box) {
        float radius=sphere.type.wMeters()*0.5f;
        float dx=sphere.x-box.x,dy=sphere.y-box.y,dz=sphere.z-box.z;
        float cs=(float)Math.cos(box.yaw),sn=(float)Math.sin(box.yaw);

        float lx= cs*dx+sn*dy;
        float ly=-sn*dx+cs*dy;
        float lz=dz;

        float hx=box.halfW(),hy=box.halfD(),hz=box.type.hMeters()*0.5f;
        float qx=clamp(lx,-hx,hx);
        float qy=clamp(ly,-hy,hy);
        float qz=clamp(lz,-hz,hz);

        float ex=qx-lx,ey=qy-ly,ez=qz-lz; // sphere -> closest box point
        float d2=ex*ex+ey*ey+ez*ez;

        float nlx,nly,nlz,penetration;
        if(d2>1e-12f){
            float d=(float)Math.sqrt(d2);
            if(d>radius+CONTACT_SLOP)return null;
            nlx=ex/d;nly=ey/d;nlz=ez/d;
            penetration=Math.max(0f,radius-d);
        }else{
            float toX=hx-Math.abs(lx);
            float toY=hy-Math.abs(ly);
            float toZ=hz-Math.abs(lz);
            if(toX<=toY && toX<=toZ){
                float outward=lx>=0?1f:-1f;
                nlx=-outward;nly=0f;nlz=0f;
                qx=outward*hx;qy=ly;qz=lz;
                penetration=radius+toX;
            }else if(toY<=toZ){
                float outward=ly>=0?1f:-1f;
                nlx=0f;nly=-outward;nlz=0f;
                qx=lx;qy=outward*hy;qz=lz;
                penetration=radius+toY;
            }else{
                float outward=lz>=0?1f:-1f;
                nlx=0f;nly=0f;nlz=-outward;
                qx=lx;qy=ly;qz=outward*hz;
                penetration=radius+toZ;
            }
        }

        float nx=cs*nlx-sn*nly;
        float ny=sn*nlx+cs*nly;
        float nz=nlz;

        Contact contact=new Contact();
        contact.a=sphere;contact.b=box;
        contact.nx=nx;contact.ny=ny;contact.nz=nz;
        contact.penetration=penetration;
        contact.px=box.x+cs*qx-sn*qy;
        contact.py=box.y+sn*qx+cs*qy;
        contact.pz=box.z+qz;
        contact.friction=(float)Math.sqrt(friction(sphere.material)*friction(box.material));
        contact.restitution=(float)Math.sqrt(restitution(sphere.material)*restitution(box.material));
        return contact;
    }

    private Contact stairContact(Prop a,Prop b,Hit2 h) {
        if(a.type==PropType.STAIRS && b.z>=a.z){
            float top=stairTopAt(a,b.x,b.y,b.radius()*0.20f);
            if(!Float.isNaN(top) && b.bottom()<=top+CONTACT_SLOP && b.bottom()>=a.bottom()-0.02f){
                Contact c=new Contact();
                c.a=a;c.b=b;c.nx=0;c.ny=0;c.nz=1f;
                c.px=b.x;c.py=b.y;c.pz=top;
                c.penetration=Math.max(0f,top-b.bottom());
                c.friction=(float)Math.sqrt(friction(a.material)*friction(b.material));
                c.restitution=(float)Math.sqrt(restitution(a.material)*restitution(b.material));
                return c;
            }
        }
        if(b.type==PropType.STAIRS && a.z>=b.z){
            float top=stairTopAt(b,a.x,a.y,a.radius()*0.20f);
            if(!Float.isNaN(top) && a.bottom()<=top+CONTACT_SLOP && a.bottom()>=b.bottom()-0.02f){
                Contact c=new Contact();
                c.a=a;c.b=b;c.nx=0;c.ny=0;c.nz=-1f;
                c.px=a.x;c.py=a.y;c.pz=top;
                c.penetration=Math.max(0f,top-a.bottom());
                c.friction=(float)Math.sqrt(friction(a.material)*friction(b.material));
                c.restitution=(float)Math.sqrt(restitution(a.material)*restitution(b.material));
                return c;
            }
        }
        return null;
    }

    private float stairTopAt(Prop stairs,float wx,float wy,float margin) {
        float dx=wx-stairs.x,dy=wy-stairs.y;
        float cs=(float)Math.cos(stairs.yaw),sn=(float)Math.sin(stairs.yaw);
        float lx= cs*dx+sn*dy;
        float ly=-sn*dx+cs*dy;
        float hx=stairs.halfW(),hy=stairs.halfD();
        if(lx<-hx-margin||lx>hx+margin||ly<-hy-margin||ly>hy+margin)return Float.NaN;
        int steps=8;
        float u=clamp((lx+hx)/(2f*hx),0f,0.9999f);
        int index=clampInt((int)(u*steps),0,steps-1);
        return stairs.bottom()+stairs.type.hMeters()*(index+1)/(float)steps;
    }

    private Hit2 horizontalHit(Prop a,Prop b) {
        if(a.circleFootprint() && b.circleFootprint())return circleCircle(a,b);
        if(a.circleFootprint() && !b.circleFootprint())return circleBox(a,b);
        if(!a.circleFootprint() && b.circleFootprint()){
            Hit2 h=circleBox(b,a);
            if(h.hit){h.nx=-h.nx;h.ny=-h.ny;}
            return h;
        }
        return boxBox(a,b);
    }

    private Hit2 circleCircle(Prop a,Prop b) {
        Hit2 h=new Hit2();
        float dx=b.x-a.x,dy=b.y-a.y;
        float r=a.radius()+b.radius();
        float d2=dx*dx+dy*dy;
        if(d2>=r*r)return h;
        float d=(float)Math.sqrt(Math.max(d2,1e-10f));
        h.hit=true;
        if(d<0.00001f){h.nx=1f;h.ny=0f;}
        else{h.nx=dx/d;h.ny=dy/d;}
        h.penetration=r-d;
        h.px=a.x+h.nx*(a.radius()-h.penetration*0.5f);
        h.py=a.y+h.ny*(a.radius()-h.penetration*0.5f);
        return h;
    }

    private Hit2 circleBox(Prop circle,Prop box) {
        Hit2 h=new Hit2();
        float dx=circle.x-box.x,dy=circle.y-box.y;
        float cs=(float)Math.cos(box.yaw),sn=(float)Math.sin(box.yaw);
        float lx= cs*dx+sn*dy;
        float ly=-sn*dx+cs*dy;
        float hx=box.halfW(),hy=box.halfD();
        float qx=clamp(lx,-hx,hx);
        float qy=clamp(ly,-hy,hy);
        float ex=qx-lx,ey=qy-ly;
        float d2=ex*ex+ey*ey;
        float r=circle.radius();

        float nlx,nly,penetration;
        if(d2>1e-10f){
            float d=(float)Math.sqrt(d2);
            if(d>=r)return h;
            nlx=ex/d;nly=ey/d; // circle -> box
            penetration=r-d;
        }else{
            float toX=hx-Math.abs(lx);
            float toY=hy-Math.abs(ly);
            if(toX<toY){
                float outward=lx>=0?1f:-1f;
                nlx=-outward;nly=0f;
                penetration=r+toX;
                qx=outward*hx;qy=ly;
            }else{
                float outward=ly>=0?1f:-1f;
                nlx=0f;nly=-outward;
                penetration=r+toY;
                qx=lx;qy=outward*hy;
            }
        }

        h.hit=true;
        h.nx=cs*nlx-sn*nly;
        h.ny=sn*nlx+cs*nly;
        h.penetration=penetration;
        h.px=box.x+cs*qx-sn*qy;
        h.py=box.y+sn*qx+cs*qy;
        return h;
    }

    private Hit2 boxBox(Prop a,Prop b) {
        Hit2 h=new Hit2();
        float ca=(float)Math.cos(a.yaw),sa=(float)Math.sin(a.yaw);
        float cb=(float)Math.cos(b.yaw),sb=(float)Math.sin(b.yaw);
        float[][] axes={{ca,sa},{-sa,ca},{cb,sb},{-sb,cb}};
        float dx=b.x-a.x,dy=b.y-a.y;
        float best=Float.POSITIVE_INFINITY,bnx=0,bny=0;

        for(int i=0;i<4;i++){
            float ax=axes[i][0],ay=axes[i][1];
            float ra=footprintExtent(a,ax,ay);
            float rb=footprintExtent(b,ax,ay);
            float dist=dx*ax+dy*ay;
            float overlap=ra+rb-Math.abs(dist);
            if(overlap<=0f)return h;
            if(overlap<best){
                best=overlap;
                float sign=dist>=0?1f:-1f;
                bnx=ax*sign;bny=ay*sign;
            }
        }

        float[] saPoint=supportPoint(a,bnx,bny);
        float[] sbPoint=supportPoint(b,-bnx,-bny);
        h.hit=true;h.nx=bnx;h.ny=bny;h.penetration=best;
        h.px=(saPoint[0]+sbPoint[0])*0.5f;
        h.py=(saPoint[1]+sbPoint[1])*0.5f;
        return h;
    }

    private float footprintExtent(Prop p,float ax,float ay) {
        if(p.circleFootprint())return p.radius();
        float cs=(float)Math.cos(p.yaw),sn=(float)Math.sin(p.yaw);
        float ux=cs,uy=sn;
        float vx=-sn,vy=cs;
        return p.halfW()*Math.abs(ux*ax+uy*ay)+p.halfD()*Math.abs(vx*ax+vy*ay);
    }

    private float[] supportPoint(Prop p,float dx,float dy) {
        if(p.circleFootprint()){
            float len=(float)Math.sqrt(Math.max(dx*dx+dy*dy,1e-10f));
            return new float[]{p.x+dx/len*p.radius(),p.y+dy/len*p.radius()};
        }
        float cs=(float)Math.cos(p.yaw),sn=(float)Math.sin(p.yaw);
        float ux=cs,uy=sn,vx=-sn,vy=cs;
        float su=(ux*dx+uy*dy)>=0?1f:-1f;
        float sv=(vx*dx+vy*dy)>=0?1f:-1f;
        return new float[]{
                p.x+ux*p.halfW()*su+vx*p.halfD()*sv,
                p.y+uy*p.halfW()*su+vy*p.halfD()*sv
        };
    }

    private void solveVelocity(Contact c) {
        Prop a=c.a,b=c.b;
        float ia=a.invMass(),ib=b==null?0f:b.invMass();
        float iia=a.invInertia(),iib=b==null?0f:b.invInertia();

        float rax=c.px-a.x,ray=c.py-a.y;
        float rbx=b==null?0f:c.px-b.x,rby=b==null?0f:c.py-b.y;

        float vax=a.vx-a.spin*ray;
        float vay=a.vy+a.spin*rax;
        float vaz=a.vz;
        float vbx=b==null?0f:b.vx-b.spin*rby;
        float vby=b==null?0f:b.vy+b.spin*rbx;
        float vbz=b==null?0f:b.vz;

        float rvx=vbx-vax,rvy=vby-vay,rvz=vbz-vaz;
        float vn=rvx*c.nx+rvy*c.ny+rvz*c.nz;
        if(vn>=0f)return;

        float ran=rax*c.ny-ray*c.nx;
        float rbn=rbx*c.ny-rby*c.nx;
        float k=ia+ib+ran*ran*iia+rbn*rbn*iib;
        if(k<1e-8f)return;

        float impactSpeed=-vn;
        float e=impactSpeed>RESTITUTION_VELOCITY_THRESHOLD?c.restitution:0f;
        float j=-(1f+e)*vn/k;
        applyImpulse(a,-j*c.nx,-j*c.ny,-j*c.nz,rax,ray);
        if(b!=null)applyImpulse(b,j*c.nx,j*c.ny,j*c.nz,rbx,rby);

        if(j>0.015f && impactSpeed>0.25f){
            a.wake();
            if(b!=null)b.wake();
        }

        // Recompute relative velocity after the normal impulse.
        vax=a.vx-a.spin*ray;vay=a.vy+a.spin*rax;vaz=a.vz;
        vbx=b==null?0f:b.vx-b.spin*rby;
        vby=b==null?0f:b.vy+b.spin*rbx;
        vbz=b==null?0f:b.vz;
        rvx=vbx-vax;rvy=vby-vay;rvz=vbz-vaz;

        float vtX=rvx-c.nx*(rvx*c.nx+rvy*c.ny+rvz*c.nz);
        float vtY=rvy-c.ny*(rvx*c.nx+rvy*c.ny+rvz*c.nz);
        float vtZ=rvz-c.nz*(rvx*c.nx+rvy*c.ny+rvz*c.nz);
        float vtLen=(float)Math.sqrt(vtX*vtX+vtY*vtY+vtZ*vtZ);
        if(vtLen<1e-6f)return;
        float tx=vtX/vtLen,ty=vtY/vtLen,tz=vtZ/vtLen;

        float rat=rax*ty-ray*tx;
        float rbt=rbx*ty-rby*tx;
        float kt=ia+ib+rat*rat*iia+rbt*rbt*iib;
        if(kt<1e-8f)return;
        float jt=-(rvx*tx+rvy*ty+rvz*tz)/kt;
        float dynamicLimit=c.friction*j;
        float staticLimit=Math.min(1.35f,c.friction*1.24f)*j;
        if(Math.abs(jt)>staticLimit)jt=Math.copySign(dynamicLimit,jt);
        else jt=clamp(jt,-staticLimit,staticLimit);

        applyImpulse(a,-jt*tx,-jt*ty,-jt*tz,rax,ray);
        if(b!=null)applyImpulse(b,jt*tx,jt*ty,jt*tz,rbx,rby);

        // A single-point 2.5D floor contact otherwise has no torsional friction.
        // Approximate the distributed contact patch with a bounded angular impulse.
        if(b==null && Math.abs(c.nz)>0.5f && Math.abs(a.spin)>0.0001f){
            float invI=a.invInertia();
            if(invI>1e-8f){
                float desired=-a.spin/invI;
                float limit=c.friction*j*Math.max(0.01f,a.radius()*0.35f);
                float angularImpulse=clamp(desired,-limit,limit);
                a.spin+=angularImpulse*invI;
            }
        }
    }

    private void applyImpulse(Prop p,float ix,float iy,float iz,float rx,float ry) {
        if(p==null||p.frozen)return;
        p.vx+=ix*p.invMass();
        p.vy+=iy*p.invMass();
        p.vz+=iz*p.invMass();
        p.spin+=(rx*iy-ry*ix)*p.invInertia();
    }

    private void solvePosition(Contact c) {
        float penetration=Math.max(0f,c.penetration-CONTACT_SLOP);
        if(penetration<=0f)return;
        Prop a=c.a,b=c.b;
        float ia=a.invMass(),ib=b==null?0f:b.invMass();
        float sum=ia+ib;
        if(sum<=1e-8f)return;
        float correction=POSITION_BETA*penetration/sum;

        if(!a.frozen){
            a.x-=c.nx*correction*ia;
            a.y-=c.ny*correction*ia;
            a.z-=c.nz*correction*ia;
        }
        if(b!=null&&!b.frozen){
            b.x+=c.nx*correction*ib;
            b.y+=c.ny*correction*ib;
            b.z+=c.nz*correction*ib;
        }
    }

    private boolean supported(Prop p) {
        if(p.bottom()<=CONTACT_SLOP*2f)return true;
        for(int i=0;i<props.size();i++){
            Prop q=props.get(i);
            if(q==p)continue;
            float dx=q.x-p.x,dy=q.y-p.y;
            float rr=p.broadRadius()+q.broadRadius()+CONTACT_SLOP*4f;
            if(dx*dx+dy*dy>rr*rr)continue;
            Hit2 h=horizontalHit(p,q);
            if(!h.hit)continue;
            float top=q.type==PropType.STAIRS?stairTopAt(q,p.x,p.y,p.radius()*0.15f):q.top();
            if(Float.isNaN(top))continue;
            if(Math.abs(p.bottom()-top)<=CONTACT_SLOP*4f)return true;
        }
        return false;
    }

    private void updateSleeping(float dt) {
        for(int i=0;i<props.size();i++){
            Prop p=props.get(i);
            if(p.frozen||p==grabbed){
                p.sleepTimer=0f;
                if(p==grabbed)p.sleeping=false;
                continue;
            }

            float horizontal=p.vx*p.vx+p.vy*p.vy;
            boolean slow=horizontal<SLEEP_LINEAR*SLEEP_LINEAR
                    && Math.abs(p.vz)<SLEEP_VERTICAL
                    && Math.abs(p.spin)<SLEEP_ANGULAR
                    && supported(p);

            if(slow){
                p.sleepTimer+=dt;
                if(p.sleepTimer>=SLEEP_TIME){
                    p.sleeping=true;
                    p.vx=p.vy=p.vz=p.spin=0f;
                }
            }else{
                p.sleepTimer=0f;
                p.sleeping=false;
            }
        }
    }

    private float restitution(MaterialKind m) {
        return m==MaterialKind.RUBBER?0.72f:(m==MaterialKind.METAL?0.08f:0.14f);
    }

    private float friction(MaterialKind m) {
        return m==MaterialKind.RUBBER?0.96f:(m==MaterialKind.METAL?0.36f:0.64f);
    }

    private float rollingResistance(MaterialKind m) {
        // Effective rolling coefficient for the reduced-DOF 2.5D ball model.
        return m==MaterialKind.RUBBER?0.055f:(m==MaterialKind.METAL?0.018f:0.032f);
    }

    private void updateGrab(float dt) {
        if(mode!=Mode.GRAB||grabbed==null||grabbed.frozen)return;
        grabbed.wake();

        float cs=(float)Math.cos(grabbed.yaw),sn=(float)Math.sin(grabbed.yaw);
        float rx=cs*grabLocalX-sn*grabLocalY;
        float ry=sn*grabLocalX+cs*grabLocalY;
        float rz=grabLocalZ;

        float grabX=grabbed.x+rx;
        float grabY=grabbed.y+ry;
        float grabZ=grabbed.z+rz;

        float pointVx=grabbed.vx-grabbed.spin*ry;
        float pointVy=grabbed.vy+grabbed.spin*rx;
        float pointVz=grabbed.vz;

        float mass=Math.max(0.25f,grabbed.type.mass);
        float ms=(float)Math.pow(mass,0.34f);
        float kp=64f*ms;
        float kd=2f*(float)Math.sqrt(kp*mass)*0.82f;

        float fx=(targetX-grabX)*kp-pointVx*kd;
        float fy=(targetY-grabY)*kp-pointVy*kd;
        float fz=(targetZ-grabZ)*kp-pointVz*kd;
        float cap=120f*(float)Math.pow(Math.max(1f,mass),0.58f);
        float len=(float)Math.sqrt(fx*fx+fy*fy+fz*fz);
        if(len>cap){
            float s=cap/len;fx*=s;fy*=s;fz*=s;
        }

        grabbed.vx+=fx*grabbed.invMass()*dt;
        grabbed.vy+=fy*grabbed.invMass()*dt;
        grabbed.vz+=fz*grabbed.invMass()*dt;

        float torque=rx*fy-ry*fx;

        // The gesture controls an angular spring, not the angle directly.
        float angleError=wrapAngle(targetYaw-grabbed.yaw);
        float inertia=grabbed.inertia();
        float angularAccel=42f*angleError-11f*grabbed.spin;
        torque+=inertia*angularAccel;
        float maxTorque=Math.max(0.02f,inertia*85f);
        torque=clamp(torque,-maxTorque,maxTorque);
        grabbed.spin+=torque*grabbed.invInertia()*dt;
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
        p.wake();

        PointF w=unprojectAtZ(sx,sy,p.z);
        float dx=w.x-p.x,dy=w.y-p.y;
        float cs=(float)Math.cos(p.yaw),sn=(float)Math.sin(p.yaw);
        grabLocalX= cs*dx+sn*dy;
        grabLocalY=-sn*dx+cs*dy;
        grabLocalZ=0f;

        targetX=w.x;
        targetY=w.y;
        targetZ=p.z;
        targetYaw=p.yaw;
    }

    private void updateGrabTarget(float sx,float sy) {
        if(grabbed==null)return;
        PointF w=unprojectAtZ(sx,sy,targetZ);
        float margin=Math.max(0.02f,grabbed.radius()*0.35f);
        targetX=clamp(w.x,-ROOM+margin,ROOM-margin);
        targetY=clamp(w.y,-ROOM+margin,ROOM-margin);
    }

    private void endGrab() {
        if(grabbed!=null)grabbed.wake();
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
                    targetYaw=wrapAngle(targetYaw+da);
                    grabbed.wake();
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
