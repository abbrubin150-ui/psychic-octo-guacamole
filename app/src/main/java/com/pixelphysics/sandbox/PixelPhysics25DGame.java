package com.pixelphysics.sandbox;

import com.badlogic.gdx.*;
import com.badlogic.gdx.graphics.*;
import com.badlogic.gdx.graphics.g2d.*;
import com.badlogic.gdx.graphics.glutils.FrameBuffer;
import com.badlogic.gdx.graphics.glutils.ShapeRenderer;
import com.badlogic.gdx.math.*;
import com.badlogic.gdx.utils.*;

import java.util.Comparator;
import java.util.Locale;

public class PixelPhysics25DGame extends ApplicationAdapter implements InputProcessor {
    private static final int W = 480;
    private static final int H = 270;
    private static final int FRAME = 80;
    private static final float ROOM = 4.3f;
    private static final float WALL_H = 3.6f;
    private static final float ORIGIN_X = 240f;
    private static final float ORIGIN_Y = 88f;
    private static final float ISO_X = 26f;
    private static final float ISO_Y = 10f;
    private static final float Z_PX = 24f;
    private static final float FIXED_DT = 1f / 60f;
    private static final float MAX_ACCUM = 0.12f;
    private static final float GRAVITY = 9.81f;
    private static final int MAX_PROPS = 100;

    private SpriteBatch batch;
    private ShapeRenderer shapes;
    private BitmapFont font;
    private FrameBuffer buffer;
    private TextureRegion bufferRegion;
    private final Matrix4 projection = new Matrix4();

    private float presentScale = 1f;
    private float presentX, presentY, presentW = W, presentH = H;

    private final Array<Prop> props = new Array<>();
    private final Array<Prop> drawOrder = new Array<>();
    private final ObjectMap<Integer, Prop> byId = new ObjectMap<>();
    private final ObjectMap<String, Texture[]> sprites = new ObjectMap<>();
    private final Array<UndoAction> undo = new Array<>();
    private final Array<Platform> platforms = new Array<>();
    private final Array<RampSurface> ramps = new Array<>();
    private int nextId = 1;

    private Preferences prefs;
    private Json json;
    private float accumulator;
    private float autosaveClock;
    private boolean haptics = true;

    private enum Mode { IDLE, GRAB, CONTEXT, SPAWN, SETTINGS }
    private Mode mode = Mode.IDLE;

    private final int[] pointerX = new int[10];
    private final int[] pointerY = new int[10];
    private int primaryPointer = -1;
    private int secondPointer = -1;
    private float downX, downY;
    private long downNanos;
    private Prop pressed;
    private Prop grabbed;
    private Prop contextProp;
    private boolean contextTriggered;
    private float contextX, contextY;

    private float grabOffsetX, grabOffsetY;
    private float targetX, targetY, targetZ;
    private float lastTwoDistance, lastTwoAngle;

    private final Vector2 tmp2 = new Vector2();
    private final Vector2 tmp3 = new Vector2();
    private final Vector2 tmp4 = new Vector2();

    private enum MaterialKind { MAHOGANY, METAL, RUBBER }

    private enum PropType {
        CUBE("CUBE", 1.0f, 1.0f, 1.0f, 1.0f, MaterialKind.MAHOGANY, false),
        BEAM_SHORT("BEAM S", 2.2f, 0.55f, 0.50f, 1.3f, MaterialKind.MAHOGANY, false),
        BEAM_LONG("BEAM L", 3.3f, 0.55f, 0.50f, 2.0f, MaterialKind.MAHOGANY, false),
        PLANK("PLANK", 2.7f, 1.0f, 0.28f, 1.2f, MaterialKind.MAHOGANY, false),
        WOOD_BALL("WOOD BALL", 0.9f, 0.9f, 0.9f, 0.7f, MaterialKind.MAHOGANY, true),
        METAL_BALL("METAL BALL", 0.9f, 0.9f, 0.9f, 4.0f, MaterialKind.METAL, true),
        WEIGHT("WEIGHT", 1.0f, 1.0f, 1.05f, 10.0f, MaterialKind.METAL, false),
        WHEEL("WHEEL", 1.25f, 0.42f, 1.25f, 1.1f, MaterialKind.MAHOGANY, false),
        RUBBER_BALL("RUBBER", 1.0f, 1.0f, 1.0f, 0.85f, MaterialKind.RUBBER, true),
        BARREL("BARREL", 1.1f, 1.1f, 1.5f, 3.6f, MaterialKind.METAL, false),
        CRATE("CRATE", 1.4f, 1.4f, 1.4f, 1.8f, MaterialKind.MAHOGANY, false),
        RAMP("RAMP", 2.4f, 1.4f, 0.8f, 1.8f, MaterialKind.MAHOGANY, false),
        SPRING("SPRING", 0.8f, 0.8f, 2.0f, 1.5f, MaterialKind.METAL, false);

        final String label;
        final float w, d, h, mass;
        final MaterialKind defaultMaterial;
        final boolean sphere;

        PropType(String label, float w, float d, float h, float mass, MaterialKind defaultMaterial, boolean sphere) {
            this.label = label;
            this.w = w;
            this.d = d;
            this.h = h;
            this.mass = mass;
            this.defaultMaterial = defaultMaterial;
            this.sphere = sphere;
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

    private class Prop {
        int id;
        PropType type;
        MaterialKind material;
        float x, y, z;
        float vx, vy, vz;
        float yaw, spin;
        boolean frozen;

        float radius() { return Math.max(type.w, type.d) * 0.5f; }
        float bottom() { return z - type.h * 0.5f; }
        float top() { return z + type.h * 0.5f; }

        SaveState snapshot() {
            SaveState s = new SaveState();
            s.id = id;
            s.type = type.name();
            s.material = material.name();
            s.x = x; s.y = y; s.z = z; s.yaw = yaw; s.frozen = frozen;
            return s;
        }
    }

    private static class SaveState {
        public int id;
        public String type;
        public String material;
        public float x, y, z, yaw;
        public boolean frozen;
        public SaveState() {}
    }

    private static class Platform {
        float cx, cy, w, d, top;
        Platform(float cx, float cy, float w, float d, float top) {
            this.cx = cx; this.cy = cy; this.w = w; this.d = d; this.top = top;
        }
        boolean contains(float x, float y, float margin) {
            return Math.abs(x-cx) <= w*0.5f + margin && Math.abs(y-cy) <= d*0.5f + margin;
        }
    }

    private static class RampSurface {
        float x0, x1, y0, y1, z0, z1;
        RampSurface(float x0, float x1, float y0, float y1, float z0, float z1) {
            this.x0=x0; this.x1=x1; this.y0=y0; this.y1=y1; this.z0=z0; this.z1=z1;
        }
        boolean contains(float x,float y,float margin) {
            return x>=Math.min(x0,x1)-margin && x<=Math.max(x0,x1)+margin &&
                    y>=Math.min(y0,y1)-margin && y<=Math.max(y0,y1)+margin;
        }
        float heightAt(float x) {
            float t = MathUtils.clamp((x-x0)/(x1-x0),0f,1f);
            return MathUtils.lerp(z0,z1,t);
        }
        float downhillSign() { return z1>z0 ? -1f : 1f; }
    }

    private interface UndoAction { void undo(); }

    @Override
    public void create() {
        Locale.setDefault(Locale.US);
        prefs = Gdx.app.getPreferences("pixel-physics-25d-world-v1");
        json = new Json();

        batch = new SpriteBatch();
        shapes = new ShapeRenderer();
        font = new BitmapFont();
        font.getData().setScale(0.70f);
        font.getRegion().getTexture().setFilter(Texture.TextureFilter.Nearest, Texture.TextureFilter.Nearest);
        projection.setToOrtho2D(0,0,W,H);

        createSpriteSets();
        buildStaticPhysics();
        recreateBuffer();
        if (!restoreWorld()) createStarterSet();

        Gdx.input.setInputProcessor(this);
    }

    private Palette palette(MaterialKind k) {
        switch(k) {
            case METAL: return new Palette("101319","272D34","4A5660","73818B","AAB4BA","E0E4E6");
            case RUBBER: return new Palette("0F1511","1A241C","2B3B2E","435746","627A63","A5B39D");
            default: return new Palette("24100D","421813","6F271A","963B21","C25D30","E88A49");
        }
    }

    private void createSpriteSets() {
        for (PropType type : PropType.values()) {
            for (MaterialKind material : MaterialKind.values()) {
                Texture[] frames = new Texture[4];
                for (int dir=0; dir<4; dir++) frames[dir] = makeSprite(type,material,dir);
                sprites.put(spriteKey(type,material),frames);
            }
        }
    }

    private String spriteKey(PropType type, MaterialKind material) {
        return type.name()+":"+material.name();
    }

    private Texture makeSprite(PropType type, MaterialKind material, int dir) {
        Pixmap p = new Pixmap(FRAME,FRAME,Pixmap.Format.RGBA8888);
        p.setBlending(Pixmap.Blending.None);
        p.setColor(0,0,0,0); p.fill();
        p.setBlending(Pixmap.Blending.SourceOver);
        Palette pal = palette(material);

        if (type == PropType.WOOD_BALL || type == PropType.METAL_BALL || type == PropType.RUBBER_BALL) {
            drawBallSprite(p,type,pal,dir);
        } else if (type == PropType.WHEEL) {
            drawWheelSprite(p,pal,dir);
        } else if (type == PropType.RAMP) {
            drawRampSprite(p,pal,dir);
        } else if (type == PropType.SPRING) {
            drawSpringSprite(p,pal,dir);
        } else if (type == PropType.BARREL) {
            drawBarrelSprite(p,pal,dir);
        } else {
            drawCuboidSprite(p,type,pal,dir);
        }

        Texture t = new Texture(p);
        p.dispose();
        t.setFilter(Texture.TextureFilter.Nearest,Texture.TextureFilter.Nearest);
        return t;
    }

    private void drawCuboidSprite(Pixmap p, PropType type, Palette pal, int dir) {
        float rw = (dir%2==0?type.w:type.d);
        float rd = (dir%2==0?type.d:type.w);
        int hx = Math.max(5,Math.round(rw*8.5f));
        int hy = Math.max(3,Math.round(rd*4.0f));
        int hh = Math.max(5,Math.round(type.h*15f));
        hx = Math.min(28,hx); hy = Math.min(14,hy); hh = Math.min(34,hh);

        int cx=FRAME/2, cy=18+hy;
        int ax=cx, ay=cy-hy;
        int bx=cx+hx, by=cy;
        int cxp=cx, cyp=cy+hy;
        int dx=cx-hx, dy=cy;
        int b2x=bx,b2y=by+hh,c2x=cxp,c2y=cyp+hh,d2x=dx,d2y=dy+hh;

        Color top = dir==0?pal.light:(dir==1?pal.mid:(dir==2?pal.base:pal.hi));
        Color left = dir==3?pal.light:pal.base;
        Color right = dir==1?pal.light:pal.dark;

        p.setColor(left);
        fillQuad(p,dx,dy,cxp,cyp,c2x,c2y,d2x,d2y);
        p.setColor(right);
        fillQuad(p,bx,by,cxp,cyp,c2x,c2y,b2x,b2y);
        p.setColor(top);
        fillQuad(p,ax,ay,bx,by,cxp,cyp,dx,dy);

        p.setColor(pal.outline);
        lineLoop(p,new int[]{ax,ay,bx,by,cxp,cyp,dx,dy});
        p.drawLine(dx,dy,d2x,d2y); p.drawLine(cxp,cyp,c2x,c2y); p.drawLine(bx,by,b2x,b2y);
        p.drawLine(d2x,d2y,c2x,c2y); p.drawLine(c2x,c2y,b2x,b2y);

        p.setColor(pal.mid);
        if(type==PropType.CRATE) {
            p.drawLine(dx+4,dy+5,c2x-4,c2y-5);
            p.drawLine(d2x+4,d2y-5,cxp-4,cyp+5);
            p.drawLine(bx-4,by+5,c2x+4,c2y-5);
            p.drawLine(b2x-4,b2y-5,cxp+4,cyp+5);
        } else {
            int stripeY = Math.min(c2y-4,dy+7+dir*2);
            p.drawLine(dx+3,stripeY,cxp-3,stripeY+hy/2);
        }

        p.setColor(pal.hi);
        int markX = (dir==1||dir==2)?bx-4:dx+3;
        int markY = (dir>=2)?Math.min(d2y-5,dy+hh-5):dy+5;
        p.fillRectangle(markX,markY,2,2);
    }

    private void drawBallSprite(Pixmap p, PropType type, Palette pal, int dir) {
        int r=Math.max(8,Math.round(type.w*13f));
        int cx=FRAME/2,cy=FRAME/2+4;
        p.setColor(pal.outline); p.fillCircle(cx,cy,r+1);
        p.setColor(pal.dark); p.fillCircle(cx,cy,r);
        p.setColor(pal.base); p.fillCircle(cx-1,cy-1,r-2);
        p.setColor(pal.mid); p.fillCircle(cx-2,cy-2,Math.max(3,r-5));
        int[][] h={{-5,-5},{5,-5},{5,5},{-5,5}};
        p.setColor(pal.light); p.fillRectangle(cx+h[dir][0]-2,cy+h[dir][1]-2,5,4);
        p.setColor(pal.hi); p.fillRectangle(cx+h[dir][0]-1,cy+h[dir][1]-1,2,2);
        p.setColor(pal.outline);
        if(dir%2==0) p.drawLine(cx-r+3,cy+3,cx+r-3,cy-3);
        else p.drawLine(cx-3,cy-r+3,cx+3,cy+r-3);
    }

    private void drawWheelSprite(Pixmap p, Palette pal, int dir) {
        int cx=FRAME/2,cy=FRAME/2+4,r=17;
        p.setColor(pal.outline); p.fillCircle(cx,cy,r+1);
        p.setColor(pal.dark); p.fillCircle(cx,cy,r);
        p.setColor(pal.base); p.fillCircle(cx,cy,r-4);
        p.setColor(pal.dark); p.fillCircle(cx,cy,5);
        p.setColor(pal.hi); p.fillCircle(cx,cy,2);
        p.setColor(pal.light);
        float phase=dir*MathUtils.PI/8f;
        for(int i=0;i<8;i++){
            float a=phase+i*MathUtils.PI/4f;
            int x1=cx+Math.round(MathUtils.cos(a)*5),y1=cy+Math.round(MathUtils.sin(a)*5);
            int x2=cx+Math.round(MathUtils.cos(a)*(r-5)),y2=cy+Math.round(MathUtils.sin(a)*(r-5));
            p.drawLine(x1,y1,x2,y2);
        }
        p.setColor(pal.hi); p.fillRectangle(cx-7+(dir%2)*10,cy-r+4,3,2);
    }

    private void drawRampSprite(Pixmap p, Palette pal, int dir) {
        int cx=FRAME/2, baseY=57;
        boolean flip=dir==2||dir==3;
        int left=cx-28,right=cx+28,topY=25;
        p.setColor(pal.outline);
        if(!flip) p.fillTriangle(left-1,baseY+1,right+1,baseY+1,right+1,topY-1);
        else p.fillTriangle(left-1,baseY+1,right+1,baseY+1,left-1,topY-1);
        p.setColor(pal.base);
        if(!flip) p.fillTriangle(left,baseY,right,baseY,right,topY);
        else p.fillTriangle(left,baseY,right,baseY,left,topY);
        p.setColor(pal.light);
        if(!flip) p.drawLine(left+4,baseY-3,right-3,topY+4);
        else p.drawLine(left+3,topY+4,right-4,baseY-3);
        p.setColor(pal.dark); p.drawLine(left+4,baseY-5,right-4,baseY-5);
        p.setColor(pal.hi);
        int mx=flip?left+7:right-9;
        p.fillRectangle(mx,baseY-10,3,3);
    }

    private void drawSpringSprite(Pixmap p, Palette pal, int dir) {
        int cx=FRAME/2;
        p.setColor(pal.outline); p.fillRectangle(cx-9,11,18,5); p.fillRectangle(cx-11,64,22,5);
        p.setColor(pal.base); p.fillRectangle(cx-7,12,14,3); p.fillRectangle(cx-9,65,18,3);
        int shift=(dir%2==0)?0:2;
        for(int y=20;y<60;y+=5){
            p.setColor((y/5)%2==0?pal.light:pal.mid);
            if(((y/5)+dir)%2==0) p.drawLine(cx-9+shift,y,cx+9-shift,y+3);
            else p.drawLine(cx+9-shift,y,cx-9+shift,y+3);
            p.setColor(pal.outline); p.drawPixel(cx-9+shift,y); p.drawPixel(cx+9-shift,y+3);
        }
        p.setColor(pal.hi); p.fillRectangle(cx-2+(dir-1),14,3,2);
    }

    private void drawBarrelSprite(Pixmap p, Palette pal, int dir) {
        int cx=FRAME/2,cy=39,w=23,h=36;
        p.setColor(pal.outline);
        p.fillRectangle(cx-w/2-1,cy-h/2+3,w+2,h-6);
        p.fillCircle(cx,cy-h/2+4,w/2+1); p.fillCircle(cx,cy+h/2-4,w/2+1);
        p.setColor(pal.base); p.fillRectangle(cx-w/2,cy-h/2+4,w,h-8);
        p.setColor(pal.mid); p.fillCircle(cx,cy-h/2+4,w/2-1); p.fillCircle(cx,cy+h/2-4,w/2-1);
        p.setColor(pal.dark); p.fillRectangle(cx-w/2,cy-8,w,3); p.fillRectangle(cx-w/2,cy+6,w,3);
        p.setColor(pal.light);
        int lx=dir<2?cx-w/2+4:cx+w/2-6;
        p.fillRectangle(lx,cy-h/2+7,3,h-14);
        p.setColor(pal.hi); p.fillRectangle(lx,cy-h/2+8,2,3);
    }

    private void fillQuad(Pixmap p,int x1,int y1,int x2,int y2,int x3,int y3,int x4,int y4) {
        p.fillTriangle(x1,y1,x2,y2,x3,y3);
        p.fillTriangle(x1,y1,x3,y3,x4,y4);
    }

    private void lineLoop(Pixmap p,int[] a) {
        for(int i=0;i<a.length;i+=2){
            int j=(i+2)%a.length;
            p.drawLine(a[i],a[i+1],a[j],a[j+1]);
        }
    }

    private void buildStaticPhysics() {
        platforms.add(new Platform(2.1f,-0.5f,1.5f,1.5f,0.72f));
        platforms.add(new Platform(-2.8f,1.0f,1.7f,1.4f,0.80f));
        ramps.add(new RampSurface(-3.3f,0.1f,-1.5f,-0.25f,1.35f,0f));
    }

    private void createStarterSet() {
        spawnInternal(PropType.WOOD_BALL,MaterialKind.MAHOGANY,-2.4f,-0.85f,1.9f,0f,false,nextId++);
        spawnInternal(PropType.CRATE,MaterialKind.MAHOGANY,-3.0f,1.0f,1.55f,0f,false,nextId++);
        spawnInternal(PropType.WHEEL,MaterialKind.MAHOGANY,-1.7f,-2.4f,0.75f,0.1f,false,nextId++);
        spawnInternal(PropType.CUBE,MaterialKind.METAL,-0.6f,-2.5f,0.55f,0f,false,nextId++);
        spawnInternal(PropType.CUBE,MaterialKind.MAHOGANY,2.1f,-0.5f,1.72f,0f,false,nextId++);
        spawnInternal(PropType.CUBE,MaterialKind.RUBBER,2.1f,-0.5f,2.72f,0f,false,nextId++);
        spawnInternal(PropType.CUBE,MaterialKind.MAHOGANY,2.1f,-0.5f,3.72f,0f,false,nextId++);
        spawnInternal(PropType.BEAM_SHORT,MaterialKind.MAHOGANY,1.0f,-2.3f,0.55f,0.15f,false,nextId++);
        spawnInternal(PropType.WEIGHT,MaterialKind.METAL,0.7f,2.6f,2.2f,0f,true,nextId++);
        saveWorld();
    }

    private Prop spawnInternal(PropType type, MaterialKind material, float x,float y,float z,float yaw,boolean frozen,int id) {
        if(props.size>=MAX_PROPS)return null;
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
        props.removeValue(p,true);byId.remove(p.id);
    }

    private void deleteWithUndo(Prop p) {
        if(p==null)return;
        SaveState s=p.snapshot();
        pushUndo(()->spawnFromState(s));
        removeProp(p);feedback();saveWorld();
    }

    private void duplicateWithUndo(Prop source) {
        if(source==null)return;
        Prop p=spawnInternal(source.type,source.material,source.x+0.35f,source.y-0.35f,source.z+0.45f,
                source.yaw,source.frozen,nextId++);
        if(p==null)return;
        final int id=p.id;
        pushUndo(()->{Prop q=byId.get(id);if(q!=null)removeProp(q);});
        feedback();saveWorld();
    }

    private void spawnWithUndo(PropType type) {
        float n=(props.size%5)-2;
        Prop p=spawnInternal(type,type.defaultMaterial,n*0.35f,0.1f,3.3f,0f,false,nextId++);
        if(p==null)return;
        final int id=p.id;
        pushUndo(()->{Prop q=byId.get(id);if(q!=null)removeProp(q);});
        feedback();saveWorld();
    }

    private void setFrozen(Prop p,boolean frozen,boolean record) {
        if(p==null||p.frozen==frozen)return;
        final int id=p.id;final boolean prior=p.frozen;
        if(record)pushUndo(()->{Prop q=byId.get(id);if(q!=null)setFrozen(q,prior,false);});
        p.frozen=frozen;p.vx=p.vy=p.vz=p.spin=0f;
        feedback();saveWorld();
    }

    private void cycleMaterial(Prop p) {
        if(p==null)return;
        final int id=p.id;final MaterialKind prior=p.material;
        MaterialKind next=prior==MaterialKind.MAHOGANY?MaterialKind.METAL:
                (prior==MaterialKind.METAL?MaterialKind.RUBBER:MaterialKind.MAHOGANY);
        pushUndo(()->{Prop q=byId.get(id);if(q!=null)q.material=prior;});
        p.material=next;feedback();saveWorld();
    }

    private void pushUndo(UndoAction a) {
        undo.add(a);while(undo.size>32)undo.removeIndex(0);
    }

    private void doUndo() {
        if(undo.size==0)return;
        undo.pop().undo();feedback();saveWorld();
    }

    private void physicsStep(float dt) {
        for(Prop p:props) {
            if(p.frozen)continue;

            if(p!=grabbed)p.vz-=GRAVITY*dt;
            p.x+=p.vx*dt;p.y+=p.vy*dt;p.z+=p.vz*dt;
            p.yaw+=p.spin*dt;
            p.spin*=0.992f;

            float r=p.radius();
            if(p.x-r<-ROOM){p.x=-ROOM+r;p.vx=Math.abs(p.vx)*restitution(p.material);}
            if(p.x+r> ROOM){p.x= ROOM-r;p.vx=-Math.abs(p.vx)*restitution(p.material);}
            if(p.y-r<-ROOM){p.y=-ROOM+r;p.vy=Math.abs(p.vy)*restitution(p.material);}
            if(p.y+r> ROOM){p.y= ROOM-r;p.vy=-Math.abs(p.vy)*restitution(p.material);}

            float support=supportHeight(p);
            if(p.bottom()<support) {
                p.z=support+p.type.h*0.5f;
                if(p.vz<0)p.vz=-p.vz*restitution(p.material);
                if(Math.abs(p.vz)<0.35f)p.vz=0f;
                float fr=friction(p.material);
                p.vx*=Math.max(0f,1f-fr*dt*2.4f);
                p.vy*=Math.max(0f,1f-fr*dt*2.4f);

                RampSurface ramp=rampUnder(p);
                if(ramp!=null&&Math.abs(p.vz)<0.5f)p.vx+=ramp.downhillSign()*3.0f*dt;
            }

            if(p.z<-3f){
                p.x=0;p.y=0;p.z=4f;p.vx=p.vy=p.vz=0;
            }
        }

        solveHorizontalCollisions();
    }

    private float restitution(MaterialKind m) {
        return m==MaterialKind.RUBBER?0.72f:(m==MaterialKind.METAL?0.12f:0.16f);
    }

    private float friction(MaterialKind m) {
        return m==MaterialKind.RUBBER?0.85f:(m==MaterialKind.METAL?0.45f:0.68f);
    }

    private float supportHeight(Prop p) {
        float best=0f;
        float margin=p.radius()*0.35f;

        for(Platform platform:platforms)
            if(platform.contains(p.x,p.y,margin))best=Math.max(best,platform.top);

        for(RampSurface ramp:ramps)
            if(ramp.contains(p.x,p.y,margin))best=Math.max(best,ramp.heightAt(p.x));

        for(Prop q:props) {
            if(q==p)continue;
            float dx=p.x-q.x,dy=p.y-q.y;
            float rr=(p.radius()+q.radius())*0.72f;
            if(dx*dx+dy*dy>rr*rr)continue;
            float top=q.top();
            if(top<=p.z+0.18f&&top>best)best=top;
        }
        return best;
    }

    private RampSurface rampUnder(Prop p) {
        for(RampSurface r:ramps)if(r.contains(p.x,p.y,p.radius()*0.2f))return r;
        return null;
    }

    private void solveHorizontalCollisions() {
        for(int i=0;i<props.size;i++){
            Prop a=props.get(i);
            for(int j=i+1;j<props.size;j++){
                Prop b=props.get(j);
                if(a.frozen&&b.frozen)continue;
                if(a.top()<b.bottom()+0.04f||b.top()<a.bottom()+0.04f)continue;

                float dx=b.x-a.x,dy=b.y-a.y;
                float min=(a.radius()+b.radius())*0.72f;
                float d2=dx*dx+dy*dy;
                if(d2>=min*min)continue;
                float d=(float)Math.sqrt(Math.max(d2,0.0001f));
                float nx=dx/d,ny=dy/d,penetration=min-d;

                float wa=a.frozen?0f:1f,wb=b.frozen?0f:1f,sum=wa+wb;
                if(sum<=0)continue;
                if(!a.frozen){a.x-=nx*penetration*(wa/sum);a.y-=ny*penetration*(wa/sum);}
                if(!b.frozen){b.x+=nx*penetration*(wb/sum);b.y+=ny*penetration*(wb/sum);}

                float rvx=b.vx-a.vx,rvy=b.vy-a.vy;
                float rel=rvx*nx+rvy*ny;
                if(rel<0){
                    float e=Math.min(restitution(a.material),restitution(b.material));
                    float invA=a.frozen?0f:1f/a.type.mass,invB=b.frozen?0f:1f/b.type.mass;
                    float impulse=-(1f+e)*rel/Math.max(0.0001f,invA+invB);
                    if(!a.frozen){a.vx-=impulse*nx*invA;a.vy-=impulse*ny*invA;}
                    if(!b.frozen){b.vx+=impulse*nx*invB;b.vy+=impulse*ny*invB;}
                }
            }
        }
    }

    private void updateGrab(float dt) {
        if(mode!=Mode.GRAB||grabbed==null||grabbed.frozen)return;
        float mass=Math.max(0.25f,grabbed.type.mass);
        float ms=(float)Math.pow(mass,0.34);
        float kp=55f*ms,kd=9f*(float)Math.sqrt(ms);
        float fx=(targetX-grabbed.x)*kp-grabbed.vx*kd;
        float fy=(targetY-grabbed.y)*kp-grabbed.vy*kd;
        float fz=(targetZ-grabbed.z)*kp-grabbed.vz*kd;
        float cap=95f*(float)Math.pow(Math.max(1f,mass),0.58);
        float len=(float)Math.sqrt(fx*fx+fy*fy+fz*fz);
        if(len>cap){float s=cap/len;fx*=s;fy*=s;fz*=s;}
        grabbed.vx+=fx/mass*dt;grabbed.vy+=fy/mass*dt;grabbed.vz+=fz/mass*dt;
    }

    private Vector2 project(float x,float y,float z,Vector2 out) {
        out.x=ORIGIN_X+(x-y)*ISO_X;
        out.y=ORIGIN_Y+(x+y)*ISO_Y+z*Z_PX;
        return out;
    }

    private Vector2 unprojectAtZ(float sx,float sy,float z,Vector2 out) {
        float a=(sx-ORIGIN_X)/ISO_X;
        float b=(sy-ORIGIN_Y-z*Z_PX)/ISO_Y;
        out.x=(a+b)*0.5f;
        out.y=(b-a)*0.5f;
        return out;
    }

    private int directionIndex(float yaw) {
        int q=Math.round(yaw*MathUtils.radiansToDegrees/90f)%4;
        if(q<0)q+=4;
        return q;
    }

    private void recreateBuffer() {
        if(buffer!=null)buffer.dispose();
        buffer=new FrameBuffer(Pixmap.Format.RGBA8888,W,H,false);
        buffer.getColorBufferTexture().setFilter(Texture.TextureFilter.Nearest,Texture.TextureFilter.Nearest);
        bufferRegion=new TextureRegion(buffer.getColorBufferTexture());
        bufferRegion.flip(false,true);
        updatePresentation();
    }

    private void updatePresentation() {
        float raw=Math.min(Gdx.graphics.getWidth()/(float)W,Gdx.graphics.getHeight()/(float)H);
        float integer=(float)Math.floor(raw);
        presentScale=integer>=1f?integer:raw;
        presentW=W*presentScale;presentH=H*presentScale;
        presentX=(Gdx.graphics.getWidth()-presentW)*0.5f;
        presentY=(Gdx.graphics.getHeight()-presentH)*0.5f;
    }

    private boolean insidePresentation(int sx,int syTop) {
        float by=Gdx.graphics.getHeight()-syTop;
        return sx>=presentX&&sx<=presentX+presentW&&by>=presentY&&by<=presentY+presentH;
    }

    private int virtualX(int sx) {
        return MathUtils.clamp(Math.round((sx-presentX)/presentScale),0,W-1);
    }

    private int virtualYBottom(int syTop) {
        float by=Gdx.graphics.getHeight()-syTop;
        return MathUtils.clamp(Math.round((by-presentY)/presentScale),0,H-1);
    }

    @Override
    public void render() {
        float delta=Math.min(Gdx.graphics.getDeltaTime(),0.05f);

        if(pressed!=null&&primaryPointer>=0&&!contextTriggered&&(mode==Mode.GRAB||mode==Mode.IDLE)){
            float moved=Vector2.dst(downX,downY,pointerX[primaryPointer],pointerY[primaryPointer]);
            if(moved<6f&&(TimeUtils.nanoTime()-downNanos)>550_000_000L){
                contextTriggered=true;
                contextProp=pressed;
                contextX=downX;contextY=downY;
                if(mode==Mode.GRAB)endGrab();
                primaryPointer=-1;
                mode=Mode.CONTEXT;
            }
        }

        accumulator=Math.min(MAX_ACCUM,accumulator+delta);
        while(accumulator>=FIXED_DT){
            updateGrab(FIXED_DT);
            physicsStep(FIXED_DT);
            accumulator-=FIXED_DT;
        }

        buffer.begin();
        Gdx.gl.glViewport(0,0,W,H);
        Gdx.gl.glClearColor(0.035f,0.043f,0.055f,1f);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT);
        drawRoom();
        drawProps();
        drawHud();
        buffer.end();

        updatePresentation();
        Gdx.gl.glViewport(0,0,Gdx.graphics.getWidth(),Gdx.graphics.getHeight());
        Gdx.gl.glClearColor(0.012f,0.015f,0.021f,1f);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT);
        batch.setProjectionMatrix(new Matrix4().setToOrtho2D(0,0,Gdx.graphics.getWidth(),Gdx.graphics.getHeight()));
        batch.begin();
        batch.draw(bufferRegion,presentX,presentY,presentW,presentH);
        batch.end();

        autosaveClock+=delta;
        if(autosaveClock>3f){autosaveClock=0;saveWorld();}
    }

    private void drawRoom() {
        shapes.setProjectionMatrix(projection);
        shapes.begin(ShapeRenderer.ShapeType.Filled);

        // Dark exterior vignette field.
        shapes.setColor(Color.valueOf("111822"));shapes.rect(0,0,W,H);

        Vector2 front=project(-ROOM,-ROOM,0,tmp2);
        Vector2 right=project(ROOM,-ROOM,0,tmp3);
        Vector2 back=project(ROOM,ROOM,0,tmp4);
        Vector2 left=project(-ROOM,ROOM,0,new Vector2());

        // Back walls.
        Vector2 rightTop=project(ROOM,-ROOM,WALL_H,new Vector2());
        Vector2 backTop=project(ROOM,ROOM,WALL_H,new Vector2());
        Vector2 leftTop=project(-ROOM,ROOM,WALL_H,new Vector2());

        fillQuadShape(right,back,backTop,rightTop,Color.valueOf("56433F"));
        fillQuadShape(left,back,backTop,leftTop,Color.valueOf("5D4640"));

        // Wall panel blocks.
        shapes.setColor(Color.valueOf("6B5650"));
        for(float z=0.45f;z<WALL_H;z+=0.72f){
            Vector2 a=project(ROOM,-ROOM,z,new Vector2());
            Vector2 b=project(ROOM,ROOM,z,new Vector2());
            drawPixelLine(a,b,1f);
            Vector2 c=project(-ROOM,ROOM,z,new Vector2());
            drawPixelLine(c,b,1f);
        }

        // Floor diamond.
        shapes.setColor(Color.valueOf("AD7A4E"));
        shapes.triangle(front.x,front.y,right.x,right.y,back.x,back.y);
        shapes.triangle(front.x,front.y,back.x,back.y,left.x,left.y);

        // Floor tiles.
        shapes.setColor(Color.valueOf("7D583C"));
        for(float v=-ROOM;v<=ROOM+0.01f;v+=0.86f){
            drawPixelLine(project(v,-ROOM,0.01f,new Vector2()),project(v,ROOM,0.01f,new Vector2()),1f);
            drawPixelLine(project(-ROOM,v,0.01f,new Vector2()),project(ROOM,v,0.01f,new Vector2()),1f);
        }

        // Heavy wood wall rails.
        drawWallRailY(ROOM,3.18f,0.26f,Color.valueOf("6A2F1D"));
        drawWallRailX(ROOM,3.18f,0.26f,Color.valueOf("6A2F1D"));
        drawWallRailY(ROOM,0.25f,0.22f,Color.valueOf("4B241A"));
        drawWallRailX(ROOM,0.25f,0.22f,Color.valueOf("4B241A"));

        // Corner post.
        drawVerticalPost(ROOM,ROOM,0,WALL_H,0.32f,Color.valueOf("4A2118"));
        drawVerticalPost(-ROOM,ROOM,0,WALL_H,0.26f,Color.valueOf("6B301D"));
        drawVerticalPost(ROOM,-ROOM,0,WALL_H,0.26f,Color.valueOf("6B301D"));

        // Window on left wall y=ROOM.
        drawWallRectY(ROOM,-3.35f,-1.95f,1.35f,2.65f,Color.valueOf("3A231C"));
        drawWallRectY(ROOM,-3.22f,-2.08f,1.48f,2.52f,Color.valueOf("F0B85D"));
        drawWallLineY(ROOM,-2.65f,1.48f,-2.65f,2.52f,3f,Color.valueOf("5A2A1B"));
        drawWallLineY(ROOM,-3.22f,2.0f,-2.08f,2.0f,3f,Color.valueOf("5A2A1B"));

        // Chalkboard.
        drawWallRectY(ROOM,-1.55f,0.45f,1.12f,2.35f,Color.valueOf("202B3A"));
        drawWallLineY(ROOM,-1.35f,1.34f,0.18f,2.05f,2f,Color.valueOf("8A8580"));
        drawWallLineY(ROOM,-1.20f,1.48f,-0.55f,1.22f,2f,Color.valueOf("8A8580"));

        // Warm shelf on left wall.
        drawShelfY(ROOM,-0.8f,0.75f,2.65f);

        // Shelf on right wall.
        drawShelfX(ROOM,-1.4f,0.6f,2.45f);

        // Hanging lamp near left wall.
        Vector2 lampTop=project(-1.95f,3.9f,3.45f,new Vector2());
        shapes.setColor(Color.valueOf("2D241F"));shapes.rect(lampTop.x-1,lampTop.y-24,2,24);
        shapes.setColor(Color.valueOf("C06C24"));shapes.rect(lampTop.x-7,lampTop.y-28,14,5);
        shapes.setColor(Color.valueOf("FFD66B"));shapes.rect(lampTop.x-5,lampTop.y-39,10,11);
        shapes.setColor(Color.valueOf("FFF0A8"));shapes.rect(lampTop.x-2,lampTop.y-37,4,7);

        // Wall spring, pixel-stepped.
        Vector2 s=project(3.92f,1.55f,2.95f,new Vector2());
        shapes.setColor(Color.valueOf("20252B"));shapes.rect(s.x-5,s.y-4,10,5);
        for(int i=0;i<8;i++){
            shapes.setColor(i%2==0?Color.valueOf("AEB8BF"):Color.valueOf("59636B"));
            float yy=s.y-9-i*6;
            if(i%2==0)shapes.rect(s.x-7,yy,14,3);
            else shapes.rect(s.x-4,yy,8,3);
        }
        shapes.setColor(Color.valueOf("20252B"));shapes.rect(s.x-1,s.y-58,2,10);

        // Static platforms and ramp.
        for(Platform p:platforms)drawPlatform(p);
        for(RampSurface r:ramps)drawStaticRamp(r);

        // Hanging rope + frozen weight reference.
        Vector2 ropeTop=project(0.7f,2.6f,WALL_H,new Vector2());
        Vector2 ropeBottom=project(0.7f,2.6f,2.72f,new Vector2());
        shapes.setColor(Color.valueOf("9B5D2B"));drawPixelLine(ropeTop,ropeBottom,3f);

        shapes.end();
    }

    private void drawPlatform(Platform p) {
        float z=p.top;
        Vector2 a=project(p.cx-p.w/2,p.cy-p.d/2,z,new Vector2());
        Vector2 b=project(p.cx+p.w/2,p.cy-p.d/2,z,new Vector2());
        Vector2 c=project(p.cx+p.w/2,p.cy+p.d/2,z,new Vector2());
        Vector2 d=project(p.cx-p.w/2,p.cy+p.d/2,z,new Vector2());
        fillQuadShape(a,b,c,d,Color.valueOf("7A3C23"));
        shapes.setColor(Color.valueOf("D17A3C"));drawPixelLine(d,c,2f);
        Vector2 base=project(p.cx,p.cy,0,new Vector2());
        Vector2 top=project(p.cx,p.cy,z,new Vector2());
        shapes.setColor(Color.valueOf("4A2419"));
        shapes.rect(base.x-4,base.y,8,Math.max(2,top.y-base.y));
    }

    private void drawStaticRamp(RampSurface r) {
        Vector2 a=project(r.x0,r.y0,r.z0,new Vector2());
        Vector2 b=project(r.x1,r.y0,r.z1,new Vector2());
        Vector2 c=project(r.x1,r.y1,r.z1,new Vector2());
        Vector2 d=project(r.x0,r.y1,r.z0,new Vector2());
        fillQuadShape(a,b,c,d,Color.valueOf("B65E2E"));
        shapes.setColor(Color.valueOf("E58A48"));drawPixelLine(d,c,2f);
        shapes.setColor(Color.valueOf("552619"));drawPixelLine(a,b,2f);
        for(int i=1;i<5;i++){
            float t=i/5f;
            float x=MathUtils.lerp(r.x0,r.x1,t);
            float z=MathUtils.lerp(r.z0,r.z1,t);
            drawPixelLine(project(x,r.y0,z+0.01f,new Vector2()),project(x,r.y1,z+0.01f,new Vector2()),1f);
        }
    }

    private void fillQuadShape(Vector2 a,Vector2 b,Vector2 c,Vector2 d,Color color) {
        shapes.setColor(color);
        shapes.triangle(a.x,a.y,b.x,b.y,c.x,c.y);
        shapes.triangle(a.x,a.y,c.x,c.y,d.x,d.y);
    }

    private void drawPixelLine(Vector2 a,Vector2 b,float width) {
        float dx=b.x-a.x,dy=b.y-a.y;
        float len=(float)Math.sqrt(dx*dx+dy*dy);
        if(len<0.01f)return;
        float angle=MathUtils.atan2(dy,dx)*MathUtils.radiansToDegrees;
        shapes.rect(a.x,a.y-width*0.5f,0f,width*0.5f,len,width,1f,1f,angle);
    }

    private void drawVerticalPost(float x,float y,float z0,float z1,float width,Color color) {
        Vector2 a=project(x,y,z0,new Vector2()),b=project(x,y,z1,new Vector2());
        shapes.setColor(color);shapes.rect(a.x-width*ISO_X/2,a.y,width*ISO_X,b.y-a.y);
        shapes.setColor(Color.valueOf("A55730"));shapes.rect(a.x-width*ISO_X/2+2,a.y,2,b.y-a.y);
    }

    private void drawWallRailY(float y,float z,float thickness,Color color) {
        drawWallLineY(y,-ROOM,z,ROOM,z,thickness*ISO_Y,color);
    }

    private void drawWallRailX(float x,float z,float thickness,Color color) {
        Vector2 a=project(x,-ROOM,z,new Vector2()),b=project(x,ROOM,z,new Vector2());
        shapes.setColor(color);drawPixelLine(a,b,Math.max(2f,thickness*ISO_Y));
    }

    private void drawWallRectY(float y,float x0,float x1,float z0,float z1,Color color) {
        fillQuadShape(project(x0,y,z0,new Vector2()),project(x1,y,z0,new Vector2()),
                project(x1,y,z1,new Vector2()),project(x0,y,z1,new Vector2()),color);
    }

    private void drawWallLineY(float y,float x0,float z0,float x1,float z1,float width,Color color) {
        shapes.setColor(color);drawPixelLine(project(x0,y,z0,new Vector2()),project(x1,y,z1,new Vector2()),width);
    }

    private void drawShelfY(float y,float x0,float x1,float z) {
        Vector2 a=project(x0,y-0.05f,z,new Vector2()),b=project(x1,y-0.05f,z,new Vector2());
        shapes.setColor(Color.valueOf("B65D2D"));drawPixelLine(a,b,6f);
        shapes.setColor(Color.valueOf("E28A49"));drawPixelLine(project(x0,y-0.05f,z+0.08f,new Vector2()),
                project(x1,y-0.05f,z+0.08f,new Vector2()),2f);
    }

    private void drawShelfX(float x,float y0,float y1,float z) {
        Vector2 a=project(x-0.05f,y0,z,new Vector2()),b=project(x-0.05f,y1,z,new Vector2());
        shapes.setColor(Color.valueOf("B65D2D"));drawPixelLine(a,b,6f);
        shapes.setColor(Color.valueOf("E28A49"));drawPixelLine(project(x-0.05f,y0,z+0.08f,new Vector2()),
                project(x-0.05f,y1,z+0.08f,new Vector2()),2f);
    }

    private void drawProps() {
        drawOrder.clear();drawOrder.addAll(props);
        drawOrder.sort(new Comparator<Prop>() {
            @Override public int compare(Prop a,Prop b) {
                float da=a.x+a.y+a.z*0.05f,db=b.x+b.y+b.z*0.05f;
                return Float.compare(db,da);
            }
        });

        batch.setProjectionMatrix(projection);
        batch.begin();
        for(Prop p:drawOrder){
            Vector2 s=project(p.x,p.y,p.z,tmp2);
            Texture tex=sprites.get(spriteKey(p.type,p.material))[directionIndex(p.yaw)];
            float scale=1f;
            batch.draw(tex,Math.round(s.x-FRAME*0.5f*scale),Math.round(s.y-FRAME*0.5f*scale),FRAME*scale,FRAME*scale);
        }
        batch.end();

        shapes.setProjectionMatrix(projection);
        shapes.begin(ShapeRenderer.ShapeType.Filled);
        for(Prop p:props){
            if(p!=grabbed&&!p.frozen)continue;
            Vector2 s=project(p.x,p.y,p.z,tmp2);
            Color c=p==grabbed?Color.valueOf("F4D35E"):Color.valueOf("70D6FF");
            shapes.setColor(Color.valueOf("09090B"));
            shapes.rect(s.x-15,s.y-15,31,2);shapes.rect(s.x-15,s.y+14,31,2);
            shapes.rect(s.x-15,s.y-15,2,31);shapes.rect(s.x+14,s.y-15,2,31);
            shapes.setColor(c);
            shapes.rect(s.x-14,s.y-14,7,1);shapes.rect(s.x+7,s.y-14,7,1);
            shapes.rect(s.x-14,s.y+13,7,1);shapes.rect(s.x+7,s.y+13,7,1);
        }
        shapes.end();
    }

    private void drawHud() {
        shapes.setProjectionMatrix(projection);
        shapes.begin(ShapeRenderer.ShapeType.Filled);
        panel(7,7,38,28,"25303A");
        panel(W-66,7,59,28,"25303A");
        panel(W-61,H-35,54,28,"25303A");

        if(mode==Mode.SPAWN){
            panel(0,0,W,116,"15171D");
            float cell=W/4f;
            for(int i=0;i<PropType.values().length;i++){
                int row=i/4,col=i%4;
                panel(col*cell+3,5+(3-row)*27,cell-6,24,(i&1)==0?"2C211D":"22252A");
            }
        }

        if(mode==Mode.SETTINGS)panel(W-158,H-149,151,142,"14181E");

        if(mode==Mode.CONTEXT){
            float cx=MathUtils.clamp(contextX,68,W-68),cy=MathUtils.clamp(contextY,48,H-48);
            panel(cx-66,cy-41,132,82,"17151A");
            shapes.setColor(Color.valueOf("5A321F"));shapes.rect(cx-1,cy-40,2,80);shapes.rect(cx-65,cy-1,130,2);
        }
        shapes.end();

        batch.setProjectionMatrix(projection);
        batch.begin();
        pixelText("+",20,28,Color.WHITE);
        pixelText("UNDO",W-60,26,Color.WHITE);
        pixelText("MENU",W-55,H-18,Color.WHITE);
        pixelText("PIXEL PHYSICS 2.5D",8,H-9,Color.valueOf("F0C06A"));
        pixelText("ISOMETRIC / NO 3D ENGINE",8,H-22,Color.valueOf("8FA3AD"));

        if(grabbed!=null){
            pixelText(grabbed.type.label+" / "+(directionIndex(grabbed.yaw)*90)+" DEG",8,H-36,Color.WHITE);
        }
        if(mode==Mode.SPAWN)drawSpawnDrawer();
        if(mode==Mode.SETTINGS)drawSettings();
        if(mode==Mode.CONTEXT)drawContext();
        batch.end();
    }

    private void panel(float x,float y,float w,float h,String fill) {
        shapes.setColor(Color.valueOf("09090C"));shapes.rect(x-2,y-2,w+4,h+4);
        shapes.setColor(Color.valueOf(fill));shapes.rect(x,y,w,h);
        shapes.setColor(Color.valueOf("6A5A50"));shapes.rect(x,y+h-1,w,1);shapes.rect(x,y,1,h);
        shapes.setColor(Color.valueOf("101218"));shapes.rect(x,y,w,1);shapes.rect(x+w-1,y,1,h);
    }

    private void drawSpawnDrawer() {
        pixelText("SPAWN / 4 AUTHORED ANGLES",8,109,Color.valueOf("F0C06A"));
        int preview=(int)(TimeUtils.millis()/650L)%4;
        float cell=W/4f;
        PropType[] values=PropType.values();
        for(int i=0;i<values.length;i++){
            int row=i/4,col=i%4;
            float x=col*cell+5,y=8+(3-row)*27;
            Texture tex=sprites.get(spriteKey(values[i],values[i].defaultMaterial))[preview];
            batch.draw(tex,x,y-2,24,24);
            pixelText(values[i].label,x+28,y+16,Color.WHITE);
            pixelText((preview*90)+"",x+28,y+5,Color.valueOf("8C9AA7"));
        }
    }

    private void drawSettings() {
        float x=W-149,y=H-19;
        pixelText("TRUE 2.5D PIXEL",x,y,Color.valueOf("F0C06A"));
        pixelText("RESET WORLD",x,y-32,Color.WHITE);
        pixelText("HAPTICS: "+(haptics?"ON":"OFF"),x,y-62,Color.WHITE);
        pixelText("XYZ CUSTOM PHYSICS",x,y-92,Color.valueOf("70D6FF"));
        pixelText("NO MESH / NO CAMERA3D",x,y-108,Color.valueOf("70D6FF"));
        pixelText("MENU = CLOSE",x,y-127,Color.valueOf("776F72"));
    }

    private void drawContext() {
        float cx=MathUtils.clamp(contextX,68,W-68),cy=MathUtils.clamp(contextY,48,H-48);
        pixelText(contextProp!=null&&contextProp.frozen?"UNFREEZE":"FREEZE",cx-58,cy+23,Color.valueOf("70D6FF"));
        pixelText("DELETE",cx+10,cy+23,Color.valueOf("FF806C"));
        pixelText("DUPLICATE",cx-58,cy-19,Color.WHITE);
        pixelText("MATERIAL",cx+10,cy-19,Color.valueOf("F0C06A"));
    }

    private void pixelText(String s,float x,float y,Color color) {
        font.setColor(Color.valueOf("08080A"));font.draw(batch,s,x+1,y-1);
        font.setColor(color);font.draw(batch,s,x,y);
    }

    private Prop pick(float sx,float sy) {
        drawOrder.clear();drawOrder.addAll(props);
        drawOrder.sort(new Comparator<Prop>() {
            @Override public int compare(Prop a,Prop b) {
                float da=a.x+a.y+a.z*0.05f,db=b.x+b.y+b.z*0.05f;
                return Float.compare(da,db);
            }
        });
        for(int i=drawOrder.size-1;i>=0;i--){
            Prop p=drawOrder.get(i);
            Vector2 s=project(p.x,p.y,p.z,tmp2);
            float hw=Math.max(14f,p.type.w*10f),hh=Math.max(14f,p.type.h*13f);
            if(Math.abs(sx-s.x)<=hw&&Math.abs(sy-s.y)<=hh)return p;
        }
        return null;
    }

    private void startGrab(Prop p,float sx,float sy) {
        if(p==null||p.frozen)return;
        grabbed=p;mode=Mode.GRAB;
        Vector2 w=unprojectAtZ(sx,sy,p.z,tmp2);
        grabOffsetX=w.x-p.x;grabOffsetY=w.y-p.y;
        targetX=p.x;targetY=p.y;targetZ=p.z;
    }

    private void updateGrabTarget(float sx,float sy) {
        if(grabbed==null)return;
        Vector2 w=unprojectAtZ(sx,sy,targetZ,tmp2);
        targetX=w.x-grabOffsetX;targetY=w.y-grabOffsetY;
        targetX=MathUtils.clamp(targetX,-ROOM+grabbed.radius(),ROOM-grabbed.radius());
        targetY=MathUtils.clamp(targetY,-ROOM+grabbed.radius(),ROOM-grabbed.radius());
    }

    private void endGrab() {
        grabbed=null;secondPointer=-1;
        if(mode==Mode.GRAB)mode=Mode.IDLE;
    }

    private void feedback() { if(haptics)Gdx.input.vibrate(12); }

    private void saveWorld() {
        Array<SaveState> states=new Array<>();
        for(Prop p:props)states.add(p.snapshot());
        prefs.putString("world",json.toJson(states,Array.class,SaveState.class));
        prefs.putBoolean("haptics",haptics);prefs.flush();
    }

    @SuppressWarnings("unchecked")
    private boolean restoreWorld() {
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

    private void spawnFromState(SaveState s) {
        spawnInternal(PropType.valueOf(s.type),MaterialKind.valueOf(s.material),s.x,s.y,s.z,s.yaw,s.frozen,s.id);
    }

    private void resetWorld() {
        endGrab();props.clear();byId.clear();undo.clear();nextId=1;
        pressed=null;contextProp=null;mode=Mode.IDLE;
        createStarterSet();feedback();saveWorld();
    }

    @Override
    public boolean touchDown(int screenX,int screenY,int pointer,int button) {
        if(pointer>=pointerX.length||!insidePresentation(screenX,screenY))return false;
        int vx=virtualX(screenX),vy=virtualYBottom(screenY);
        pointerX[pointer]=vx;pointerY[pointer]=vy;

        if(mode==Mode.SPAWN){
            if(vy<116){
                int col=MathUtils.clamp((int)(vx/(W/4f)),0,3);
                int row=MathUtils.clamp(3-(int)(vy/27f),0,3);
                int idx=row*4+col;
                if(idx>=0&&idx<PropType.values().length)spawnWithUndo(PropType.values()[idx]);
            }
            mode=Mode.IDLE;return true;
        }

        if(mode==Mode.SETTINGS){
            if(vx<W-158||vy<H-149){mode=Mode.IDLE;return true;}
            if(vy>H-61&&vy<H-25){resetWorld();mode=Mode.IDLE;return true;}
            if(vy>H-94&&vy<=H-61){haptics=!haptics;feedback();saveWorld();return true;}
            return true;
        }

        if(mode==Mode.CONTEXT){
            float cx=MathUtils.clamp(contextX,68,W-68),cy=MathUtils.clamp(contextY,48,H-48);
            boolean left=vx<cx,top=vy>cy;
            Prop target=contextProp;mode=Mode.IDLE;contextProp=null;
            if(target!=null){
                if(top&&left)setFrozen(target,!target.frozen,true);
                else if(top)deleteWithUndo(target);
                else if(left)duplicateWithUndo(target);
                else cycleMaterial(target);
            }
            return true;
        }

        if(vx<52&&vy<42){mode=Mode.SPAWN;feedback();return true;}
        if(vx>W-73&&vy<42){doUndo();return true;}
        if(vx>W-70&&vy>H-42){mode=Mode.SETTINGS;feedback();return true;}

        if(mode==Mode.GRAB&&pointer!=primaryPointer&&secondPointer<0){
            secondPointer=pointer;
            lastTwoDistance=pointerDistance(primaryPointer,secondPointer);
            lastTwoAngle=pointerAngle(primaryPointer,secondPointer);
            return true;
        }

        if(pointer==0&&mode==Mode.IDLE){
            primaryPointer=pointer;secondPointer=-1;
            downX=vx;downY=vy;downNanos=TimeUtils.nanoTime();contextTriggered=false;
            pressed=pick(vx,vy);
            if(pressed!=null&&!pressed.frozen)startGrab(pressed,vx,vy);
            return true;
        }
        return false;
    }

    @Override
    public boolean touchDragged(int screenX,int screenY,int pointer) {
        if(pointer>=pointerX.length)return false;
        int vx=virtualX(screenX),vy=virtualYBottom(screenY);
        pointerX[pointer]=vx;pointerY[pointer]=vy;

        if(mode==Mode.GRAB&&grabbed!=null){
            if(pointer==primaryPointer)updateGrabTarget(vx,vy);
            if(secondPointer>=0){
                float dist=pointerDistance(primaryPointer,secondPointer);
                float dd=dist-lastTwoDistance;
                targetZ=MathUtils.clamp(targetZ+dd*0.022f,grabbed.type.h*0.5f,5.0f);
                float a=pointerAngle(primaryPointer,secondPointer);
                float da=wrapAngle(a-lastTwoAngle);
                grabbed.yaw+=da;
                grabbed.spin+=da*2.5f;
                lastTwoDistance=dist;lastTwoAngle=a;
                updateGrabTarget(pointerX[primaryPointer],pointerY[primaryPointer]);
            }
            return true;
        }
        return false;
    }

    @Override
    public boolean touchUp(int screenX,int screenY,int pointer,int button) {
        if(pointer<pointerX.length){
            pointerX[pointer]=virtualX(screenX);
            pointerY[pointer]=virtualYBottom(screenY);
        }
        if(mode==Mode.GRAB){
            if(pointer==secondPointer){secondPointer=-1;return true;}
            if(pointer==primaryPointer){
                endGrab();pressed=null;primaryPointer=-1;return true;
            }
        }
        if(pointer==primaryPointer){primaryPointer=-1;pressed=null;}
        return true;
    }

    private float pointerDistance(int a,int b){return Vector2.dst(pointerX[a],pointerY[a],pointerX[b],pointerY[b]);}
    private float pointerAngle(int a,int b){return MathUtils.atan2(pointerY[b]-pointerY[a],pointerX[b]-pointerX[a]);}
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
    public void dispose() {
        saveWorld();
        for(Texture[] arr:sprites.values())for(Texture t:arr)t.dispose();
        if(buffer!=null)buffer.dispose();
        if(batch!=null)batch.dispose();
        if(shapes!=null)shapes.dispose();
        if(font!=null)font.dispose();
    }
}
