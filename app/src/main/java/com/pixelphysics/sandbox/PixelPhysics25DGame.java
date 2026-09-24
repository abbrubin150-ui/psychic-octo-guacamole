package com.pixelphysics.sandbox;

import com.badlogic.gdx.*;
import com.badlogic.gdx.graphics.*;
import com.badlogic.gdx.graphics.g2d.*;
import com.badlogic.gdx.graphics.glutils.*;
import com.badlogic.gdx.math.*;
import com.badlogic.gdx.physics.box2d.*;
import com.badlogic.gdx.utils.*;

import java.util.Locale;

public class PixelPhysics25DGame extends ApplicationAdapter implements InputProcessor {
    private static final int W = 480;
    private static final int H = 320;
    private static final float STEP = 1f/60f;
    private static final float MAX_ACCUM = 0.12f;

    // The room is simulated on a 2D floor plane. Rendering maps that plane to
    // isometric screen coordinates. Height is an independent scalar channel.
    private static final float ROOM = 8f;
    private static final float ISO_X = 21f;
    private static final float ISO_Y = 10f;
    private static final float ORIGIN_X = 240f;
    private static final float ORIGIN_Y = 46f;
    private static final float Z_SCALE = 23f;

    private static final int FRAME = 112;
    private static final int MAX_PROPS = 80;

    private SpriteBatch batch;
    private ShapeRenderer shapes;
    private BitmapFont font;
    private FrameBuffer buffer;
    private TextureRegion bufferRegion;
    private final Matrix4 projection = new Matrix4();

    private float presentScale=1f,presentX,presentY,presentW=W,presentH=H;

    private World world;
    private final Array<Body> boundaries = new Array<>();
    private final Array<Prop> props = new Array<>();
    private final ObjectMap<Integer,Prop> byId = new ObjectMap<>();
    private final ObjectMap<String,Texture[]> spriteSets = new ObjectMap<>();
    private final Array<UndoAction> undo = new Array<>();
    private int nextId=1;

    private Preferences prefs;
    private Json json;
    private float accumulator;
    private float autosaveClock;
    private boolean haptics=true;

    private enum Mode { IDLE, GRAB, CONTEXT, SPAWN, SETTINGS }
    private Mode mode=Mode.IDLE;

    private final int[] tx=new int[10];
    private final int[] ty=new int[10];
    private int primaryPointer=-1,secondPointer=-1;
    private float downX,downY;
    private long downNanos;
    private float lastTwoAngle,lastTwoDistance;
    private float accumulatedTorque;
    private Prop pressed,grabbed,contextProp;
    private final Vector2 localGrab=new Vector2();
    private final Vector2 grabTarget=new Vector2();
    private float contextX,contextY;
    private boolean contextTriggered;

    private final Vector2 tmp2a=new Vector2();
    private final Vector2 tmp2b=new Vector2();
    private final Vector2 tmpScreen=new Vector2();

    private enum MaterialKind { MAHOGANY, METAL, RUBBER }

    private enum PropType {
        CUBE("CUBE",0.82f,0.82f,26,false,1.0f,MaterialKind.MAHOGANY),
        BEAM_SHORT("BEAM S",2.0f,0.42f,14,false,1.3f,MaterialKind.MAHOGANY),
        BEAM_LONG("BEAM L",3.0f,0.42f,14,false,2.0f,MaterialKind.MAHOGANY),
        PLANK("PLANK",2.4f,0.62f,10,false,1.2f,MaterialKind.MAHOGANY),
        WOOD_BALL("WOOD BALL",0.72f,0.72f,21,true,0.7f,MaterialKind.MAHOGANY),
        METAL_BALL("METAL BALL",0.72f,0.72f,21,true,4.0f,MaterialKind.METAL),
        WEIGHT("WEIGHT",0.9f,0.9f,31,false,8.0f,MaterialKind.METAL),
        WHEEL("WHEEL",0.92f,0.92f,29,true,1.1f,MaterialKind.MAHOGANY),
        RUBBER_BALL("RUBBER",0.78f,0.78f,22,true,0.8f,MaterialKind.RUBBER),
        BARREL("BARREL",0.86f,0.86f,36,false,3.1f,MaterialKind.METAL),
        CRATE("CRATE",1.05f,1.05f,32,false,1.7f,MaterialKind.MAHOGANY),
        RAMP("RAMP",2.1f,1.0f,22,false,1.5f,MaterialKind.MAHOGANY);

        final String label;
        final float w,d;
        final int visualH;
        final boolean circle;
        final float mass;
        final MaterialKind defaultMaterial;

        PropType(String label,float w,float d,int visualH,boolean circle,float mass,MaterialKind m){
            this.label=label;this.w=w;this.d=d;this.visualH=visualH;
            this.circle=circle;this.mass=mass;this.defaultMaterial=m;
        }
    }

    private static class Palette {
        final Color outline,dark,base,mid,light,hi;
        Palette(String a,String b,String c,String d,String e,String f){
            outline=Color.valueOf(a);dark=Color.valueOf(b);base=Color.valueOf(c);
            mid=Color.valueOf(d);light=Color.valueOf(e);hi=Color.valueOf(f);
        }
    }

    private class Prop {
        int id;
        PropType type;
        MaterialKind material;
        Body body;
        boolean frozen;
        float z;
        float vz;

        SaveState save(){
            SaveState s=new SaveState();
            s.id=id;s.type=type.name();s.material=material.name();
            s.x=body.getPosition().x;s.y=body.getPosition().y;
            s.angle=body.getAngle();s.z=z;s.vz=vz;s.frozen=frozen;
            return s;
        }
    }

    private static class SaveState {
        public int id;
        public String type,material;
        public float x,y,angle,z,vz;
        public boolean frozen;
        public SaveState(){}
    }

    private interface UndoAction { void undo(); }

    @Override
    public void create(){
        Locale.setDefault(Locale.US);
        prefs=Gdx.app.getPreferences("pixel-physics-isometric-v1");
        json=new Json();

        batch=new SpriteBatch();
        shapes=new ShapeRenderer();
        font=new BitmapFont();
        font.getData().setScale(0.72f);
        font.getRegion().getTexture().setFilter(Texture.TextureFilter.Nearest,Texture.TextureFilter.Nearest);
        projection.setToOrtho2D(0,0,W,H);

        // Top-down floor physics. The apparent third dimension is not a 3D engine.
        world=new World(new Vector2(0,0),true);
        buildRoomBounds();
        createAllSprites();
        recreateBuffer();

        if(!restoreWorld())createStarterSet();
        Gdx.input.setInputProcessor(this);
    }

    private Palette palette(MaterialKind k){
        switch(k){
            case METAL:return new Palette("12151B","2B3038","525C66","74808A","A7B0B7","E3E5E7");
            case RUBBER:return new Palette("111713","1C271F","304133","475A49","6A806B","A2B4A0");
            default:return new Palette("26100C","451710","702519","98391F","C45A2C","EE8744");
        }
    }

    private void buildRoomBounds(){
        createBound(-0.12f,ROOM/2f,0.24f,ROOM+0.4f);
        createBound(ROOM+0.12f,ROOM/2f,0.24f,ROOM+0.4f);
        createBound(ROOM/2f,-0.12f,ROOM+0.4f,0.24f);
        createBound(ROOM/2f,ROOM+0.12f,ROOM+0.4f,0.24f);
    }

    private void createBound(float x,float y,float w,float h){
        BodyDef bd=new BodyDef();bd.type=BodyDef.BodyType.StaticBody;bd.position.set(x,y);
        Body b=world.createBody(bd);
        PolygonShape ps=new PolygonShape();ps.setAsBox(w/2f,h/2f);
        FixtureDef fd=new FixtureDef();fd.shape=ps;fd.friction=0.7f;fd.restitution=0.35f;
        b.createFixture(fd);ps.dispose();boundaries.add(b);
    }

    private void createAllSprites(){
        for(PropType t:PropType.values()){
            for(MaterialKind m:MaterialKind.values()){
                Texture[] four=new Texture[4];
                for(int d=0;d<4;d++)four[d]=makeSprite(t,m,d);
                spriteSets.put(spriteKey(t,m),four);
            }
        }
    }

    private String spriteKey(PropType t,MaterialKind m){return t.name()+":"+m.name();}

    private Texture makeSprite(PropType type,MaterialKind material,int dir){
        Pixmap p=new Pixmap(FRAME,FRAME,Pixmap.Format.RGBA8888);
        p.setBlending(Pixmap.Blending.None);p.setColor(0,0,0,0);p.fill();
        p.setBlending(Pixmap.Blending.SourceOver);
        Palette pal=palette(material);

        if(type==PropType.WOOD_BALL||type==PropType.METAL_BALL||type==PropType.RUBBER_BALL){
            drawIsoBall(p,type,pal,dir);
        }else if(type==PropType.WHEEL){
            drawIsoWheel(p,pal,dir);
        }else if(type==PropType.CRATE){
            drawIsoCrate(p,pal,dir);
        }else if(type==PropType.BARREL){
            drawIsoBarrel(p,pal,dir);
        }else if(type==PropType.RAMP){
            drawIsoRamp(p,pal,dir);
        }else{
            drawIsoPrism(p,type,pal,dir);
        }

        Texture tex=new Texture(p);p.dispose();
        tex.setFilter(Texture.TextureFilter.Nearest,Texture.TextureFilter.Nearest);
        return tex;
    }

    private void drawIsoPrism(Pixmap p,PropType type,Palette pal,int dir){
        int cx=FRAME/2;
        int baseY=78;
        int h=type.visualH;
        int sx=Math.max(7,Math.round(type.w*13f));
        int sy=Math.max(4,Math.round(type.d*7f));
        if((dir&1)==1){int q=sx;sx=Math.max(7,Math.round(type.d*13f));sy=Math.max(4,Math.round(type.w*7f));}
        int topY=baseY-h;

        // left visible side
        fillQuad(p,cx-sx,topY,cx,topY+sy,cx,baseY+sy,cx-sx,baseY,pal.dark);
        // right visible side
        fillQuad(p,cx,topY+sy,cx+sx,topY,cx+sx,baseY,cx,baseY+sy,pal.base);
        // top face
        fillQuad(p,cx,topY-sy,cx+sx,topY,cx,topY+sy,cx-sx,topY,pal.mid);

        p.setColor(pal.outline);
        p.drawLine(cx,topY-sy,cx+sx,topY);p.drawLine(cx+sx,topY,cx,topY+sy);
        p.drawLine(cx,topY+sy,cx-sx,topY);p.drawLine(cx-sx,topY,cx,topY-sy);
        p.drawLine(cx-sx,topY,cx-sx,baseY);p.drawLine(cx+sx,topY,cx+sx,baseY);
        p.drawLine(cx,topY+sy,cx,baseY+sy);

        // Direction-specific authored face marks. They are not image rotation.
        p.setColor(pal.light);
        if(dir==0){p.drawLine(cx-sx+3,topY+3,cx-4,topY+sy-2);p.fillRectangle(cx+4,topY-sy+3,3,2);}
        if(dir==1){p.drawLine(cx+3,topY+sy-3,cx+sx-3,topY+2);p.fillRectangle(cx-sx+4,topY+4,3,2);}
        if(dir==2){p.drawLine(cx+4,topY+sy-3,cx+sx-3,topY+3);p.fillRectangle(cx-6,topY-sy+3,3,2);}
        if(dir==3){p.drawLine(cx-sx+3,topY+2,cx-3,topY+sy-3);p.fillRectangle(cx+sx-7,topY+4,3,2);}

        if(type==PropType.BEAM_SHORT||type==PropType.BEAM_LONG||type==PropType.PLANK){
            p.setColor(pal.hi);
            for(int i=-sx+5;i<sx-4;i+=8)p.drawPixel(cx+i,topY+(Math.abs(i)%Math.max(2,sy)));
            p.setColor(pal.dark);
            for(int yy=topY+6;yy<baseY-2;yy+=5)p.drawLine(cx-sx+2,yy,cx-2,yy+sy-2);
        }

        if(type==PropType.WEIGHT){
            p.setColor(pal.outline);p.fillRectangle(cx-4,topY+4,8,8);
            p.setColor(pal.hi);p.fillRectangle(cx-2,topY+6,3,3);
        }
    }

    private void drawIsoCrate(Pixmap p,Palette pal,int dir){
        PropType fake=PropType.CRATE;
        drawIsoPrism(p,fake,pal,dir);
        int cx=FRAME/2,top=78-fake.visualH,sx=14,sy=7;
        p.setColor(pal.outline);
        if((dir&1)==0){
            p.drawLine(cx-sx+3,top+5,cx-2,75);
            p.drawLine(cx-2,top+6,cx-sx+4,75);
            p.drawLine(cx+3,top+sy+4,cx+sx-3,75);
            p.drawLine(cx+sx-3,top+5,cx+3,75);
        }else{
            p.drawLine(cx-sx+3,top+4,cx-2,75);
            p.drawLine(cx+2,top+sy+4,cx+sx-3,75);
        }
    }

    private void drawIsoBall(Pixmap p,PropType type,Palette pal,int dir){
        int cx=FRAME/2,cy=67,r=Math.max(8,Math.round(type.w*14f));
        p.setColor(pal.outline);p.fillCircle(cx,cy,r+1);
        p.setColor(pal.dark);p.fillCircle(cx,cy,r);
        p.setColor(pal.base);p.fillCircle(cx-1,cy-2,r-2);
        p.setColor(pal.mid);p.fillCircle(cx-2,cy-4,Math.max(3,r-5));
        int[][] off={{-4,-5},{4,-5},{5,3},{-5,3}};
        p.setColor(pal.light);p.fillRectangle(cx+off[dir][0]-2,cy+off[dir][1]-1,5,3);
        p.setColor(pal.hi);p.fillRectangle(cx+off[dir][0]-1,cy+off[dir][1],2,1);
        p.setColor(pal.outline);
        if((dir&1)==0)p.drawLine(cx-r+4,cy+3,cx+r-4,cy+3);
        else p.drawLine(cx,cy-r+4,cx,cy+r-4);
    }

    private void drawIsoWheel(Pixmap p,Palette pal,int dir){
        int cx=FRAME/2,cy=67,r=15;
        p.setColor(pal.outline);p.fillCircle(cx,cy,r+1);
        p.setColor(pal.dark);p.fillCircle(cx,cy,r);
        p.setColor(pal.base);p.fillCircle(cx,cy,r-4);
        p.setColor(pal.outline);p.fillCircle(cx,cy,4);
        p.setColor(pal.hi);p.fillCircle(cx,cy,2);
        p.setColor(pal.light);
        if((dir&1)==0){
            p.fillRectangle(cx-r+5,cy-1,(r-5)*2,3);
            p.fillRectangle(cx-1,cy-r+5,3,(r-5)*2);
        }else{
            for(int i=-9;i<=9;i++){
                p.drawPixel(cx+i,cy+i/2);
                p.drawPixel(cx+i,cy-i/2);
            }
        }
    }

    private void drawIsoBarrel(Pixmap p,Palette pal,int dir){
        int cx=FRAME/2,top=44,w=(dir&1)==0?17:14,h=36;
        p.setColor(pal.outline);p.fillRectangle(cx-w-1,top+4,w*2+2,h-8);
        p.fillCircle(cx,top+4,w+1);p.fillCircle(cx,top+h-4,w+1);
        p.setColor(pal.base);p.fillRectangle(cx-w,top+4,w*2,h-8);
        p.fillCircle(cx,top+4,w);p.fillCircle(cx,top+h-4,w);
        p.setColor(pal.dark);p.fillRectangle(cx-w,top+9,w*2,3);p.fillRectangle(cx-w,top+h-12,w*2,3);
        p.setColor(pal.light);
        if(dir<2)p.fillRectangle(cx-w+4,top+6,3,h-12);
        else p.fillRectangle(cx+w-7,top+6,3,h-12);
        p.setColor(pal.hi);p.drawLine(cx-w+4,top+2,cx+w-4,top+2);
    }

    private void drawIsoRamp(Pixmap p,Palette pal,int dir){
        int cx=FRAME/2,cy=77;
        p.setColor(pal.outline);
        if((dir&1)==0){
            p.fillTriangle(cx-31,cy+7,cx+31,cy+7,cx+23,cy-20);
            p.setColor(pal.base);p.fillTriangle(cx-29,cy+5,cx+29,cy+5,cx+22,cy-17);
            p.setColor(pal.light);p.drawLine(cx-24,cy+1,cx+20,cy-15);
        }else{
            p.fillTriangle(cx-22,cy+15,cx+22,cy-15,cx+22,cy+15);
            p.setColor(pal.base);p.fillTriangle(cx-19,cy+13,cx+19,cy-13,cx+19,cy+13);
            p.setColor(pal.light);p.drawLine(cx-14,cy+8,cx+16,cy-11);
        }
        p.setColor(pal.hi);
        int hx=(dir==0||dir==3)?cx-17:cx+12;
        p.fillRectangle(hx,cy-1,3,2);
    }

    private void fillQuad(Pixmap p,int x1,int y1,int x2,int y2,int x3,int y3,int x4,int y4,Color c){
        p.setColor(c);
        p.fillTriangle(x1,y1,x2,y2,x3,y3);
        p.fillTriangle(x1,y1,x3,y3,x4,y4);
    }

    private void createStarterSet(){
        spawnInternal(PropType.CUBE,MaterialKind.MAHOGANY,1.8f,2.0f,0f,0f,false,nextId++);
        spawnInternal(PropType.CRATE,MaterialKind.MAHOGANY,5.5f,2.0f,0f,0f,false,nextId++);
        spawnInternal(PropType.WHEEL,MaterialKind.MAHOGANY,2.8f,1.1f,0f,0.15f,false,nextId++);
        spawnInternal(PropType.METAL_BALL,MaterialKind.METAL,4.0f,4.8f,0f,1.0f,false,nextId++);
        spawnInternal(PropType.BEAM_SHORT,MaterialKind.MAHOGANY,4.6f,1.3f,0.2f,0f,false,nextId++);
        spawnInternal(PropType.RUBBER_BALL,MaterialKind.RUBBER,6.2f,4.1f,0f,0.65f,false,nextId++);
        saveWorld();
    }

    private Prop spawnInternal(PropType type,MaterialKind mat,float x,float y,float angle,float z,boolean frozen,int requestedId){
        if(props.size>=MAX_PROPS)return null;
        Prop p=new Prop();
        p.id=requestedId>0?requestedId:nextId++;nextId=Math.max(nextId,p.id+1);
        p.type=type;p.material=mat;p.frozen=frozen;p.z=Math.max(0,z);

        BodyDef bd=new BodyDef();
        bd.type=frozen?BodyDef.BodyType.StaticBody:BodyDef.BodyType.DynamicBody;
        bd.position.set(MathUtils.clamp(x,0.3f,ROOM-0.3f),MathUtils.clamp(y,0.3f,ROOM-0.3f));
        bd.angle=angle;bd.linearDamping=1.2f;bd.angularDamping=1.8f;
        p.body=world.createBody(bd);p.body.setUserData(p.id);

        Shape shape;
        float area;
        if(type.circle){
            CircleShape cs=new CircleShape();float r=Math.max(type.w,type.d)/2f;
            cs.setRadius(r);shape=cs;area=MathUtils.PI*r*r;
        }else{
            PolygonShape ps=new PolygonShape();
            ps.setAsBox(type.w/2f,type.d/2f);shape=ps;area=type.w*type.d;
        }
        FixtureDef fd=new FixtureDef();fd.shape=shape;
        fd.density=Math.max(0.2f,type.mass/Math.max(0.08f,area));
        fd.friction=mat==MaterialKind.RUBBER?0.9f:(mat==MaterialKind.METAL?0.5f:0.72f);
        fd.restitution=mat==MaterialKind.RUBBER?0.62f:(mat==MaterialKind.METAL?0.15f:0.22f);
        p.body.createFixture(fd);shape.dispose();

        props.add(p);byId.put(p.id,p);return p;
    }

    private void updateHeight(float dt){
        for(Prop p:props){
            if(p.frozen)continue;
            if(p==grabbed&&mode==Mode.GRAB){p.vz=0f;continue;}
            if(p.z>0f||p.vz>0f){
                p.vz-=10.5f*dt;
                p.z+=p.vz*dt;
                if(p.z<=0f){
                    p.z=0f;
                    float bounce=p.material==MaterialKind.RUBBER?0.55f:(p.material==MaterialKind.METAL?0.24f:0.18f);
                    if(Math.abs(p.vz)>1.0f){
                        p.vz=-p.vz*bounce;
                        if(haptics&&Math.abs(p.vz)>1.2f)Gdx.input.vibrate(8);
                    }else p.vz=0f;
                }
            }
        }
    }

    private void removeProp(Prop p){
        if(p==null)return;
        if(grabbed==p)endGrab();
        props.removeValue(p,true);byId.remove(p.id);world.destroyBody(p.body);
    }

    private void pushUndo(UndoAction a){undo.add(a);while(undo.size>32)undo.removeIndex(0);}
    private void doUndo(){if(undo.size==0)return;undo.pop().undo();feedback();saveWorld();}

    private void spawnWithUndo(PropType t){
        float a=(props.size*0.73f)%5.5f;
        Prop p=spawnInternal(t,t.defaultMaterial,1.3f+a,1.4f+((props.size*0.49f)%4.8f),0f,1.1f,false,nextId++);
        if(p==null)return;final int id=p.id;
        pushUndo(()->{Prop q=byId.get(id);if(q!=null)removeProp(q);});
        feedback();saveWorld();
    }

    private void deleteWithUndo(Prop p){
        if(p==null)return;SaveState s=p.save();
        pushUndo(()->spawnFromState(s));removeProp(p);feedback();saveWorld();
    }

    private void duplicateWithUndo(Prop s){
        if(s==null)return;Vector2 q=s.body.getPosition();
        Prop p=spawnInternal(s.type,s.material,q.x+0.35f,q.y+0.35f,s.body.getAngle(),s.z+0.35f,s.frozen,nextId++);
        if(p==null)return;final int id=p.id;
        pushUndo(()->{Prop z=byId.get(id);if(z!=null)removeProp(z);});
        feedback();saveWorld();
    }

    private void setFrozen(Prop p,boolean frozen,boolean record){
        if(p==null||p.frozen==frozen)return;final int id=p.id;final boolean prior=p.frozen;
        if(record)pushUndo(()->{Prop q=byId.get(id);if(q!=null)setFrozen(q,prior,false);});
        p.frozen=frozen;p.body.setLinearVelocity(0,0);p.body.setAngularVelocity(0);
        p.body.setType(frozen?BodyDef.BodyType.StaticBody:BodyDef.BodyType.DynamicBody);
        p.vz=0;feedback();saveWorld();
    }

    private void cycleMaterial(Prop p){
        if(p==null)return;final int id=p.id;final MaterialKind prior=p.material;
        MaterialKind next=prior==MaterialKind.MAHOGANY?MaterialKind.METAL:
                (prior==MaterialKind.METAL?MaterialKind.RUBBER:MaterialKind.MAHOGANY);
        pushUndo(()->{Prop q=byId.get(id);if(q!=null)applyMaterial(q,prior);});
        applyMaterial(p,next);feedback();saveWorld();
    }

    private void applyMaterial(Prop p,MaterialKind k){
        p.material=k;
        for(Fixture f:p.body.getFixtureList()){
            f.setFriction(k==MaterialKind.RUBBER?0.9f:(k==MaterialKind.METAL?0.5f:0.72f));
            f.setRestitution(k==MaterialKind.RUBBER?0.62f:(k==MaterialKind.METAL?0.15f:0.22f));
        }
    }

    private void startGrab(Prop p,float screenX,float screenY){
        if(p==null||p.frozen)return;
        mode=Mode.GRAB;grabbed=p;
        Vector2 w=screenToWorld(screenX,screenY,p.z,tmp2a);
        localGrab.set(p.body.getLocalPoint(w));
        grabTarget.set(w);p.body.setAwake(true);
    }

    private void endGrab(){
        grabbed=null;secondPointer=-1;accumulatedTorque=0;
        if(mode==Mode.GRAB)mode=Mode.IDLE;
    }

    private void updateGrab(){
        if(mode!=Mode.GRAB||grabbed==null||grabbed.frozen)return;
        Vector2 current=grabbed.body.getWorldPoint(localGrab);
        Vector2 vel=grabbed.body.getLinearVelocityFromWorldPoint(current);
        float mass=Math.max(0.2f,grabbed.body.getMass());
        float ms=(float)Math.pow(mass,0.32);
        float kp=34f*ms,kd=7.2f*(float)Math.sqrt(ms),maxF=55f*(float)Math.pow(Math.max(1f,mass),0.52);
        tmp2a.set(grabTarget).sub(current).scl(kp).mulAdd(vel,-kd);
        if(tmp2a.len2()>maxF*maxF)tmp2a.nor().scl(maxF);
        grabbed.body.applyForce(tmp2a,current,true);
        if(Math.abs(accumulatedTorque)>0.0001f){
            grabbed.body.applyTorque(MathUtils.clamp(accumulatedTorque*2.8f,-14f,14f),true);
            accumulatedTorque*=0.55f;
        }
    }

    private Vector2 iso(float x,float y,float z,Vector2 out){
        out.x=ORIGIN_X+(x-y)*ISO_X;
        out.y=ORIGIN_Y+(x+y)*ISO_Y+z*Z_SCALE;
        return out;
    }

    private Vector2 screenToWorld(float sx,float sy,float z,Vector2 out){
        float yy=sy-z*Z_SCALE;
        float a=(sx-ORIGIN_X)/ISO_X;
        float b=(yy-ORIGIN_Y)/ISO_Y;
        out.x=(a+b)*0.5f;
        out.y=(b-a)*0.5f;
        return out;
    }

    private int directionIndex(float radians){
        int q=Math.round(radians*MathUtils.radiansToDegrees/90f)%4;
        return q<0?q+4:q;
    }

    private float depth(Prop p){Vector2 q=p.body.getPosition();return q.x+q.y;}

    private Prop screenPick(float sx,float sy){
        Prop best=null;
        float bestDepth=-999f;
        for(Prop p:props){
            Vector2 q=p.body.getPosition();iso(q.x,q.y,p.z,tmpScreen);
            float halfW=Math.max(18f,(p.type.w+p.type.d)*11f);
            float top=tmpScreen.y+p.type.visualH+16f;
            if(sx>=tmpScreen.x-halfW&&sx<=tmpScreen.x+halfW&&sy>=tmpScreen.y-10f&&sy<=top){
                float d=depth(p)-p.z*0.05f;
                // Front-most wins.
                if(best==null||d<bestDepth){best=p;bestDepth=d;}
            }
        }
        return best;
    }

    @Override
    public void render(){
        float dt=Math.min(Gdx.graphics.getDeltaTime(),0.05f);

        if(pressed!=null&&primaryPointer>=0&&!contextTriggered&&(mode==Mode.GRAB||mode==Mode.IDLE)){
            float moved=Vector2.dst(downX,downY,tx[primaryPointer],ty[primaryPointer]);
            if(moved<7f&&(TimeUtils.nanoTime()-downNanos)>550_000_000L){
                contextTriggered=true;contextProp=pressed;contextX=downX;contextY=H-downY;
                if(mode==Mode.GRAB)endGrab();
                primaryPointer=-1;mode=Mode.CONTEXT;
            }
        }

        accumulator=Math.min(MAX_ACCUM,accumulator+dt);
        while(accumulator>=STEP){
            updateGrab();
            world.step(STEP,6,2);
            updateHeight(STEP);
            accumulator-=STEP;
        }

        buffer.begin();
        Gdx.gl.glViewport(0,0,W,H);
        Gdx.gl.glClearColor(0.035f,0.039f,0.052f,1f);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT);
        drawRoom();
        drawPropShadows();
        drawProps();
        drawHud();
        buffer.end();

        updatePresentation();
        Gdx.gl.glViewport(0,0,Gdx.graphics.getWidth(),Gdx.graphics.getHeight());
        Gdx.gl.glClearColor(0.012f,0.015f,0.022f,1f);
        Gdx.gl.glClear(GL20.GL_COLOR_BUFFER_BIT);
        batch.setProjectionMatrix(new Matrix4().setToOrtho2D(0,0,Gdx.graphics.getWidth(),Gdx.graphics.getHeight()));
        batch.begin();batch.draw(bufferRegion,presentX,presentY,presentW,presentH);batch.end();

        autosaveClock+=dt;
        if(autosaveClock>3f){autosaveClock=0;saveWorld();}
    }

    private void drawRoom(){
        shapes.setProjectionMatrix(projection);
        shapes.begin(ShapeRenderer.ShapeType.Filled);

        Vector2 front=iso(0,0,0,new Vector2());
        Vector2 right=iso(ROOM,0,0,new Vector2());
        Vector2 back=iso(ROOM,ROOM,0,new Vector2());
        Vector2 left=iso(0,ROOM,0,new Vector2());

        // Floor diamond.
        shapes.setColor(Color.valueOf("9A7558"));
        shapes.triangle(front.x,front.y,right.x,right.y,back.x,back.y);
        shapes.triangle(front.x,front.y,back.x,back.y,left.x,left.y);

        // Floor tile bands.
        for(int i=0;i<=8;i++){
            Vector2 a=iso(i,0,0,tmp2a),b=iso(i,ROOM,0,tmp2b);
            shapes.setColor(Color.valueOf("6E5547"));shapes.rectLine(a.x,a.y,b.x,b.y,1f);
            a=iso(0,i,0,tmp2a);b=iso(ROOM,i,0,tmp2b);
            shapes.rectLine(a.x,a.y,b.x,b.y,1f);
        }

        float wallH=103f;
        // Left rear wall.
        quad(left.x,left.y,back.x,back.y,back.x,back.y+wallH,left.x,left.y+wallH,Color.valueOf("5A4A46"));
        // Right rear wall.
        quad(back.x,back.y,right.x,right.y,right.x,right.y+wallH,back.x,back.y+wallH,Color.valueOf("65504A"));

        // Wall panel seams.
        shapes.setColor(Color.valueOf("473B3A"));
        for(int i=1;i<8;i+=2){
            Vector2 a=iso(0,i,0,tmp2a);
            shapes.rect(a.x-1,a.y+5,2,wallH-7);
            a=iso(i,0,0,tmp2a);
            shapes.rect(a.x-1,a.y+5,2,wallH-7);
        }

        // Heavy wood base rails.
        drawIsoBeam(left.x,left.y+3,back.x,back.y+3,9,Color.valueOf("5A2A1B"),Color.valueOf("B55B31"));
        drawIsoBeam(back.x,back.y+3,right.x,right.y+3,9,Color.valueOf("5A2A1B"),Color.valueOf("B55B31"));

        // Top rails.
        drawIsoBeam(left.x,left.y+wallH-8,back.x,back.y+wallH-8,12,Color.valueOf("642B1A"),Color.valueOf("D66B35"));
        drawIsoBeam(back.x,back.y+wallH-8,right.x,right.y+wallH-8,12,Color.valueOf("642B1A"),Color.valueOf("D66B35"));

        // Corner posts.
        shapes.setColor(Color.valueOf("512416"));
        shapes.rect(back.x-6,back.y-3,12,wallH+10);
        shapes.rect(left.x-5,left.y-2,10,wallH+4);
        shapes.rect(right.x-5,right.y-2,10,wallH+4);

        // Window on left wall.
        shapes.setColor(Color.valueOf("2B1A17"));shapes.rect(101,207,55,58);
        shapes.setColor(Color.valueOf("F1B65D"));shapes.rect(106,212,45,47);
        shapes.setColor(Color.valueOf("FFF0A4"));shapes.rect(109,215,39,41);
        shapes.setColor(Color.valueOf("6B3520"));shapes.rect(127,213,4,45);shapes.rect(107,233,43,4);

        // Chalk board.
        shapes.setColor(Color.valueOf("202835"));shapes.rect(164,205,82,51);
        shapes.setColor(Color.valueOf("7C8795"));shapes.rect(168,209,74,2);
        shapes.setColor(Color.valueOf("A2A0A0"));shapes.rectLine(178,220,224,244,2);
        shapes.rectLine(181,240,229,217,2);

        // Shelves.
        shapes.setColor(Color.valueOf("4A2317"));shapes.rect(281,234,85,8);
        shapes.setColor(Color.valueOf("C06A36"));shapes.rect(285,240,77,5);
        shapes.setColor(Color.valueOf("4A2317"));shapes.rect(336,189,70,7);
        shapes.setColor(Color.valueOf("C06A36"));shapes.rect(340,195,62,5);

        // Shelf blocks.
        drawTinyIsoCube(302,248,Color.valueOf("BD3C2E"));
        drawTinyIsoCube(326,248,Color.valueOf("4477B2"));
        drawTinyIsoCube(350,248,Color.valueOf("78A046"));

        // Hanging spring and weight.
        shapes.setColor(Color.valueOf("8E99A8"));
        for(int i=0;i<7;i++){
            float yy=213-i*6;
            shapes.circle(391,yy,6,12);
            shapes.setColor(Color.valueOf("313944"));shapes.circle(391,yy,4,12);
            shapes.setColor(Color.valueOf("8E99A8"));
        }
        shapes.setColor(Color.valueOf("2B3039"));shapes.rect(387,169,8,18);
        shapes.setColor(Color.valueOf("626C78"));shapes.rect(383,163,16,7);

        // Hanging rope + weight on back wall.
        shapes.setColor(Color.valueOf("A86B2E"));shapes.rect(270,203,4,58);
        shapes.setColor(Color.valueOf("2F3239"));shapes.rect(255,183,34,20);
        shapes.setColor(Color.valueOf("747B84"));shapes.rect(259,187,26,12);

        // Front wooden frame.
        drawIsoBeam(front.x,front.y,right.x,right.y,12,Color.valueOf("562417"),Color.valueOf("C76332"));
        drawIsoBeam(front.x,front.y,left.x,left.y,12,Color.valueOf("562417"),Color.valueOf("C76332"));

        shapes.end();

        batch.setProjectionMatrix(projection);
        batch.begin();
        pixelText("PIXEL PHYSICS / ISOMETRIC 2.5D",8,H-8,Color.valueOf("F0C06A"));
        pixelText("2D ENGINE + HEIGHT CHANNEL",8,H-21,Color.valueOf("91A3AF"));
        batch.end();
    }

    private void quad(float x1,float y1,float x2,float y2,float x3,float y3,float x4,float y4,Color c){
        shapes.setColor(c);
        shapes.triangle(x1,y1,x2,y2,x3,y3);
        shapes.triangle(x1,y1,x3,y3,x4,y4);
    }

    private void drawIsoBeam(float x1,float y1,float x2,float y2,float thick,Color dark,Color light){
        shapes.setColor(dark);shapes.rectLine(x1,y1,x2,y2,thick);
        shapes.setColor(light);shapes.rectLine(x1,y1+thick*0.28f,x2,y2+thick*0.28f,Math.max(1f,thick*0.26f));
    }

    private void drawTinyIsoCube(float x,float y,Color c){
        shapes.setColor(c);shapes.rect(x-8,y-8,16,16);
        shapes.setColor(Color.valueOf("F2A15E"));shapes.rect(x-6,y+5,11,2);
        shapes.setColor(Color.valueOf("391B18"));shapes.rect(x+5,y-6,2,11);
    }

    private void drawPropShadows(){
        shapes.setProjectionMatrix(projection);
        shapes.begin(ShapeRenderer.ShapeType.Filled);
        for(Prop p:props){
            Vector2 q=p.body.getPosition();iso(q.x,q.y,0,tmpScreen);
            float s=Math.max(7f,(p.type.w+p.type.d)*7f);
            float alpha=MathUtils.clamp(0.24f-p.z*0.035f,0.05f,0.24f);
            shapes.setColor(0.05f,0.04f,0.05f,alpha);
            shapes.ellipse(tmpScreen.x-s,tmpScreen.y-4,s*2,8,18);
        }
        shapes.end();
    }

    private void drawProps(){
        props.sort((a,b)->Float.compare(depth(b),depth(a)));
        batch.setProjectionMatrix(projection);
        batch.begin();
        for(Prop p:props){
            Vector2 q=p.body.getPosition();iso(q.x,q.y,p.z,tmpScreen);
            Texture tex=spriteSets.get(spriteKey(p.type,p.material))[directionIndex(p.body.getAngle())];
            batch.draw(tex,Math.round(tmpScreen.x)-FRAME/2f,Math.round(tmpScreen.y)-82f);
        }
        batch.end();
    }

    private void drawHud(){
        shapes.setProjectionMatrix(projection);
        shapes.begin(ShapeRenderer.ShapeType.Filled);

        panel(7,7,38,28,"26303A");
        panel(W-65,7,58,28,"26303A");
        panel(W-61,H-35,54,28,"26303A");

        if(grabbed!=null){
            Vector2 q=grabbed.body.getPosition();iso(q.x,q.y,grabbed.z,tmpScreen);
            shapes.setColor(Color.valueOf("0B0B0D"));shapes.rect(tmpScreen.x-18,tmpScreen.y-8,36,2);
            shapes.rect(tmpScreen.x-18,tmpScreen.y+20,36,2);shapes.rect(tmpScreen.x-18,tmpScreen.y-8,2,30);
            shapes.rect(tmpScreen.x+16,tmpScreen.y-8,2,30);
            shapes.setColor(Color.valueOf("F4D35E"));shapes.rect(tmpScreen.x-16,tmpScreen.y-6,7,1);
            shapes.rect(tmpScreen.x+9,tmpScreen.y-6,7,1);shapes.rect(tmpScreen.x-16,tmpScreen.y+18,7,1);
            shapes.rect(tmpScreen.x+9,tmpScreen.y+18,7,1);
        }

        if(mode==Mode.SPAWN){
            panel(0,0,W,122,"151419");
            int cols=4;float cell=W/(float)cols;
            for(int i=0;i<PropType.values().length;i++){
                int row=i/cols,col=i%cols;
                panel(col*cell+3,5+(2-row)*36,cell-6,32,(i&1)==0?"2C201C":"242027");
            }
        }

        if(mode==Mode.SETTINGS)panel(W-162,H-155,154,147,"151A20");

        if(mode==Mode.CONTEXT){
            float cx=MathUtils.clamp(contextX,68,W-68),cy=MathUtils.clamp(contextY,50,H-50);
            panel(cx-66,cy-42,132,84,"17151A");
            shapes.setColor(Color.valueOf("5A321F"));shapes.rect(cx-1,cy-41,2,82);shapes.rect(cx-65,cy-1,130,2);
        }
        shapes.end();

        batch.setProjectionMatrix(projection);
        batch.begin();
        pixelText("+",20,28,Color.WHITE);
        pixelText("UNDO",W-59,26,Color.WHITE);
        pixelText("MENU",W-55,H-18,Color.WHITE);
        if(grabbed!=null){
            pixelText(grabbed.type.label+"  Z:"+String.format(Locale.US,"%.1f",grabbed.z),8,H-34,Color.WHITE);
        }
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
        pixelText("SPAWN / 4 ISOMETRIC ANGLES PER PROP",8,115,Color.valueOf("F0C06A"));
        int preview=(int)(TimeUtils.millis()/650L)%4;
        int cols=4;float cell=W/(float)cols;
        PropType[] vals=PropType.values();
        for(int i=0;i<vals.length;i++){
            int row=i/cols,col=i%cols;
            float cx=col*cell+18,cy=7+(2-row)*36+16;
            Texture t=spriteSets.get(spriteKey(vals[i],vals[i].defaultMaterial))[preview];
            batch.draw(t,cx-18,cy-18,36,36);
            pixelText(vals[i].label,col*cell+40,cy+5,Color.WHITE);
        }
        String[] dirs={"NE","SE","SW","NW"};
        pixelText("ANGLE: "+dirs[preview],W-82,115,Color.valueOf("8FA3AD"));
    }

    private void drawSettings(){
        float x=W-153,y=H-18;
        pixelText("ISOMETRIC 2.5D",x,y,Color.valueOf("F0C06A"));
        pixelText("RESET WORLD",x,y-33,Color.WHITE);
        pixelText("HAPTICS: "+(haptics?"ON":"OFF"),x,y-64,Color.WHITE);
        pixelText("CANVAS 480x320",x,y-94,Color.valueOf("70D6FF"));
        pixelText("NO 3D ENGINE",x,y-111,Color.valueOf("70D6FF"));
        pixelText("PINCH = HEIGHT",x,y-128,Color.valueOf("8FA3AD"));
    }

    private void drawContext(){
        float cx=MathUtils.clamp(contextX,68,W-68),cy=MathUtils.clamp(contextY,50,H-50);
        pixelText(contextProp!=null&&contextProp.frozen?"UNFREEZE":"FREEZE",cx-58,cy+24,Color.valueOf("70D6FF"));
        pixelText("DELETE",cx+10,cy+24,Color.valueOf("FF806C"));
        pixelText("DUPLICATE",cx-58,cy-20,Color.WHITE);
        pixelText("MATERIAL",cx+10,cy-20,Color.valueOf("F0C06A"));
    }

    private void pixelText(String s,float x,float y,Color c){
        font.setColor(Color.valueOf("08080A"));font.draw(batch,s,x+1,y-1);
        font.setColor(c);font.draw(batch,s,x,y);
    }

    private void recreateBuffer(){
        if(buffer!=null)buffer.dispose();
        buffer=new FrameBuffer(Pixmap.Format.RGBA8888,W,H,false);
        buffer.getColorBufferTexture().setFilter(Texture.TextureFilter.Nearest,Texture.TextureFilter.Nearest);
        bufferRegion=new TextureRegion(buffer.getColorBufferTexture());bufferRegion.flip(false,true);
        updatePresentation();
    }

    private void updatePresentation(){
        float raw=Math.min(Gdx.graphics.getWidth()/(float)W,Gdx.graphics.getHeight()/(float)H);
        float integer=(float)Math.floor(raw);presentScale=integer>=1?integer:raw;
        presentW=W*presentScale;presentH=H*presentScale;
        presentX=(Gdx.graphics.getWidth()-presentW)/2f;presentY=(Gdx.graphics.getHeight()-presentH)/2f;
    }

    private boolean insidePresentation(int sx,int sy){
        return sx>=presentX&&sx<=presentX+presentW&&sy>=presentY&&sy<=presentY+presentH;
    }

    private int virtualX(int sx){return MathUtils.clamp(Math.round((sx-presentX)/presentScale),0,W-1);}
    private int virtualTopY(int sy){return MathUtils.clamp(Math.round((sy-presentY)/presentScale),0,H-1);}
    private float canvasYFromTop(int topY){return H-topY;}

    private void feedback(){if(haptics)Gdx.input.vibrate(12);}

    private void saveWorld(){
        Array<SaveState> states=new Array<>();for(Prop p:props)states.add(p.save());
        prefs.putString("world",json.toJson(states,Array.class,SaveState.class));
        prefs.putBoolean("haptics",haptics);prefs.flush();
    }

    @SuppressWarnings("unchecked")
    private boolean restoreWorld(){
        haptics=prefs.getBoolean("haptics",true);
        String data=prefs.getString("world","");
        if(data.isEmpty())return false;
        try{
            Array<SaveState> states=json.fromJson(Array.class,SaveState.class,data);
            if(states==null||states.size==0)return false;
            for(SaveState s:states)spawnFromState(s);return true;
        }catch(Exception e){return false;}
    }

    private void spawnFromState(SaveState s){
        Prop p=spawnInternal(PropType.valueOf(s.type),MaterialKind.valueOf(s.material),s.x,s.y,s.angle,s.z,s.frozen,s.id);
        if(p!=null)p.vz=s.vz;
    }

    private void resetWorld(){
        endGrab();Array<Prop> copy=new Array<>(props);for(Prop p:copy)removeProp(p);
        undo.clear();nextId=1;mode=Mode.IDLE;contextProp=null;pressed=null;createStarterSet();feedback();saveWorld();
    }

    @Override
    public boolean touchDown(int screenX,int screenY,int pointer,int button){
        if(pointer>=tx.length||!insidePresentation(screenX,screenY))return false;
        int x=virtualX(screenX),topY=virtualTopY(screenY);
        float cy=canvasYFromTop(topY);
        tx[pointer]=x;ty[pointer]=topY;

        if(mode==Mode.SPAWN){
            if(topY>H-122+18){
                int cols=4;float cell=W/(float)cols;
                int row=MathUtils.clamp((topY-(H-122+18))/36,0,2);
                int col=MathUtils.clamp((int)(x/cell),0,3);
                int idx=row*cols+col;
                if(idx>=0&&idx<PropType.values().length)spawnWithUndo(PropType.values()[idx]);
            }
            mode=Mode.IDLE;return true;
        }

        if(mode==Mode.SETTINGS){
            float uy=cy;
            if(x<W-162||uy<H-155){mode=Mode.IDLE;return true;}
            if(uy>H-68&&uy<H-27){resetWorld();mode=Mode.IDLE;return true;}
            if(uy>H-100&&uy<=H-68){haptics=!haptics;feedback();saveWorld();return true;}
            return true;
        }

        if(mode==Mode.CONTEXT){
            float cx=MathUtils.clamp(contextX,68,W-68),ccy=MathUtils.clamp(contextY,50,H-50);
            boolean left=x<cx,top=cy>ccy;Prop target=contextProp;
            mode=Mode.IDLE;contextProp=null;
            if(target!=null){
                if(top&&left)setFrozen(target,!target.frozen,true);
                else if(top)deleteWithUndo(target);
                else if(left)duplicateWithUndo(target);
                else cycleMaterial(target);
            }
            return true;
        }

        if(x<52&&topY>H-42){mode=Mode.SPAWN;feedback();return true;}
        if(x>W-72&&topY>H-42){doUndo();return true;}
        if(x>W-70&&topY<42){mode=Mode.SETTINGS;feedback();return true;}

        if(mode==Mode.GRAB&&pointer!=primaryPointer&&secondPointer<0){
            secondPointer=pointer;lastTwoAngle=pointerAngle(primaryPointer,secondPointer);
            lastTwoDistance=pointerDistance(primaryPointer,secondPointer);return true;
        }

        if(pointer==0&&mode==Mode.IDLE){
            primaryPointer=pointer;secondPointer=-1;downX=x;downY=topY;downNanos=TimeUtils.nanoTime();contextTriggered=false;
            pressed=screenPick(x,cy);
            if(pressed!=null&&!pressed.frozen)startGrab(pressed,x,cy);
            return true;
        }
        return false;
    }

    @Override
    public boolean touchDragged(int screenX,int screenY,int pointer){
        if(pointer>=tx.length)return false;
        int x=virtualX(screenX),topY=virtualTopY(screenY);float cy=canvasYFromTop(topY);
        tx[pointer]=x;ty[pointer]=topY;

        if(mode==Mode.GRAB&&grabbed!=null){
            if(pointer==primaryPointer){
                screenToWorld(x,cy,grabbed.z,grabTarget);
                grabTarget.x=MathUtils.clamp(grabTarget.x,0.25f,ROOM-0.25f);
                grabTarget.y=MathUtils.clamp(grabTarget.y,0.25f,ROOM-0.25f);
            }
            if(secondPointer>=0){
                float dist=pointerDistance(primaryPointer,secondPointer);
                float dd=dist-lastTwoDistance;
                grabbed.z=MathUtils.clamp(grabbed.z+dd*0.018f,0f,3.2f);
                float a=pointerAngle(primaryPointer,secondPointer);
                accumulatedTorque+=wrapAngle(a-lastTwoAngle);
                lastTwoAngle=a;lastTwoDistance=dist;
            }
            return true;
        }
        return false;
    }

    @Override
    public boolean touchUp(int screenX,int screenY,int pointer,int button){
        if(mode==Mode.GRAB){
            if(pointer==secondPointer){secondPointer=-1;return true;}
            if(pointer==primaryPointer){endGrab();pressed=null;primaryPointer=-1;return true;}
        }
        if(pointer==primaryPointer){primaryPointer=-1;pressed=null;return true;}
        return false;
    }

    private float pointerAngle(int a,int b){return MathUtils.atan2(ty[b]-ty[a],tx[b]-tx[a]);}
    private float pointerDistance(int a,int b){return Vector2.dst(tx[a],ty[a],tx[b],ty[b]);}
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
