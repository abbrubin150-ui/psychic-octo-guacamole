package com.pixelphysics.sandbox;

import android.content.Context;
import android.content.SharedPreferences;
import android.graphics.*;
import android.os.SystemClock;
import android.view.HapticFeedbackConstants;
import android.view.MotionEvent;
import android.view.View;

import org.json.JSONArray;
import org.json.JSONObject;

import java.util.*;

public final class PixelPhysicsView extends View {
    private static final int VW=480, VH=270, SPR=80, MAX_PROPS=100;
    private static final float ROOM=4.3f, WALL_H=3.6f;
    private static final float ORIGIN_X=240f, ORIGIN_Y=190f, ISO_X=26f, ISO_Y=10f, ZPX=24f;
    private static final float DT=1f/60f, MAX_ACCUM=0.12f, G=9.81f;

    private final Bitmap frame=Bitmap.createBitmap(VW,VH,Bitmap.Config.ARGB_8888);
    private final Canvas worldCanvas=new Canvas(frame);
    private final Paint p=new Paint();
    private final Paint nearest=new Paint();
    private final RectF dst=new RectF();
    private final SharedPreferences prefs;

    private float presentScale=1f, presentX=0, presentY=0, presentW=VW, presentH=VH;
    private boolean running=true;
    private long lastNs=0;
    private float accumulator=0, autosave=0;

    private enum MaterialKind { MAHOGANY, METAL, RUBBER }
    private enum Mode { IDLE, GRAB, CONTEXT, SPAWN, SETTINGS }

    private enum PropType {
        CUBE("CUBE",1f,1f,1f,1f,MaterialKind.MAHOGANY,false),
        BEAM_SHORT("BEAM S",2.2f,.55f,.50f,1.3f,MaterialKind.MAHOGANY,false),
        BEAM_LONG("BEAM L",3.3f,.55f,.50f,2.0f,MaterialKind.MAHOGANY,false),
        PLANK("PLANK",2.7f,1f,.28f,1.2f,MaterialKind.MAHOGANY,false),
        WOOD_BALL("WOOD BALL",.9f,.9f,.9f,.7f,MaterialKind.MAHOGANY,true),
        METAL_BALL("METAL BALL",.9f,.9f,.9f,4f,MaterialKind.METAL,true),
        WEIGHT("WEIGHT",1f,1f,1.05f,10f,MaterialKind.METAL,false),
        WHEEL("WHEEL",1.25f,.42f,1.25f,1.1f,MaterialKind.MAHOGANY,false),
        RUBBER_BALL("RUBBER",1f,1f,1f,.85f,MaterialKind.RUBBER,true),
        BARREL("BARREL",1.1f,1.1f,1.5f,3.6f,MaterialKind.METAL,false),
        CRATE("CRATE",1.4f,1.4f,1.4f,1.8f,MaterialKind.MAHOGANY,false),
        RAMP("RAMP",2.4f,1.4f,.8f,1.8f,MaterialKind.MAHOGANY,false),
        SPRING("SPRING",.8f,.8f,2f,1.5f,MaterialKind.METAL,false);
        final String label; final float w,d,h,mass; final MaterialKind def; final boolean sphere;
        PropType(String label,float w,float d,float h,float mass,MaterialKind def,boolean sphere){
            this.label=label;this.w=w;this.d=d;this.h=h;this.mass=mass;this.def=def;this.sphere=sphere;
        }
    }

    private static final class Prop {
        int id; PropType type; MaterialKind material;
        float x,y,z,vx,vy,vz,yaw,spin; boolean frozen;
        float radius(){ return Math.max(type.w,type.d)*.5f; }
        float bottom(){ return z-type.h*.5f; }
        float top(){ return z+type.h*.5f; }
    }
    private static final class Platform {
        float cx,cy,w,d,top;
        Platform(float cx,float cy,float w,float d,float top){this.cx=cx;this.cy=cy;this.w=w;this.d=d;this.top=top;}
        boolean contains(float x,float y,float m){return Math.abs(x-cx)<=w*.5f+m&&Math.abs(y-cy)<=d*.5f+m;}
    }
    private static final class RampSurface {
        float x0,x1,y0,y1,z0,z1;
        RampSurface(float x0,float x1,float y0,float y1,float z0,float z1){this.x0=x0;this.x1=x1;this.y0=y0;this.y1=y1;this.z0=z0;this.z1=z1;}
        boolean contains(float x,float y,float m){return x>=Math.min(x0,x1)-m&&x<=Math.max(x0,x1)+m&&y>=Math.min(y0,y1)-m&&y<=Math.max(y0,y1)+m;}
        float heightAt(float x){float t=clamp((x-x0)/(x1-x0),0,1);return lerp(z0,z1,t);}
        float downhill(){return z1>z0?-1f:1f;}
    }

    private final ArrayList<Prop> props=new ArrayList<>();
    private final HashMap<Integer,Prop> byId=new HashMap<>();
    private final ArrayList<Platform> platforms=new ArrayList<>();
    private final ArrayList<RampSurface> ramps=new ArrayList<>();
    private final HashMap<String,Bitmap[]> sprites=new HashMap<>();
    private final ArrayDeque<Runnable> undo=new ArrayDeque<>();
    private int nextId=1;

    private Mode mode=Mode.IDLE;
    private Prop grabbed, pressed, contextProp;
    private int primaryId=-1, secondId=-1;
    private float downX,downY,contextX,contextY;
    private long downNs;
    private boolean contextTriggered=false;
    private float grabOffsetX,grabOffsetY,targetX,targetY,targetZ,lastTwoDistance,lastTwoAngle;
    private boolean haptics=true;

    public PixelPhysicsView(Context context){
        super(context);
        setFocusable(true);
        setKeepScreenOn(true);
        p.setAntiAlias(false); p.setFilterBitmap(false); p.setDither(false);
        p.setTypeface(Typeface.create(Typeface.MONOSPACE,Typeface.BOLD));
        nearest.setAntiAlias(false); nearest.setFilterBitmap(false); nearest.setDither(false);
        prefs=context.getSharedPreferences("pixel-physics-native-25d-v1",Context.MODE_PRIVATE);
        buildStatics();
        createSprites();
        if(!restoreWorld()) createStarterSet();
        lastNs=System.nanoTime();
    }

    public void resumeLoop(){running=true;lastNs=System.nanoTime();postInvalidateOnAnimation();}
    public void pauseLoop(){saveWorld();running=false;}

    @Override protected void onSizeChanged(int w,int h,int oldw,int oldh){updatePresentation(w,h);}

    private void updatePresentation(int w,int h){
        float raw=Math.min(w/(float)VW,h/(float)VH);
        float integer=(float)Math.floor(raw);
        presentScale=integer>=1f?integer:raw;
        presentW=VW*presentScale;presentH=VH*presentScale;
        presentX=(w-presentW)*.5f;presentY=(h-presentH)*.5f;
        dst.set(presentX,presentY,presentX+presentW,presentY+presentH);
    }

    @Override protected void onDraw(Canvas canvas){
        super.onDraw(canvas);
        long now=System.nanoTime();
        float delta=Math.min(.05f,(now-lastNs)/1_000_000_000f);
        lastNs=now;
        if(running){
            if(checkLongPress(now)){}
            accumulator=Math.min(MAX_ACCUM,accumulator+delta);
            while(accumulator>=DT){updateGrab(DT);physicsStep(DT);accumulator-=DT;}
            autosave+=delta;if(autosave>3f){autosave=0;saveWorld();}
        }
        drawFrame();
        canvas.drawColor(Color.rgb(12,16,22));
        canvas.drawBitmap(frame,null,dst,nearest);
        if(running)postInvalidateOnAnimation();
    }

    private boolean checkLongPress(long now){
        if(pressed==null||primaryId<0||contextTriggered)return false;
        if(mode!=Mode.GRAB&&mode!=Mode.IDLE)return false;
        if(now-downNs<550_000_000L)return false;
        contextTriggered=true;contextProp=pressed;contextX=downX;contextY=downY;
        grabbed=null;mode=Mode.CONTEXT;primaryId=-1;
        haptic();return true;
    }

    private void drawFrame(){
        worldCanvas.drawColor(Color.rgb(17,24,34));
        drawRoom(worldCanvas);
        drawProps(worldCanvas);
        drawHud(worldCanvas);
    }

    private void drawRoom(Canvas c){
        PointF front=project(-ROOM,-ROOM,0), right=project(ROOM,-ROOM,0), back=project(ROOM,ROOM,0), left=project(-ROOM,ROOM,0);
        PointF rt=project(ROOM,-ROOM,WALL_H), bt=project(ROOM,ROOM,WALL_H), lt=project(-ROOM,ROOM,WALL_H);

        quad(c,right,back,bt,rt,0xFF55433E);
        quad(c,left,back,bt,lt,0xFF5E4842);
        quad(c,front,right,back,left,0xFFB07C4E);

        stroke(0xFF7A563B,1);
        for(float v=-ROOM;v<=ROOM+.01f;v+=.86f){
            line(c,project(v,-ROOM,.01f),project(v,ROOM,.01f));
            line(c,project(-ROOM,v,.01f),project(ROOM,v,.01f));
        }

        stroke(0xFF725B54,1);
        for(float z=.45f;z<WALL_H;z+=.72f){
            line(c,project(ROOM,-ROOM,z),project(ROOM,ROOM,z));
            line(c,project(-ROOM,ROOM,z),project(ROOM,ROOM,z));
        }

        wallRailY(c,ROOM,3.18f,7,0xFF6B301D);
        wallRailX(c,ROOM,3.18f,7,0xFF6B301D);
        wallRailY(c,ROOM,.25f,5,0xFF4A241A);
        wallRailX(c,ROOM,.25f,5,0xFF4A241A);
        verticalPost(c,ROOM,ROOM,WALL_H,8,0xFF4A2118);
        verticalPost(c,-ROOM,ROOM,WALL_H,7,0xFF6B301D);
        verticalPost(c,ROOM,-ROOM,WALL_H,7,0xFF6B301D);

        wallRectY(c,ROOM,-3.35f,-1.95f,1.35f,2.65f,0xFF3A231C);
        wallRectY(c,ROOM,-3.22f,-2.08f,1.48f,2.52f,0xFFF0B85D);
        stroke(0xFF5A2A1B,3);
        line(c,project(-2.65f,ROOM,1.48f),project(-2.65f,ROOM,2.52f));
        line(c,project(-3.22f,ROOM,2.0f),project(-2.08f,ROOM,2.0f));

        wallRectY(c,ROOM,-1.55f,.45f,1.12f,2.35f,0xFF202B3A);
        stroke(0xFF8A8580,2);
        line(c,project(-1.35f,ROOM,1.34f),project(.18f,ROOM,2.05f));
        line(c,project(-1.20f,ROOM,1.48f),project(-.55f,ROOM,1.22f));

        shelfY(c,ROOM,-.8f,.75f,2.65f);
        shelfX(c,ROOM,-1.4f,.6f,2.45f);

        PointF lamp=project(-1.95f,3.9f,3.45f);
        fillRect(c,lamp.x-1,lamp.y,2,24,0xFF2D241F);
        fillRect(c,lamp.x-7,lamp.y+24,14,5,0xFFC06C24);
        fillRect(c,lamp.x-5,lamp.y+29,10,11,0xFFFFD66B);
        fillRect(c,lamp.x-2,lamp.y+31,4,7,0xFFFFF0A8);

        PointF s=project(3.92f,1.55f,2.95f);
        fillRect(c,s.x-5,s.y-4,10,5,0xFF20252B);
        for(int i=0;i<8;i++){
            int col=i%2==0?0xFFAEB8BF:0xFF59636B; float yy=s.y+9+i*6;
            fillRect(c,s.x-(i%2==0?7:4),yy,i%2==0?14:8,3,col);
        }
        fillRect(c,s.x-1,s.y+58,2,10,0xFF20252B);

        for(Platform pl:platforms)drawPlatform(c,pl);
        for(RampSurface r:ramps)drawStaticRamp(c,r);

        stroke(0xFF9B5D2B,3);
        line(c,project(.7f,2.6f,WALL_H),project(.7f,2.6f,2.72f));
    }

    private void drawPlatform(Canvas c,Platform pl){
        float z=pl.top;
        PointF a=project(pl.cx-pl.w/2,pl.cy-pl.d/2,z), b=project(pl.cx+pl.w/2,pl.cy-pl.d/2,z),
                cc=project(pl.cx+pl.w/2,pl.cy+pl.d/2,z), d=project(pl.cx-pl.w/2,pl.cy+pl.d/2,z);
        quad(c,a,b,cc,d,0xFF7A3C23);
        stroke(0xFFD17A3C,2);line(c,d,cc);
        PointF base=project(pl.cx,pl.cy,0), top=project(pl.cx,pl.cy,z);
        fillRect(c,base.x-4,top.y,8,base.y-top.y,0xFF4A2419);
    }

    private void drawStaticRamp(Canvas c,RampSurface r){
        PointF a=project(r.x0,r.y0,r.z0),b=project(r.x1,r.y0,r.z1),cc=project(r.x1,r.y1,r.z1),d=project(r.x0,r.y1,r.z0);
        quad(c,a,b,cc,d,0xFFB65E2E);
        stroke(0xFFE58A48,2);line(c,d,cc);
        stroke(0xFF552619,2);line(c,a,b);
        for(int i=1;i<5;i++){float t=i/5f,x=lerp(r.x0,r.x1,t),z=lerp(r.z0,r.z1,t);line(c,project(x,r.y0,z+.01f),project(x,r.y1,z+.01f));}
    }

    private void drawProps(Canvas c){
        ArrayList<Prop> order=new ArrayList<>(props);
        order.sort((a,b)->Float.compare(b.x+b.y+b.z*.05f,a.x+a.y+a.z*.05f));
        for(Prop pr:order){
            PointF s=project(pr.x,pr.y,pr.z);
            Bitmap bm=sprites.get(spriteKey(pr.type,pr.material))[dir(pr.yaw)];
            c.drawBitmap(bm,Math.round(s.x-SPR/2f),Math.round(s.y-SPR/2f),nearest);
        }
        for(Prop pr:props){
            if(pr!=grabbed&&!pr.frozen)continue;
            PointF s=project(pr.x,pr.y,pr.z);
            int col=pr==grabbed?0xFFF4D35E:0xFF70D6FF;
            stroke(0xFF09090B,2);rectStroke(c,s.x-15,s.y-15,30,30);
            stroke(col,1);
            line(c,new PointF(s.x-14,s.y-14),new PointF(s.x-7,s.y-14));
            line(c,new PointF(s.x+7,s.y-14),new PointF(s.x+14,s.y-14));
            line(c,new PointF(s.x-14,s.y+14),new PointF(s.x-7,s.y+14));
            line(c,new PointF(s.x+7,s.y+14),new PointF(s.x+14,s.y+14));
        }
    }

    private void drawHud(Canvas c){
        panel(c,7,VH-35,38,28,0xFF25303A);
        panel(c,VW-66,VH-35,59,28,0xFF25303A);
        panel(c,VW-61,7,54,28,0xFF25303A);
        text(c,"+",20,VH-15,0xFFFFFFFF,12);
        text(c,"UNDO",VW-60,VH-16,0xFFFFFFFF,9);
        text(c,"MENU",VW-55,25,0xFFFFFFFF,9);
        text(c,"PIXEL PHYSICS 2.5D",8,13,0xFFF0C06A,9);
        text(c,"NATIVE ANDROID / NO ENGINE",8,25,0xFF8FA3AD,7);

        if(grabbed!=null)text(c,grabbed.type.label+" / "+(dir(grabbed.yaw)*90)+" DEG",8,39,0xFFFFFFFF,8);

        if(mode==Mode.SPAWN){
            panel(c,0,VH-116,VW,116,0xFF15171D);
            float cell=VW/4f;
            for(int i=0;i<PropType.values().length;i++){
                int row=i/4,col=i%4;
                panel(c,col*cell+3,VH-112+row*27,cell-6,24,(i&1)==0?0xFF2C211D:0xFF22252A);
            }
            text(c,"SPAWN / 4 AUTHORED ANGLES",8,VH-104,0xFFF0C06A,8);
            int preview=(int)((SystemClock.uptimeMillis()/650L)%4);
            PropType[] vals=PropType.values();
            for(int i=0;i<vals.length;i++){
                int row=i/4,col=i%4;float x=col*cell+5,y=VH-110+row*27;
                Bitmap bm=sprites.get(spriteKey(vals[i],vals[i].def))[preview];
                c.drawBitmap(bm,null,new RectF(x,y,x+24,y+24),nearest);
                text(c,vals[i].label,x+28,y+10,0xFFFFFFFF,7);
                text(c,String.valueOf(preview*90),x+28,y+20,0xFF8C9AA7,6);
            }
        }

        if(mode==Mode.SETTINGS){
            panel(c,VW-158,7,151,142,0xFF14181E);
            float x=VW-149;
            text(c,"TRUE 2.5D PIXEL",x,23,0xFFF0C06A,8);
            text(c,"RESET WORLD",x,54,0xFFFFFFFF,8);
            text(c,"HAPTICS: "+(haptics?"ON":"OFF"),x,84,0xFFFFFFFF,8);
            text(c,"JAVA CANVAS",x,114,0xFF70D6FF,8);
            text(c,"NO OPENGL / NO JNI",x,130,0xFF70D6FF,7);
        }

        if(mode==Mode.CONTEXT){
            float cx=clamp(contextX,68,VW-68),cy=clamp(contextY,48,VH-48);
            panel(c,cx-66,cy-41,132,82,0xFF17151A);
            fillRect(c,cx-1,cy-40,2,80,0xFF5A321F);fillRect(c,cx-65,cy-1,130,2,0xFF5A321F);
            text(c,contextProp!=null&&contextProp.frozen?"UNFREEZE":"FREEZE",cx-58,cy-17,0xFF70D6FF,7);
            text(c,"DELETE",cx+10,cy-17,0xFFFF806C,7);
            text(c,"DUPLICATE",cx-58,cy+24,0xFFFFFFFF,7);
            text(c,"MATERIAL",cx+10,cy+24,0xFFF0C06A,7);
        }
    }

    private void buildStatics(){
        platforms.add(new Platform(2.1f,-.5f,1.5f,1.5f,.72f));
        platforms.add(new Platform(-2.8f,1f,1.7f,1.4f,.80f));
        ramps.add(new RampSurface(-3.3f,.1f,-1.5f,-.25f,1.35f,0f));
    }

    private void createStarterSet(){
        spawn(PropType.WOOD_BALL,MaterialKind.MAHOGANY,-2.4f,-.85f,1.9f,0,false,nextId++);
        spawn(PropType.CRATE,MaterialKind.MAHOGANY,-3f,1f,1.55f,0,false,nextId++);
        spawn(PropType.WHEEL,MaterialKind.MAHOGANY,-1.7f,-2.4f,.75f,.1f,false,nextId++);
        spawn(PropType.CUBE,MaterialKind.METAL,-.6f,-2.5f,.55f,0,false,nextId++);
        spawn(PropType.CUBE,MaterialKind.MAHOGANY,2.1f,-.5f,1.72f,0,false,nextId++);
        spawn(PropType.CUBE,MaterialKind.RUBBER,2.1f,-.5f,2.72f,0,false,nextId++);
        spawn(PropType.CUBE,MaterialKind.MAHOGANY,2.1f,-.5f,3.72f,0,false,nextId++);
        spawn(PropType.BEAM_SHORT,MaterialKind.MAHOGANY,1f,-2.3f,.55f,.15f,false,nextId++);
        spawn(PropType.WEIGHT,MaterialKind.METAL,.7f,2.6f,2.2f,0,true,nextId++);
        saveWorld();
    }

    private Prop spawn(PropType type,MaterialKind m,float x,float y,float z,float yaw,boolean frozen,int id){
        if(props.size()>=MAX_PROPS)return null;
        Prop pr=new Prop();pr.id=id>0?id:nextId++;nextId=Math.max(nextId,pr.id+1);
        pr.type=type;pr.material=m;pr.x=x;pr.y=y;pr.z=z;pr.yaw=yaw;pr.frozen=frozen;
        props.add(pr);byId.put(pr.id,pr);return pr;
    }

    private void remove(Prop pr){if(pr==null)return;if(grabbed==pr)grabbed=null;props.remove(pr);byId.remove(pr.id);}
    private void pushUndo(Runnable r){undo.addLast(r);while(undo.size()>32)undo.removeFirst();}
    private void doUndo(){if(undo.isEmpty())return;undo.removeLast().run();haptic();saveWorld();}

    private void spawnWithUndo(PropType type){
        float n=(props.size()%5)-2;
        Prop pr=spawn(type,type.def,n*.35f,.1f,3.3f,0,false,nextId++);
        if(pr==null)return;int id=pr.id;pushUndo(()->remove(byId.get(id)));haptic();saveWorld();
    }

    private void deleteWithUndo(Prop pr){
        if(pr==null)return;JSONObject snap=toJson(pr);pushUndo(()->fromJson(snap));remove(pr);haptic();saveWorld();
    }

    private void duplicateWithUndo(Prop src){
        Prop pr=spawn(src.type,src.material,src.x+.35f,src.y-.35f,src.z+.45f,src.yaw,src.frozen,nextId++);
        if(pr==null)return;int id=pr.id;pushUndo(()->remove(byId.get(id)));haptic();saveWorld();
    }

    private void setFrozen(Prop pr,boolean frozen,boolean record){
        if(pr==null||pr.frozen==frozen)return;boolean prev=pr.frozen;int id=pr.id;
        if(record)pushUndo(()->{Prop q=byId.get(id);if(q!=null)setFrozen(q,prev,false);});
        pr.frozen=frozen;pr.vx=pr.vy=pr.vz=pr.spin=0;haptic();saveWorld();
    }

    private void cycleMaterial(Prop pr){
        MaterialKind prev=pr.material;int id=pr.id;
        MaterialKind next=prev==MaterialKind.MAHOGANY?MaterialKind.METAL:(prev==MaterialKind.METAL?MaterialKind.RUBBER:MaterialKind.MAHOGANY);
        pushUndo(()->{Prop q=byId.get(id);if(q!=null)q.material=prev;});
        pr.material=next;haptic();saveWorld();
    }

    private void physicsStep(float dt){
        for(Prop pr:props){
            if(pr.frozen)continue;
            if(pr!=grabbed)pr.vz-=G*dt;
            pr.x+=pr.vx*dt;pr.y+=pr.vy*dt;pr.z+=pr.vz*dt;pr.yaw+=pr.spin*dt;pr.spin*=.992f;
            float r=pr.radius(),e=restitution(pr.material);
            if(pr.x-r<-ROOM){pr.x=-ROOM+r;pr.vx=Math.abs(pr.vx)*e;}
            if(pr.x+r> ROOM){pr.x= ROOM-r;pr.vx=-Math.abs(pr.vx)*e;}
            if(pr.y-r<-ROOM){pr.y=-ROOM+r;pr.vy=Math.abs(pr.vy)*e;}
            if(pr.y+r> ROOM){pr.y= ROOM-r;pr.vy=-Math.abs(pr.vy)*e;}

            float support=supportHeight(pr);
            if(pr.bottom()<support){
                pr.z=support+pr.type.h*.5f;if(pr.vz<0)pr.vz=-pr.vz*e;if(Math.abs(pr.vz)<.35f)pr.vz=0;
                float fr=friction(pr.material);pr.vx*=Math.max(0,1-fr*dt*2.4f);pr.vy*=Math.max(0,1-fr*dt*2.4f);
                RampSurface rr=rampUnder(pr);if(rr!=null&&Math.abs(pr.vz)<.5f)pr.vx+=rr.downhill()*3f*dt;
            }
            if(pr.z<-3){pr.x=pr.y=0;pr.z=4;pr.vx=pr.vy=pr.vz=0;}
        }
        solveHorizontal();
    }

    private void updateGrab(float dt){
        if(mode!=Mode.GRAB||grabbed==null||grabbed.frozen)return;
        float mass=Math.max(.25f,grabbed.type.mass),ms=(float)Math.pow(mass,.34),kp=55f*ms,kd=9f*(float)Math.sqrt(ms);
        float fx=(targetX-grabbed.x)*kp-grabbed.vx*kd,fy=(targetY-grabbed.y)*kp-grabbed.vy*kd,fz=(targetZ-grabbed.z)*kp-grabbed.vz*kd;
        float cap=95f*(float)Math.pow(Math.max(1,mass),.58),len=(float)Math.sqrt(fx*fx+fy*fy+fz*fz);
        if(len>cap){float s=cap/len;fx*=s;fy*=s;fz*=s;}
        grabbed.vx+=fx/mass*dt;grabbed.vy+=fy/mass*dt;grabbed.vz+=fz/mass*dt;
    }

    private float supportHeight(Prop pr){
        float best=0,margin=pr.radius()*.35f;
        for(Platform pl:platforms)if(pl.contains(pr.x,pr.y,margin))best=Math.max(best,pl.top);
        for(RampSurface r:ramps)if(r.contains(pr.x,pr.y,margin))best=Math.max(best,r.heightAt(pr.x));
        for(Prop q:props){
            if(q==pr)continue;float dx=pr.x-q.x,dy=pr.y-q.y,rr=(pr.radius()+q.radius())*.72f;
            if(dx*dx+dy*dy>rr*rr)continue;float top=q.top();if(top<=pr.z+.18f&&top>best)best=top;
        }
        return best;
    }

    private RampSurface rampUnder(Prop pr){for(RampSurface r:ramps)if(r.contains(pr.x,pr.y,pr.radius()*.2f))return r;return null;}

    private void solveHorizontal(){
        for(int i=0;i<props.size();i++)for(int j=i+1;j<props.size();j++){
            Prop a=props.get(i),b=props.get(j);if(a.frozen&&b.frozen)continue;
            if(a.top()<b.bottom()+.04f||b.top()<a.bottom()+.04f)continue;
            float dx=b.x-a.x,dy=b.y-a.y,min=(a.radius()+b.radius())*.72f,d2=dx*dx+dy*dy;
            if(d2>=min*min)continue;float d=(float)Math.sqrt(Math.max(d2,.0001)),nx=dx/d,ny=dy/d,pen=min-d;
            float wa=a.frozen?0:1,wb=b.frozen?0:1,sum=wa+wb;if(sum<=0)continue;
            if(!a.frozen){a.x-=nx*pen*wa/sum;a.y-=ny*pen*wa/sum;}
            if(!b.frozen){b.x+=nx*pen*wb/sum;b.y+=ny*pen*wb/sum;}
            float rel=(b.vx-a.vx)*nx+(b.vy-a.vy)*ny;
            if(rel<0){float e=Math.min(restitution(a.material),restitution(b.material)),ia=a.frozen?0:1/a.type.mass,ib=b.frozen?0:1/b.type.mass;
                float imp=-(1+e)*rel/Math.max(.0001f,ia+ib);
                if(!a.frozen){a.vx-=imp*nx*ia;a.vy-=imp*ny*ia;}if(!b.frozen){b.vx+=imp*nx*ib;b.vy+=imp*ny*ib;}
            }
        }
    }

    private float restitution(MaterialKind m){return m==MaterialKind.RUBBER?.72f:(m==MaterialKind.METAL?.12f:.16f);}
    private float friction(MaterialKind m){return m==MaterialKind.RUBBER?.85f:(m==MaterialKind.METAL?.45f:.68f);}

    @Override public boolean onTouchEvent(MotionEvent e){
        int action=e.getActionMasked(),idx=e.getActionIndex();
        if(action==MotionEvent.ACTION_DOWN){
            if(!inside(e.getX(),e.getY()))return true;
            float vx=toVX(e.getX()),vy=toVY(e.getY());
            if(handleUiDown(vx,vy))return true;
            primaryId=e.getPointerId(0);downX=vx;downY=vy;downNs=System.nanoTime();contextTriggered=false;
            pressed=pick(vx,vy);if(pressed!=null&&!pressed.frozen)startGrab(pressed,vx,vy);return true;
        }
        if(action==MotionEvent.ACTION_POINTER_DOWN&&mode==Mode.GRAB&&secondId<0){
            secondId=e.getPointerId(idx);int a=e.findPointerIndex(primaryId),b=idx;if(a>=0){lastTwoDistance=distance(e,a,b);lastTwoAngle=angle(e,a,b);}return true;
        }
        if(action==MotionEvent.ACTION_MOVE&&mode==Mode.GRAB&&grabbed!=null){
            int a=e.findPointerIndex(primaryId);if(a>=0)updateTarget(toVX(e.getX(a)),toVY(e.getY(a)));
            if(secondId>=0){int b=e.findPointerIndex(secondId);if(a>=0&&b>=0){float d=distance(e,a,b),dd=(d-lastTwoDistance)/presentScale;targetZ=clamp(targetZ+dd*.022f,grabbed.type.h*.5f,5f);
                float an=angle(e,a,b),da=wrap(an-lastTwoAngle);grabbed.yaw+=da;grabbed.spin+=da*2.5f;lastTwoDistance=d;lastTwoAngle=an;}}
            return true;
        }
        if(action==MotionEvent.ACTION_POINTER_UP){
            int id=e.getPointerId(idx);if(id==secondId)secondId=-1;return true;
        }
        if(action==MotionEvent.ACTION_UP||action==MotionEvent.ACTION_CANCEL){
            int id=e.getPointerId(idx);if(id==primaryId){grabbed=null;pressed=null;primaryId=-1;secondId=-1;if(mode==Mode.GRAB)mode=Mode.IDLE;}return true;
        }
        return true;
    }

    private boolean handleUiDown(float x,float y){
        if(mode==Mode.SPAWN){
            if(y>VH-116){int col=clampi((int)(x/(VW/4f)),0,3),row=clampi((int)((y-(VH-116))/27f),0,3),n=row*4+col;if(n<PropType.values().length)spawnWithUndo(PropType.values()[n]);}
            mode=Mode.IDLE;return true;
        }
        if(mode==Mode.SETTINGS){
            if(x<VW-158||y>149){mode=Mode.IDLE;return true;}
            if(y>25&&y<64){resetWorld();mode=Mode.IDLE;return true;}
            if(y>=64&&y<98){haptics=!haptics;haptic();saveWorld();return true;}
            return true;
        }
        if(mode==Mode.CONTEXT){
            float cx=clamp(contextX,68,VW-68),cy=clamp(contextY,48,VH-48);boolean left=x<cx,top=y<cy;Prop t=contextProp;mode=Mode.IDLE;contextProp=null;
            if(t!=null){if(top&&left)setFrozen(t,!t.frozen,true);else if(top)deleteWithUndo(t);else if(left)duplicateWithUndo(t);else cycleMaterial(t);}return true;
        }
        if(x<52&&y>VH-42){mode=Mode.SPAWN;haptic();return true;}
        if(x>VW-73&&y>VH-42){doUndo();return true;}
        if(x>VW-70&&y<42){mode=Mode.SETTINGS;haptic();return true;}
        return false;
    }

    private Prop pick(float sx,float sy){
        ArrayList<Prop> order=new ArrayList<>(props);order.sort((a,b)->Float.compare(a.x+a.y+a.z*.05f,b.x+b.y+b.z*.05f));
        for(int i=order.size()-1;i>=0;i--){Prop pr=order.get(i);PointF s=project(pr.x,pr.y,pr.z);float hw=Math.max(14,pr.type.w*10),hh=Math.max(14,pr.type.h*13);
            if(Math.abs(sx-s.x)<=hw&&Math.abs(sy-s.y)<=hh)return pr;}return null;
    }

    private void startGrab(Prop pr,float sx,float sy){
        grabbed=pr;mode=Mode.GRAB;PointF w=unprojectAtZ(sx,sy,pr.z);grabOffsetX=w.x-pr.x;grabOffsetY=w.y-pr.y;targetX=pr.x;targetY=pr.y;targetZ=pr.z;
    }

    private void updateTarget(float sx,float sy){
        if(grabbed==null)return;PointF w=unprojectAtZ(sx,sy,targetZ);targetX=clamp(w.x-grabOffsetX,-ROOM+grabbed.radius(),ROOM-grabbed.radius());targetY=clamp(w.y-grabOffsetY,-ROOM+grabbed.radius(),ROOM-grabbed.radius());
    }

    private float distance(MotionEvent e,int a,int b){float dx=e.getX(b)-e.getX(a),dy=e.getY(b)-e.getY(a);return (float)Math.sqrt(dx*dx+dy*dy);}
    private float angle(MotionEvent e,int a,int b){return (float)Math.atan2(e.getY(b)-e.getY(a),e.getX(b)-e.getX(a));}
    private float wrap(float a){while(a>(float)Math.PI)a-=2f*(float)Math.PI;while(a<-(float)Math.PI)a+=2f*(float)Math.PI;return a;}

    private boolean inside(float sx,float sy){return sx>=presentX&&sx<=presentX+presentW&&sy>=presentY&&sy<=presentY+presentH;}
    private float toVX(float sx){return clamp((sx-presentX)/presentScale,0,VW-1);}
    private float toVY(float sy){return clamp((sy-presentY)/presentScale,0,VH-1);}

    private PointF project(float x,float y,float z){return new PointF(ORIGIN_X+(x-y)*ISO_X,ORIGIN_Y-(x+y)*ISO_Y-z*ZPX);}
    private PointF unprojectAtZ(float sx,float sy,float z){float a=(sx-ORIGIN_X)/ISO_X,b=(ORIGIN_Y-sy-z*ZPX)/ISO_Y;return new PointF((a+b)*.5f,(b-a)*.5f);}
    private int dir(float yaw){int q=Math.round((float)Math.toDegrees(yaw)/90f)%4;return q<0?q+4:q;}

    private void createSprites(){
        for(PropType t:PropType.values())for(MaterialKind m:MaterialKind.values()){Bitmap[] a=new Bitmap[4];for(int d=0;d<4;d++)a[d]=makeSprite(t,m,d);sprites.put(spriteKey(t,m),a);}
    }
    private String spriteKey(PropType t,MaterialKind m){return t.name()+":"+m.name();}

    private Bitmap makeSprite(PropType type,MaterialKind m,int dir){
        Bitmap b=Bitmap.createBitmap(SPR,SPR,Bitmap.Config.ARGB_8888);Canvas c=new Canvas(b);int[] pal=palette(m);
        if(type.sphere)drawBallSprite(c,type,pal,dir);else if(type==PropType.WHEEL)drawWheelSprite(c,pal,dir);else if(type==PropType.RAMP)drawRampSprite(c,pal,dir);
        else if(type==PropType.SPRING)drawSpringSprite(c,pal,dir);else if(type==PropType.BARREL)drawBarrelSprite(c,pal,dir);else drawCuboidSprite(c,type,pal,dir);
        return b;
    }

    private int[] palette(MaterialKind m){
        if(m==MaterialKind.METAL)return new int[]{0xFF101319,0xFF272D34,0xFF4A5660,0xFF73818B,0xFFAAB4BA,0xFFE0E4E6};
        if(m==MaterialKind.RUBBER)return new int[]{0xFF0F1511,0xFF1A241C,0xFF2B3B2E,0xFF435746,0xFF627A63,0xFFA5B39D};
        return new int[]{0xFF24100D,0xFF421813,0xFF6F271A,0xFF963B21,0xFFC25D30,0xFFE88A49};
    }

    private void drawCuboidSprite(Canvas c,PropType t,int[] pal,int d){
        float rw=d%2==0?t.w:t.d,rd=d%2==0?t.d:t.w;int hx=Math.min(28,Math.max(5,Math.round(rw*8.5f))),hy=Math.min(14,Math.max(3,Math.round(rd*4f))),hh=Math.min(34,Math.max(5,Math.round(t.h*15f)));
        int cx=SPR/2,cy=18+hy;PointF a=new PointF(cx,cy-hy),bb=new PointF(cx+hx,cy),cc=new PointF(cx,cy+hy),dd=new PointF(cx-hx,cy);
        PointF b2=new PointF(bb.x,bb.y+hh),c2=new PointF(cc.x,cc.y+hh),d2=new PointF(dd.x,dd.y+hh);
        quad(c,dd,cc,c2,d2,pal[2]);quad(c,bb,cc,c2,b2,d==1?pal[4]:pal[1]);quad(c,a,bb,cc,dd,d==0?pal[4]:(d==3?pal[5]:pal[3]));
        stroke(pal[0],1);line(c,a,bb);line(c,bb,cc);line(c,cc,dd);line(c,dd,a);line(c,dd,d2);line(c,cc,c2);line(c,bb,b2);line(c,d2,c2);line(c,c2,b2);
        if(t==PropType.CRATE){stroke(pal[3],1);line(c,new PointF(dd.x+4,dd.y+5),new PointF(c2.x-4,c2.y-5));line(c,new PointF(d2.x+4,d2.y-5),new PointF(cc.x-4,cc.y+5));}
        fillRect(c,(d==1||d==2)?bb.x-4:dd.x+3,(d>=2)?Math.min(d2.y-5,dd.y+hh-5):dd.y+5,2,2,pal[5]);
    }

    private void drawBallSprite(Canvas c,PropType t,int[] pal,int d){
        int r=Math.max(8,Math.round(t.w*13)),cx=SPR/2,cy=SPR/2+4;circle(c,cx,cy,r+1,pal[0]);circle(c,cx,cy,r,pal[1]);circle(c,cx-1,cy-1,r-2,pal[2]);circle(c,cx-2,cy-2,Math.max(3,r-5),pal[3]);
        int[][] h={{-5,-5},{5,-5},{5,5},{-5,5}};fillRect(c,cx+h[d][0]-2,cy+h[d][1]-2,5,4,pal[4]);fillRect(c,cx+h[d][0]-1,cy+h[d][1]-1,2,2,pal[5]);
        stroke(pal[0],1);if(d%2==0)line(c,new PointF(cx-r+3,cy+3),new PointF(cx+r-3,cy-3));else line(c,new PointF(cx-3,cy-r+3),new PointF(cx+3,cy+r-3));
    }

    private void drawWheelSprite(Canvas c,int[] pal,int d){
        int cx=SPR/2,cy=SPR/2+4,r=17;circle(c,cx,cy,r+1,pal[0]);circle(c,cx,cy,r,pal[1]);circle(c,cx,cy,r-4,pal[2]);circle(c,cx,cy,5,pal[1]);circle(c,cx,cy,2,pal[5]);
        stroke(pal[4],1);float phase=d*(float)Math.PI/8f;for(int i=0;i<8;i++){float a=phase+i*(float)Math.PI/4f;line(c,new PointF(cx+(float)Math.cos(a)*5,cy+(float)Math.sin(a)*5),new PointF(cx+(float)Math.cos(a)*(r-5),cy+(float)Math.sin(a)*(r-5)));}}
    private void drawRampSprite(Canvas c,int[] pal,int d){
        int cx=SPR/2,by=57,left=cx-28,right=cx+28,top=25;Path path=new Path();boolean flip=d==2||d==3;path.moveTo(left,by);path.lineTo(right,by);path.lineTo(flip?left:right,top);path.close();fillPath(c,path,pal[2]);stroke(pal[0],2);c.drawPath(path,p);
        stroke(pal[4],1);line(c,new PointF(flip?left+3:left+4,flip?top+4:by-3),new PointF(flip?right-4:right-3,flip?by-3:top+4));
    }
    private void drawSpringSprite(Canvas c,int[] pal,int d){
        int cx=SPR/2;fillRect(c,cx-9,11,18,5,pal[0]);fillRect(c,cx-7,12,14,3,pal[2]);fillRect(c,cx-11,64,22,5,pal[0]);int shift=d%2==0?0:2;
        for(int y=20;y<60;y+=5){stroke((y/5)%2==0?pal[4]:pal[3],2);if(((y/5)+d)%2==0)line(c,new PointF(cx-9+shift,y),new PointF(cx+9-shift,y+3));else line(c,new PointF(cx+9-shift,y),new PointF(cx-9+shift,y+3));}
    }
    private void drawBarrelSprite(Canvas c,int[] pal,int d){
        int cx=SPR/2,cy=39,w=23,h=36;fillRect(c,cx-w/2,cy-h/2+4,w,h-8,pal[2]);circle(c,cx,cy-h/2+4,w/2-1,pal[3]);circle(c,cx,cy+h/2-4,w/2-1,pal[3]);fillRect(c,cx-w/2,cy-8,w,3,pal[1]);fillRect(c,cx-w/2,cy+6,w,3,pal[1]);
        fillRect(c,d<2?cx-w/2+4:cx+w/2-6,cy-h/2+7,3,h-14,pal[4]);
    }

    private void saveWorld(){
        try{JSONArray arr=new JSONArray();for(Prop pr:props)arr.put(toJson(pr));prefs.edit().putString("world",arr.toString()).putBoolean("haptics",haptics).apply();}catch(Exception ignored){}
    }
    private boolean restoreWorld(){
        haptics=prefs.getBoolean("haptics",true);String s=prefs.getString("world","");if(s==null||s.isEmpty())return false;
        try{JSONArray a=new JSONArray(s);if(a.length()==0)return false;for(int i=0;i<a.length();i++)fromJson(a.getJSONObject(i));return true;}catch(Exception e){return false;}
    }
    private JSONObject toJson(Prop pr){
        JSONObject o=new JSONObject();try{o.put("id",pr.id);o.put("type",pr.type.name());o.put("material",pr.material.name());o.put("x",pr.x);o.put("y",pr.y);o.put("z",pr.z);o.put("yaw",pr.yaw);o.put("frozen",pr.frozen);}catch(Exception ignored){}return o;
    }
    private Prop fromJson(JSONObject o){
        try{return spawn(PropType.valueOf(o.getString("type")),MaterialKind.valueOf(o.getString("material")),(float)o.getDouble("x"),(float)o.getDouble("y"),(float)o.getDouble("z"),(float)o.getDouble("yaw"),o.getBoolean("frozen"),o.getInt("id"));}catch(Exception e){return null;}
    }
    private void resetWorld(){grabbed=pressed=contextProp=null;mode=Mode.IDLE;props.clear();byId.clear();undo.clear();nextId=1;createStarterSet();haptic();saveWorld();}

    private void haptic(){if(haptics)performHapticFeedback(HapticFeedbackConstants.CLOCK_TICK);}

    private void panel(Canvas c,float x,float y,float w,float h,int color){fillRect(c,x-2,y-2,w+4,h+4,0xFF09090C);fillRect(c,x,y,w,h,color);fillRect(c,x,y,w,1,0xFF6A5A50);fillRect(c,x,y,1,h,0xFF6A5A50);fillRect(c,x,y+h-1,w,1,0xFF101218);fillRect(c,x+w-1,y,1,h,0xFF101218);}
    private void text(Canvas c,String s,float x,float y,int color,float size){p.setStyle(Paint.Style.FILL);p.setColor(color);p.setTextSize(size);p.setStrokeWidth(1);c.drawText(s,x,y,p);}
    private void fillRect(Canvas c,float x,float y,float w,float h,int color){p.setStyle(Paint.Style.FILL);p.setColor(color);c.drawRect(x,y,x+w,y+h,p);}
    private void rectStroke(Canvas c,float x,float y,float w,float h){p.setStyle(Paint.Style.STROKE);c.drawRect(x,y,x+w,y+h,p);p.setStyle(Paint.Style.FILL);}
    private void circle(Canvas c,float x,float y,float r,int color){p.setStyle(Paint.Style.FILL);p.setColor(color);c.drawCircle(x,y,r,p);}
    private void stroke(int color,float width){p.setStyle(Paint.Style.STROKE);p.setStrokeWidth(width);p.setColor(color);}
    private void line(Canvas c,PointF a,PointF b){c.drawLine(Math.round(a.x),Math.round(a.y),Math.round(b.x),Math.round(b.y),p);}
    private void fillPath(Canvas c,Path path,int color){p.setStyle(Paint.Style.FILL);p.setColor(color);c.drawPath(path,p);}
    private void quad(Canvas c,PointF a,PointF b,PointF cc,PointF d,int color){Path path=new Path();path.moveTo(a.x,a.y);path.lineTo(b.x,b.y);path.lineTo(cc.x,cc.y);path.lineTo(d.x,d.y);path.close();fillPath(c,path,color);}
    private void wallRectY(Canvas c,float y,float x0,float x1,float z0,float z1,int color){quad(c,project(x0,y,z0),project(x1,y,z0),project(x1,y,z1),project(x0,y,z1),color);}
    private void wallRailY(Canvas c,float y,float z,float width,int color){stroke(color,width);line(c,project(-ROOM,y,z),project(ROOM,y,z));}
    private void wallRailX(Canvas c,float x,float z,float width,int color){stroke(color,width);line(c,project(x,-ROOM,z),project(x,ROOM,z));}
    private void verticalPost(Canvas c,float x,float y,float z1,float width,int color){stroke(color,width);line(c,project(x,y,0),project(x,y,z1));}
    private void shelfY(Canvas c,float y,float x0,float x1,float z){stroke(0xFFB65D2D,6);line(c,project(x0,y-.05f,z),project(x1,y-.05f,z));stroke(0xFFE28A49,2);line(c,project(x0,y-.05f,z+.08f),project(x1,y-.05f,z+.08f));}
    private void shelfX(Canvas c,float x,float y0,float y1,float z){stroke(0xFFB65D2D,6);line(c,project(x-.05f,y0,z),project(x-.05f,y1,z));stroke(0xFFE28A49,2);line(c,project(x-.05f,y0,z+.08f),project(x-.05f,y1,z+.08f));}

    private static float clamp(float v,float lo,float hi){return Math.max(lo,Math.min(hi,v));}
    private static int clampi(int v,int lo,int hi){return Math.max(lo,Math.min(hi,v));}
    private static float lerp(float a,float b,float t){return a+(b-a)*t;}
}
