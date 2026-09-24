package com.pixelphysics.sandbox;

import static org.junit.Assert.*;
import org.junit.Test;

public class PhysicsMath25DTest {
    private static final float E = 1e-4f;

    @Test
    public void separatedBoxesDoNotCollide() {
        PhysicsMath25D.Shape a=new PhysicsMath25D.Shape().setBox(0,0,0,0.10f,0.05f);
        PhysicsMath25D.Shape b=new PhysicsMath25D.Shape().setBox(0.30f,0,0,0.10f,0.05f);
        PhysicsMath25D.Manifold m=new PhysicsMath25D.Manifold();
        assertFalse(PhysicsMath25D.collide(a,b,m));
    }

    @Test
    public void axisAlignedBoxesProduceCorrectNormalAndDepth() {
        PhysicsMath25D.Shape a=new PhysicsMath25D.Shape().setBox(0,0,0,0.10f,0.05f);
        PhysicsMath25D.Shape b=new PhysicsMath25D.Shape().setBox(0.17f,0,0,0.10f,0.05f);
        PhysicsMath25D.Manifold m=new PhysicsMath25D.Manifold();
        assertTrue(PhysicsMath25D.collide(a,b,m));
        assertEquals(1f,m.nx,E);
        assertEquals(0f,m.ny,E);
        assertEquals(0.03f,m.penetration,2e-3f);
    }

    @Test
    public void rotatedObbCollisionWorks() {
        PhysicsMath25D.Shape a=new PhysicsMath25D.Shape().setBox(0,0,(float)Math.toRadians(45),0.14f,0.035f);
        PhysicsMath25D.Shape b=new PhysicsMath25D.Shape().setBox(0.11f,0.02f,(float)Math.toRadians(-20),0.09f,0.045f);
        PhysicsMath25D.Manifold m=new PhysicsMath25D.Manifold();
        assertTrue(PhysicsMath25D.collide(a,b,m));
        assertTrue(m.penetration>0f);
        float n=(float)Math.sqrt(m.nx*m.nx+m.ny*m.ny);
        assertEquals(1f,n,2e-3f);
    }

    @Test
    public void circleCircleDepthIsMetric() {
        PhysicsMath25D.Shape a=new PhysicsMath25D.Shape().setCircle(0,0,0.10f);
        PhysicsMath25D.Shape b=new PhysicsMath25D.Shape().setCircle(0.15f,0,0.10f);
        PhysicsMath25D.Manifold m=new PhysicsMath25D.Manifold();
        assertTrue(PhysicsMath25D.collide(a,b,m));
        assertEquals(0.05f,m.penetration,2e-3f);
        assertEquals(1f,m.nx,E);
    }

    @Test
    public void circleBoxNormalPointsFromAtoB() {
        PhysicsMath25D.Shape circle=new PhysicsMath25D.Shape().setCircle(0.13f,0,0.05f);
        PhysicsMath25D.Shape box=new PhysicsMath25D.Shape().setBox(0,0,0,0.10f,0.10f);
        PhysicsMath25D.Manifold m=new PhysicsMath25D.Manifold();
        assertTrue(PhysicsMath25D.collide(circle,box,m));
        assertTrue("normal must point circle -> box",m.nx<0f);

        PhysicsMath25D.Manifold reverse=new PhysicsMath25D.Manifold();
        assertTrue(PhysicsMath25D.collide(box,circle,reverse));
        assertEquals(-m.nx,reverse.nx,2e-3f);
        assertEquals(-m.ny,reverse.ny,2e-3f);
    }

    @Test
    public void rotatedAabbExpansionIsCorrect() {
        PhysicsMath25D.Shape box=new PhysicsMath25D.Shape().setBox(0,0,(float)Math.toRadians(45),0.10f,0.10f);
        float expected=(float)(Math.sqrt(2)*0.10);
        assertEquals(expected,PhysicsMath25D.aabbHalfX(box),2e-3f);
        assertEquals(expected,PhysicsMath25D.aabbHalfY(box),2e-3f);
    }

    @Test
    public void pointInsideRespectsBoxRotation() {
        PhysicsMath25D.Shape box=new PhysicsMath25D.Shape().setBox(0,0,(float)Math.toRadians(45),0.20f,0.04f);
        assertTrue(PhysicsMath25D.pointInside(box,0.08f,0.08f,0f));
        assertFalse(PhysicsMath25D.pointInside(box,0.14f,-0.14f,0f));
    }
}
