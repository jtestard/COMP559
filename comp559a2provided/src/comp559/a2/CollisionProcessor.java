package comp559.a2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Matrix3d;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import no.uib.cipr.matrix.DenseVector;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.CollapsiblePanel;
import mintools.swing.VerticalFlowPanel;

/**
 * Class for detecting and resolving collisions.  Currently this class uses penalty forces between rigid bodies.
 * @author kry
 */
public class CollisionProcessor {

    private List<RigidBody> bodies; 
    /**
     * The current contacts that resulted in the last call to process collisions
     */
    public ArrayList<Contact> contacts = new ArrayList<Contact>();
    
    /**
     * Creates this collision processor with the provided set of bodies
     * @param bodies
     */
    public CollisionProcessor( List<RigidBody> bodies ) {
        this.bodies = bodies;
    }
    
    /** keeps track of the time used for collision detection on the last call */
    double collisionDetectTime = 0;
    
    /** keeps track of the time used to solve the LCP based velocity update on the last call */
    double collisionSolveTime = 0;
    
    Vector2d vunconstrained = new Vector2d();
    
    Vector2d vnew = new Vector2d();
    
    Vector2d fa = new Vector2d();
    
    Matrix3d jacobian = new Matrix3d();
    
    //The Mass matrix is a diagonal matrix, so we only need 
    //an array to represent it.
    double[] massMatrix = new double[6];
    double[] massInverse = new double[6];
    
    DenseVector u = new DenseVector(6);
    DenseVector fexternal = new DenseVector(6);    
    /**
     * Processes all collisions 
     * @param dt time step
     */
    public void processCollisions( double dt ) {
        contacts.clear();
        Contact.nextContactIndex = 0;
        long now = System.nanoTime();
        broadPhase();
        collisionDetectTime = ( System.nanoTime() - now ) * 1e-9;
                
        if ( contacts.size() > 0  && doLCP.getValue() ) {
            now = System.nanoTime();
            double bounce = restitution.getValue();
            double mu = friction.getValue();
            //int n = contacts
            
            double lambda=0;
            for (Contact c : contacts) {
                //u = v1,w1,v2,w2 linear and angular velocity
                u.set(0,c.body1.v.x);
                u.set(1,c.body1.v.y);
                u.set(2,c.body1.omega);
                u.set(3,c.body2.v.x);
                u.set(4,c.body2.v.y);
                u.set(5,c.body2.omega);                
                
                fexternal.set(0,c.body1.force.x*c.body1.massLinear);
                fexternal.set(1,c.body1.force.y*c.body1.massLinear);
                fexternal.set(2,c.body1.torque*c.body1.massAngular);
                fexternal.set(3,c.body2.force.x*c.body2.massLinear);
                fexternal.set(4,c.body2.force.y*c.body2.massLinear);
                fexternal.set(5,c.body2.torque*c.body2.massAngular);
                
                fexternal.scale(dt);
                DenseVector b = new DenseVector(2);
                //b = J (u + h * M -1 fext)
                u.add(fexternal);
                c.jacobian.mult(u, b);
                
                //  the projected Guass-Seidel complementarity solver.
                //A = J M -1 J^t h lambda



                //velocity update
                
            	
            	
            	
            	
            	
            	
            	
            	//Compute lambda max.
            	double lambda_max = contact.normal.length();//normal force
            	//Compute mass matrix and it's inverse.
            	//Compute tentative velocity
                vunconstrained.set(contact.body1.v);
                fillMassMatrices(contact.body1,contact.body2);
//                contact.body1.force
            	//Compute mc
            	
            	//Comute lambda
            	
            	//Constrain lambda;
            	constrainLambda(lambda,mu*lambda_max);
            }
            collisionSolveTime = (System.nanoTime() - now) * 1e-9;
        }
    }
    
    private void fillMassMatrices(RigidBody body1, RigidBody body2) {
    	massMatrix[0] = body1.massLinear;
    	massMatrix[1] = 1;
    	massMatrix[2] = 1;
    	massMatrix[3] = body2.massLinear;
    	massMatrix[4] = 1;
    	massMatrix[5] = 1;
    	massInverse[0] = body1.minv;
    	massInverse[1] = 1;
    	massInverse[2] = 1;
    	massInverse[3] = body2.minv;
    	massInverse[4] = 1;
    	massInverse[5] = 1;
    }
    
    /**
     * Returns new constrained lambda.
     * @param lambda
     * @param mu
     * @return
     */
    private double constrainLambda(double lambda, double bound) {
    	if (lambda > bound) {
    		return bound;
    	} else if (lambda < -bound) {
    		return -bound;
    	} else {
    		return lambda;
    	}
    }
    
    /**
     * Checks for collisions between bodies.  Note that you can optionaly implement some broad
     * phase test such as spatial hashing to reduce the n squared body-body tests.
     * Currently this does the naive n squared collision check.
     */
    private void broadPhase() {
        // Naive n squared body test.. might not be that bad for small number of bodies 
        visitID++;
        for ( RigidBody b1 : bodies ) {
            for ( RigidBody b2 : bodies ) {
                if ( b1.index >= b2.index ) continue;
                if ( b1.pinned && b2.pinned ) continue;                
                narrowPhase( b1, b2 );                
            }
        }        
    }
    
    
    /**
     * Checks for collision between boundary blocks on two rigid bodies.
     * @param body1
     * @param body2
     */    
    private void narrowPhase( RigidBody body1, RigidBody body2 ) {
        if ( ! useBVTree.getValue() ) {
            for ( Block b1 : body1.blocks ) {
                for ( Block b2 : body2.blocks ) {
                    processCollision( body1, b1, body2, b2 );
                }
            }
        } else {
        	detectCollision(body1,body2,body1.root,body2.root);
        }
    }
	
    /** This method cheaply computes the square of the distance between two discs.
     * We only need to find which distance is smallest, so there is no point in computing 
     * the square root.
     */ 
    private static double distanceSquared(Disc disk1, Disc disk2) {
    	return disk1.cW.distanceSquared(disk2.cW);
    }
	
    /**
     * This method returns true if there exists a collision and false otherwise. Note that the
     * two block arguments are only assigned if the method returns true.
     * @param left
     * @param right
     * @param b1
     * @param b2
     * @return
     */
    public void detectCollision(RigidBody body1, RigidBody body2, BVNode left, BVNode right) {
    	Disc leftDisk = left.boundingDisc;    	
    	Disc rightDisk = right.boundingDisc;
    	left.updateBoundingDisk(visitID);
    	right.updateBoundingDisk(visitID);
    	if (leftDisk.intersects(rightDisk)) {
    		//If nodes are leaves we can process the collisions directly.
    		if (left.isLeaf()&&right.isLeaf()) {
    			Block b1 = left.leafBlock;
    			Block b2 = right.leafBlock;    			
    			processCollision(body1,b1,body2,b2);
    		} else if (left.isLeaf()) {
    			//Right block is not a leaf.
    			detectCollision(body1,body2,left,right.child1);
    			detectCollision(body1,body2,left,right.child2);
    		} else if (right.isLeaf()) {
    			//Left block is not a leaf.
    			detectCollision(body1,body2,left.child1,right);
    			detectCollision(body1,body2,left.child2,right);
    		} else {
    			detectCollision(body1,body2,left.child1,right.child1);
    			detectCollision(body1,body2,left.child1,right.child2);
    			detectCollision(body1,body2,left.child2,right.child1);
    			detectCollision(body1,body2,left.child2,right.child2);    			
    		}
    	}
    }	    

    
    /** 
     * The visitID is used to tag boundary volumes that are visited in 
     * a given collision detection pass.  A boolean flag
     * at each BVNode would be inefficient as we would need
     * to visit all nodes to clear the flag.
     */
    int visitID = 0;
    
    /**
     * Resets the state of the collision processor by clearing all
     * currently identified contacts, and reseting the visitID for
     * tracking the bounding volumes used
     */
    public void reset() {
        contacts.clear();
        Contact.nextContactIndex = 0;
        visitID = 0;            
    }
    
    // some working variables for processing collisions
    private Point2d tmp1 = new Point2d();
    private Point2d tmp2 = new Point2d();
    private Point2d contactW = new Point2d();
    private Vector2d force = new Vector2d();
    private Vector2d contactV1 = new Vector2d();
    private Vector2d contactV2 = new Vector2d();
    private Vector2d relativeVelocity = new Vector2d();
    private Vector2d normal = new Vector2d();
        
    /**
     * Processes a collision between two bodies for two given blocks that are colliding.
     * Currently this implements a penalty force
     * @param body1
     * @param b1
     * @param body2
     * @param b2
     */
    private void processCollision( RigidBody body1, Block b1, RigidBody body2, Block b2) {        
        double k = contactSpringStiffness.getValue();
        double c1 = contactSpringDamping.getValue();
        double threshold = separationVelocityThreshold.getValue();
        boolean useSpring = enableContactSpring.getValue();
        boolean useDamping = enableContactDamping.getValue();
        body1.transformB2W.transform( b1.pB, tmp1 );
        body2.transformB2W.transform( b2.pB, tmp2 );
        double distance = tmp1.distance(tmp2);
        if ( distance < Block.radius * 2 ) {
        	//log.info("Collision : " + tmp1 + "," + tmp2);
            // contact point at halfway between points 
            // NOTE: this assumes that the two blocks have the same radius!
            contactW.interpolate( tmp1, tmp2, .5 );
            // contact normal
            normal.sub( tmp2, tmp1 );
            normal.normalize();
            // create the contact
            Contact contact = new Contact( body1, body2, contactW, normal);
            // simple option... add to contact list...
            contacts.add( contact );
            if ( ! doLCP.getValue() ) {
                // compute relative body velocity at contact point
                body1.getSpatialVelocity( contactW, contactV1 );
                body2.getSpatialVelocity( contactW, contactV2 );
                relativeVelocity.sub( contactV1, contactV2 );
                if ( -relativeVelocity.dot( normal ) < threshold ) {
                    if ( useSpring ) {
                        // spring force
                        double interpenetration = distance - Block.radius * 2; // a negative quantity
                        force.scale( -interpenetration * k, normal );
                        body2.applyContactForceW(contactW, force);
                        force.scale(-1);
                        body1.applyContactForceW(contactW, force);
                    }
                    if ( useDamping ) {
                        // spring damping forces!
                        // vertical
                        force.scale( relativeVelocity.dot(normal) * c1, normal );                    
                        body2.applyContactForceW( contactW, force );
                        force.scale(-1);
                        body1.applyContactForceW( contactW, force );
                    }
                }
            }
        }
    }
   
    /** Stiffness of the contact penalty spring */
    private DoubleParameter contactSpringStiffness = new DoubleParameter("penalty contact stiffness", 1e3, 1, 1e5 );
    
    /** Viscous damping coefficient for the contact penalty spring */
    private DoubleParameter contactSpringDamping = new DoubleParameter("penalty contact damping", 10, 1, 1e4 );
    
    /** Threshold for the relative velocity in the normal direction, for determining if spring force will be applied. */
    private DoubleParameter separationVelocityThreshold = new DoubleParameter( "penalty separation velocity threshold (controls bounce)", 1e-9, 1e-9, 1e3 );
    
    /** Enables the contact penalty spring */
    private BooleanParameter enableContactSpring = new BooleanParameter("enable penalty contact spring", true );
    
    /** Enables damping of the contact penalty spring */
    private BooleanParameter enableContactDamping = new BooleanParameter("enable penalty contact damping", true );
    
    /** Restitution parameter for contact constraints */
    public DoubleParameter restitution = new DoubleParameter( "restitution (bounce)", 0, 0, 1 );
    
    /** Coulomb friction coefficient for contact constraint */
    public DoubleParameter friction = new DoubleParameter("Coulomb friction", 0.33, 0, 2 );
    
    /** Number of iterations to use in projected Gauss Seidel solve */
    public IntParameter iterations = new IntParameter("iterations for GS solve", 10, 1, 500);
    
    /** Flag for switching between penalty based contact and contact constraints */
    private BooleanParameter doLCP = new BooleanParameter( "do LCP solve", false );
    
    /** Flag for enabling the use of hierarchical collision detection for body pairs */
    private BooleanParameter useBVTree = new BooleanParameter( "use BVTree", false );
    
    /**
     * @return controls for the collision processor
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.setBorder( new TitledBorder("Collision Processing Controls") );
        vfp.add( useBVTree.getControls() );
        vfp.add( doLCP.getControls() );
        vfp.add( iterations.getSliderControls() );
        vfp.add( restitution.getSliderControls(false) );
        vfp.add( friction.getSliderControls(false) );
        
        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder( new TitledBorder("penalty method controls") );
        vfp2.add( contactSpringStiffness.getSliderControls(true) );
        vfp2.add( contactSpringDamping.getSliderControls(true) );
        vfp2.add( separationVelocityThreshold.getSliderControls( true ) );
        vfp2.add( enableContactDamping.getControls() );
        vfp2.add( enableContactSpring.getControls() );
        
        CollapsiblePanel cp = new CollapsiblePanel(vfp2.getPanel());
        cp.collapse();
        vfp.add( cp );        
        return vfp.getPanel();
    }
    
}
