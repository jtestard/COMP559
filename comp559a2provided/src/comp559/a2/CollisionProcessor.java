package comp559.a2;

import java.util.ArrayList;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

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
    
    /** Used during collision detection*/
    public double distances[];    
    
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
        distances = new double[4];
    }
    
    /** keeps track of the time used for collision detection on the last call */
    double collisionDetectTime = 0;
    
    /** keeps track of the time used to solve the LCP based velocity update on the last call */
    double collisionSolveTime = 0;
    
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
            // TODO: Compute velocity update with iterative solve of contact constraint matrix.

            
            collisionSolveTime = (System.nanoTime() - now) * 1e-9;
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
    
    private Block b1=null,b2=null;
    
    /**
     * Checks for collision between boundary blocks on two rigid bodies.
     * TODO: This needs to be improved as the n-squared block test is too slow!
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
        	System.out.println("Next Collision detection phase : " + body1.root.boundingDisc.cW + "," + body2.root.boundingDisc.cW);
        	if (hasCollision(body1.root,body2.root)) {
        		processCollision(body1,b1,body2,b2);
        	}
        }
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
    public boolean hasCollision(BVNode left, BVNode right) {
    	Disc leftDisk = left.boundingDisc;    	
    	Disc rightDisk = right.boundingDisc;
    	leftDisk.updatecW();
    	rightDisk.updatecW();
    	if (leftDisk.intersects(rightDisk)) {
    		//If nodes are leaves we can process the collisions directly.
    		if (left.isLeaf()&&right.isLeaf()) {
    			b1 = left.leafBlock;
    			b2 = right.leafBlock;
    			return true;
    		} else if (left.isLeaf()) {
    			//Right block is not a leaf.
    			Disc rightleftDisk = right.child1.boundingDisc;
    			if (right.child2==null) {
    				return hasCollision(left,right.child1);
    			}
    			Disc rightrightDisk = right.child2.boundingDisc;
    			rightleftDisk.updatecW();
    			rightrightDisk.updatecW();
    			distances[0] = distanceSquared(leftDisk,rightleftDisk);
    			distances[1] = distanceSquared(leftDisk,rightrightDisk);
    			switch (smallestDistance(true)) {
    			case 0:
    				if (leftDisk.intersects(rightleftDisk))
    					return hasCollision(left,right.child1);
    				break;
    			default:
    				if (leftDisk.intersects(rightrightDisk))
    					return hasCollision(left,right.child2);
    			}
    		} else if (right.isLeaf()) {
    			//Left block is not a leaf.
    			Disc leftleftDisk = left.child1.boundingDisc;
    			if (left.child2==null) {
    				return hasCollision(left.child1,right);
    			}
    			Disc leftrightDisk = left.child2.boundingDisc;
    			leftleftDisk.updatecW();
    			leftrightDisk.updatecW();
    			distances[0] = distanceSquared(rightDisk,leftleftDisk);
    			distances[1] = distanceSquared(rightDisk,leftrightDisk);
    			switch (smallestDistance(true)) {
    			case 0:
    				if (rightDisk.intersects(leftleftDisk))
    					return hasCollision(left.child1,right);
    				break;
    			default:
    				if (rightDisk.intersects(leftrightDisk))
    					return hasCollision(left.child2,right);
    			}
    		} else {
        		//Else we need to do a narrower search.
        		Disc leftleftDisk = left.child1.boundingDisc;
        		Disc leftrightDisk = left.child2.boundingDisc;
        		Disc rightleftDisk = right.child1.boundingDisc;
        		Disc rightrightDisk = right.child2.boundingDisc;
        		leftleftDisk.updatecW();
        		leftrightDisk.updatecW();
        		rightleftDisk.updatecW();
        		rightrightDisk.updatecW();        		
        		//Compute the distances between the centers of the discs in world coordinates.
        		distances[0] = distanceSquared(leftleftDisk,rightleftDisk);
        		distances[1] = distanceSquared(leftleftDisk,rightrightDisk);
        		distances[2] = distanceSquared(leftrightDisk,rightleftDisk);
        		distances[3] = distanceSquared(leftrightDisk,rightrightDisk);
        		int smallest = smallestDistance(false);
        		System.out.println("smallest distance :" + smallest);
        		switch(smallest) {
        		case 0:
        			if (leftleftDisk.intersects(rightleftDisk))
        				return hasCollision(left.child1,right.child1);
        			break;
        		case 1:
        			if (leftleftDisk.intersects(rightrightDisk))
        				return hasCollision(left.child1,right.child2);
        			break;
        		case 2:
        			if (leftrightDisk.intersects(rightleftDisk))
        				return hasCollision(left.child2,right.child1);
        			break;
        		default :
        			if (leftrightDisk.intersects(rightrightDisk))
        				return hasCollision(left.child2,right.child2);
        		}    			
    		}
    	}
    	System.out.println("Disks not intersecting : " + leftDisk.cW + "," + leftDisk.r + "," + leftDisk.body.index + ";"+ rightDisk.cW + "," + rightDisk.r + "," +  rightDisk.body.index);
    	return false;
    }
    
    /** This method cheaply computes the square of the distance between two discs.
     * We only need to find which distance is smallest, so there is no point in computing 
     * the square root.
     */ 
    public double distanceSquared(Disc disk1, Disc disk2) {
    	return disk1.cW.x*disk2.cW.x + disk1.cW.y*disk2.cW.y;
    }
    
    /**
     * Small function that returns the index in the distances array of the smallest distance
     * between two disk when solving for collisions.
     * @return
     */
    public int smallestDistance(boolean size2) {
    	// We have only two distances to cover (only one node isn't a leaf).
    	if (size2) {
    		return (distances[0]>distances[1])?1:0;
    	} else {
        	int first, second;
        	double d1,d2;
        	if (distances[0]>distances[1]) {
        		d1 = distances[1];
        		first = 1;
        	} else {
        		d1 = distances[0];
        		first = 0;
        	}
        	if (distances[2]>distances[3]) {
        		d2 = distances[3];
        		second = 3;
        	} else {
        		d2 = distances[2];
        		second = 2;
        	}
        	return (d1>d2)?second:first;    		
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
    private void processCollision( RigidBody body1, Block b1, RigidBody body2, Block b2 ) {        
        double k = contactSpringStiffness.getValue();
        double c1 = contactSpringDamping.getValue();
        double threshold = separationVelocityThreshold.getValue();
        boolean useSpring = enableContactSpring.getValue();
        boolean useDamping = enableContactDamping.getValue();
        body1.transformB2W.transform( b1.pB, tmp1 );
        body2.transformB2W.transform( b2.pB, tmp2 );
        double distance = tmp1.distance(tmp2);
        if ( distance < Block.radius * 2 ) {
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
