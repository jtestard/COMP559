/**
 * 
 */
package comp559.a2;

/**
 * @author julestestard
 *
 */
public class CollisionDetectionHelper {
	private BVNode smallest;
	private BVNode largest;
	private BVNode closest;
	private BVNode furthest;
	private int visitID;
	private CollisionProcessor cp;
	
	/**
	 * @param left
	 * @param right
	 */
	public CollisionDetectionHelper(CollisionProcessor cp) {
		this.cp = cp;
		visitID=  cp.visitID;
	}
	
    /** This method cheaply computes the square of the distance between two discs.
     * We only need to find which distance is smallest, so there is no point in computing 
     * the square root.
     */ 
    private static double distanceSquared(Disc disk1, Disc disk2) {
    	return disk1.cW.x*disk2.cW.x + disk1.cW.y*disk2.cW.y;
    }
    
    
    private void largest(BVNode left,BVNode right) {
    	if (left.boundingDisc.r > right.boundingDisc.r) {
    		largest = left;
    		smallest = right;
    	} else {
    		largest = right;
    		smallest = left;
    	}
    }
    
    private void closest(BVNode node, BVNode leftchild, BVNode rightchild) {
    	if (distanceSquared(node.boundingDisc,leftchild.boundingDisc)>
    		distanceSquared(node.boundingDisc,rightchild.boundingDisc)){
    		//rightchild is closest to the node
    		closest = rightchild;
    		furthest = leftchild;
    	} else {
    		//leftchild is closest to the node
    		closest = leftchild;
    		furthest = rightchild;
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
    public void detectCollision(BVNode left, BVNode right) {
    	Disc leftDisk = left.boundingDisc;    	
    	Disc rightDisk = right.boundingDisc;
    	if (left.alreadyVisited(visitID) && right.alreadyVisited(visitID))
    		return;
    	left.updateBoundingDisk(visitID);
    	right.updateBoundingDisk(visitID);	
    	if (leftDisk.intersects(rightDisk)) {
    		//If nodes are leaves we can process the collisions directly.
    		if (left.isLeaf()&&right.isLeaf()) {
    			Block b1 = left.leafBlock;
    			Block b2 = right.leafBlock;
    		}
    		if (left.isLeaf()) {
    			//Right block is not a leaf.
    			right.child1.boundingDisc.updatecW();
    			right.child2.boundingDisc.updatecW();
    			closest(left,right.child1,right.child2);
    		} else if (right.isLeaf()) {
    			//Left block is not a leaf.
    			left.child1.boundingDisc.updatecW();
    			left.child2.boundingDisc.updatecW();    			
    			closest(right,left.child1,left.child2);	
    		} else {
    			largest(left,right);
    			largest.child1.updateBoundingDisk(visitID);
    			largest.child2.updateBoundingDisk(visitID);
    			closest(smallest,largest.child1,largest.child2);
    		}
    	}
    }	
}
