package comp559.a2;

import javax.media.opengl.GL;
import javax.media.opengl.GL2;
import javax.media.opengl.GLAutoDrawable;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

/**
 * Implementation of a contact constraint.
 * @author kry
 */
public class Contact {

    /** Next available contact index, used for determining which rows of the jacobian a contact uses */
    static public int nextContactIndex = 0;
    
    /** Index of this contact, determines its rows in the jacobian */
    int index;
    
    /** First RigidBody in contact */
    RigidBody body1;
    
    /** Second RigidBody in contact */
    RigidBody body2;
    
    /** Contact normal in world coordinates */
    Vector2d normal = new Vector2d();
    
    /** Position of contact point in world coordinates */
    Point2d contactW = new Point2d();
    
    /** Contact jacobian*/
    double jacobian[][]=new double[2][6];
    Vector2d ri1 = new Vector2d();
    Vector2d ri2 = new Vector2d();
    private Vector2d tangent = new Vector2d();
    
    /**
     * Creates a new contact, and assigns it an index
     * @param body1
     * @param body2
     * @param contactW
     * @param normal
     */
    public Contact( RigidBody body1, RigidBody body2, Point2d contactW, Vector2d normal ) {
        this.body1 = body1;
        this.body2 = body2;
        this.contactW.set( contactW );
        this.normal.set( normal );
        this.tangent.set(-normal.y,normal.x);
        index = nextContactIndex++;
        ri1.set(contactW.x-body1.x.x,contactW.y-body1.x.y);
        ri2.set(contactW.x-body2.x.x,contactW.y-body2.x.y);        
        
        jacobian[0][0]=-normal.x;
        jacobian[0][1]=-normal.y;
        jacobian[0][2]=-(ri1.x*normal.y-normal.x*ri1.y);
        jacobian[0][3]=normal.x;
        jacobian[0][4]=normal.y;
        jacobian[0][5]=ri2.x*normal.y-normal.x*ri2.y;
        
        jacobian[1][0]=-tangent.x;
        jacobian[1][1]=-tangent.y;
        jacobian[1][2]=-(ri1.x*tangent.y-tangent.x*ri1.y);
        jacobian[1][3]=tangent.x;
        jacobian[1][4]=tangent.y;
        jacobian[1][5]=ri2.x*tangent.y-tangent.x*ri2.y;        
    }
    
    /**
     * Draws the contact points
     * @param drawable
     */
    public void display( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        gl.glPointSize(3);
        gl.glColor3f(.7f,0,0);
        gl.glBegin( GL.GL_POINTS );
        gl.glVertex2d(contactW.x, contactW.y);
        gl.glEnd();
    }
    
    /**
     * Draws the connections between bodies to visualize the 
     * the adjacency structure of the matrix as a graph.
     * @param drawable
     */
    public void displayConnection( GLAutoDrawable drawable ) {
        GL2 gl = drawable.getGL().getGL2();
        // draw a line between the two bodies but only if they're both not pinned
        if ( !body1.pinned && ! body2.pinned ) {
            gl.glLineWidth(2);
            gl.glColor4f(0,.3f,0, 0.5f);
            gl.glBegin( GL.GL_LINES );
            gl.glVertex2d(body1.x.x, body1.x.y);
            gl.glVertex2d(body2.x.x, body2.x.y);
            gl.glEnd();
        }
    }
    
}
