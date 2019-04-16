/**
 * 
 */
package extending.aoi.clothmaker;

import java.util.Vector;

import artofillusion.math.Vec3;

/**
 * Object to represent a mass at the end of a spring used in cloth simulation.
 * Each vertex of the cloth object is considered a Mass and Springs are used
 * to connect them.
 * 
 * @author Timothy Fish
 *
 */
public class Mass {
	/**
	 * Constructor
	 * @param pos
	 * @param index
	 * @param mass
	 * @param vel
	 */
  public Mass(Vec3 pos, int index, double mass, Vec3 vel) {
    vertexIndex = index;
    vertexMass = mass;
    position = pos;
    velocity = vel;
    springRef = new Vector<Spring>();
  }
  
  /**
   * Constructor
   * @param mass
   */
  public Mass(Mass mass) {
    vertexIndex = mass.vertexIndex;
    vertexMass = mass.vertexMass;
    position = mass.position;
    velocity = mass.velocity;
    springRef = new Vector<Spring>(mass.springRef);
  }
  
  /**
   * Create another Mass like this one.
   * @return
   */
  public Mass duplicate() {
    return new Mass(this);
  }
  
  /**
   * Put the Mass at the end of a Spring
   * @param S
   */
  public void connectToSpring(Spring S) {
    
    for(Spring spg: springRef) {
      if(spg == S) {
        return;
      }
    }
    
    springRef.addElement(S);
  }
  
  /**
   * Remove the Mass from the Spring
   * @param S
   */
  public void disconnectFromSpring(Spring S) {
    
        springRef.remove(S);
  }
  
  /**
   * Get the location of the Mass/Vertex
   * @return
   */
  public Vec3 getPosition() {
    return position;
  }
  
  /**
   * Get the Springs the Mass is connected to.
   * @return
   */
  public Vector<Spring> getSprings(){
    return springRef;
  }
  
  private Vec3 position;
  private int vertexIndex; // index of this mass in the mesh
  private double vertexMass; 
  private Vec3 velocity; // velocity now, whatever "now" means
  private Vector<Spring> springRef; // references to the springs that are connect to this mass
  
  /**
   * Get how heavy the mass is.
   * @return
   */
	public double getWeight() {
		return vertexMass;
	}

	/**
	 * Get the direction and speed the mass is moving.
	 * @return
	 */
	public Vec3 getVelocity() {
		return velocity;
	}

	/**
	 * Set the position of the Mass.
	 * @param v
	 */
	public void setPosition(Vec3 v) {
		position = v;		
	}
	
	/**
	 * Set the direction and speed of the mass.
	 * @param v
	 */
	public void setVelocity(Vec3 v) {
		velocity = v;
	}
}
