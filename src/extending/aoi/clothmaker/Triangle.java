/**
 *
 */

package extending.aoi.clothmaker;

import artofillusion.math.Vec3;

/**
 * A simple class to hold the points of a triangle and to 
 * compute their normal.
 * @author Timothy Fish
 *
 */
public class Triangle {
  private Vec3 p0;
  private Vec3 p1;
  private Vec3 p2;
  private Vec3 normal;

  /**
   * Constructor
   * @param p0
   * @param p1
   * @param p2
   */
  public Triangle(final Vec3 p0, final Vec3 p1, final Vec3 p2){
    this.p0 = p0;
    this.p1 = p1;
    this.p2 = p2;
    this.setNormal();
  }

  /**
   * Find the Normal of the Triangle
   */
  private void setNormal() {
    Vec3 n = p1.minus(p0).cross(p2.minus(p0));
    this.normal = n.times(1/Math.sqrt(n.dot(n)));
    this.normal.normalize();
  }

  /**
   * Return the normal of the Triangle
   * @return
   */
  public Vec3 getNormal(){
    return normal;
  }

  /**
   * Set one of vertex positions
   * @param p0
   */
  public void setP0(Vec3 p0) {
    this.p0 = p0;
    this.setNormal();
  }

  /**
   * Get one of vertex positions
   * @return
   */
  public Vec3 getP0() {
    return p0;
  }

  /**
   * Set one of vertex positions
   * @param p1
   */
  public void setP1(Vec3 p1) {
    this.p1 = p1;
    this.setNormal();
  }

  /**
   * Get one of vertex positions
   * @return
   */
  public Vec3 getP1() {
    return p1;
  }

  /**
   * Set one of vertex positions
   * @param p2
   */
  public void setP2(Vec3 p2) {
    this.p2 = p2;
    this.setNormal();
  }

  /**
   * Get one of vertex positions
   * @return
   */
  public Vec3 getP2() {
    return p2;
  }


}
