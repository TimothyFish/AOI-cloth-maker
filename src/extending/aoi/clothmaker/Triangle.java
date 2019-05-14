/**
    Cloth Maker Plugin from Chapter 10 of the book "Extending Art of Illusion: Scripting for 3D Artists"
    Copyright (C) 2019, 2011  Timothy Fish

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>
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
