/**
    Cloth Maker Plugin from Chapter 10 of the book "Extending Art of Illusion: Scripting 3D Scene Creation"
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
 * Provides a connector that pulls two Masses toward each other during
 * a cloth simulation.
 * @author Timothy Fish
 *
 */
public class Spring {
	/**
	 * Constructor
	 * @param A
	 * @param B
	 * @param restLen
	 * @param spgConst
	 */
  public Spring(Mass A, Mass B, double restLen, double spgConst) {
    massA = A;
    massB = B;
    restingLength = restLen;
    k = spgConst;
    
    massA.connectToSpring(this);
    massB.connectToSpring(this);
  }
  
  /** 
   * Constructor
   * @param obj
   */
  public Spring(Spring obj) {
    massA = obj.massA;
    massB = obj.massB;
    restingLength = obj.restingLength;
  }

  /**
   * The force pulling on the two Masses
   * @return
   */
  Vec3 getForce() {
    Vec3 F = new Vec3();
    
    Vec3 p = massA.getPosition();
    Vec3 q = massB.getPosition();
    Vec3 d = q.minus(p);
    double x = d.length();
    d.normalize();
    F = F.plus(d.times(-k * (restingLength - x)));
    
    return F;
  
  }
  
  private Mass massA; // first mass spring is connected to
  private Mass massB; // second mass spring is connected to
  private double restingLength; // length of the spring when no force is applied
  private double k;
	
  /**
   * Get one of the masses
   * @return
   */
  public Mass getMassA() {
		return massA;
	}

  /**
   * Get one of the masses
   * @return
   */
	public Mass getMassB() {
		return massB;
	}

	/**
	 * Get the resting lenghth of the Spring
	 * @return
	 */
	public double getRestLength() {
		return restingLength;
	}
}
