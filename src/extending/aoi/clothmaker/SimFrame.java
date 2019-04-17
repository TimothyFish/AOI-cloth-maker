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

/**
 * Object for storing a version of the cloth in a list from 
 * which we can retreave the cloth at a given frame number.
 * @author Timothy Fish
 *
 */
public class SimFrame{
	public int frameNumber;
	public double time;
	public Cloth M;

	/** 
	 * Constructor
	 * @param frame
	 * @param M
	 */
	public SimFrame(int frame, Cloth M) {
    this.frameNumber = frame;
    this.time = (double)frame*(1.0/ClothSimEditorWindow.getFPS());
    this.M = (Cloth) M.duplicate();
  }
}