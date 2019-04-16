/**
 * 
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