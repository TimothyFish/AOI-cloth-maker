/**
 * 
 */
package extending.aoi.clothmaker;

import artofillusion.LayoutWindow;
import artofillusion.object.ObjectInfo;
import buoy.widget.BMenuItem;

/**
 * Thread that activates/deactivates menu items based on whether a
 * cloth object is selected.
 * 
 * @author Timothy Fish
 *
 */
public class ClothMenuItemActivator extends Thread {

  private BMenuItem theConvertMenuItem;
  private BMenuItem theGenerateMenuItem;
  private LayoutWindow theLayout;
  private CollisionDetector theCollisionDetector;
  private static final long inverseRefreshRate = 250; // refresh menuItem every quarter-second

  /**
   * Constructor
   * @param convertMenuItem Non-null menu item
   * @param generateMenuItem Non-null menu item 
   * @param layout non-null layout window
   */
  public ClothMenuItemActivator(LayoutWindow layout, BMenuItem convertMenuItem, BMenuItem generateMenuItem) {
    theConvertMenuItem = convertMenuItem;
    theGenerateMenuItem = generateMenuItem;
    theLayout = layout;
    theCollisionDetector = new CollisionDetector(layout.getScene());

    theConvertMenuItem.setEnabled(false);
    theGenerateMenuItem.setEnabled(false);
  }

  /**
   * Method that does checks the selection status.
   */
  @Override
  public void run() {
    while((theConvertMenuItem != null || theGenerateMenuItem != null) && theLayout != null){ // do while menu item exists
      // Check for object selection
      if(theGenerateMenuItem != null) {
        if(oneClothObjectSelected()){
          theGenerateMenuItem.setEnabled(true);
        }
        else {
          theGenerateMenuItem.setEnabled(false);        
        }

      }

      if(theConvertMenuItem != null) {
        if(oneNonClothObjectSelected()) {
          theConvertMenuItem.setEnabled(true);
        }
        else{
          theConvertMenuItem.setEnabled(false);
        }        
      }

      try {
        // Sleep most of the time so processor can do other things
        // Updating every quarter-second is more than sufficient since
        // the user can't select the menu item that quickly after object selection.
        Thread.sleep(inverseRefreshRate );
      } catch (InterruptedException e) {
        // This should never happen.
        e.printStackTrace();
      }
    }    
  }

  /**
   * Returns true when one object is selected that can be converted to cloth.
   * @return
   */
  private boolean oneNonClothObjectSelected() {
    int numSelected = theLayout.getSelectedObjects().size();
    if(numSelected == 1) {
      ObjectInfo ref = (ObjectInfo)theLayout.getSelectedObjects().toArray()[0];
      if(theCollisionDetector.isSpecial(ref) || (ref.getObject() instanceof Cloth)){
        return false;  
      }
      else {
        return true;
      }
    }

    return false;
  }

  /**
   * Returns true when a single Cloth object is selected.
   * @return
   */
  private boolean oneClothObjectSelected() {
    int numSelected = theLayout.getSelectedObjects().size();
    if(numSelected == 1) {
      ObjectInfo ref = (ObjectInfo)theLayout.getSelectedObjects().toArray()[0];
      if(ref.getObject() instanceof Cloth){
        return true;  
      }
    }

    return false;
  }

}
