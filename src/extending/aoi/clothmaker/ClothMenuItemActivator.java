/**
    Cloth Maker Plugin from Chapter 10 of the book "Extending Art of Illusion: Scripting for 3D Artists"
    Copyright (C) 2019, 2011  Timothy Fish
    Changes copyright (C) 2019 by Maksim Khramov

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
public class ClothMenuItemActivator implements Runnable {

  private final BMenuItem theConvertMenuItem;
  private final BMenuItem theCopyToTriMenuItem;
  private final BMenuItem theGenerateMenuItem;
  private final LayoutWindow theLayout;


  /**
   * Constructor
   * @param convertMenuItem Non-null menu item
   * @param generateMenuItem Non-null menu item 
   * @param layout non-null layout window
   */
  public ClothMenuItemActivator(LayoutWindow layout, 
  		                          BMenuItem convertMenuItem, 
  		                          BMenuItem copyToTriMenuItem,
  		                          BMenuItem generateMenuItem) {
    theConvertMenuItem = convertMenuItem;
    theCopyToTriMenuItem = copyToTriMenuItem;
    theGenerateMenuItem = generateMenuItem;
    theLayout = layout;

    theConvertMenuItem.setEnabled(false);
    theCopyToTriMenuItem.setEnabled(false);
    theGenerateMenuItem.setEnabled(false);
  }

  /**
   * Method that does checks the selection status.
   */
  @Override
  public void run() {

    boolean singleClothObjectSelected = oneClothObjectSelected();
    theGenerateMenuItem.setEnabled(singleClothObjectSelected);
    theCopyToTriMenuItem.setEnabled(singleClothObjectSelected);

    theConvertMenuItem.setEnabled(oneNonClothObjectSelected());
  }

  /**
   * Returns true when one object is selected that can be converted to cloth.
   * @return
   */
  private boolean oneNonClothObjectSelected() {
    ObjectInfo[] selection = theLayout.getSelectedObjects().toArray(new ObjectInfo[0]);
    if(selection.length != 1) return false;

    ObjectInfo ref = selection[0];
    return !(CollisionDetector.isSpecial(ref) || (ref.getObject() instanceof Cloth));
    
  }

  /**
   * Returns true when a single Cloth object is selected.
   * @return
   */
  private boolean oneClothObjectSelected() {
    ObjectInfo[] selection = theLayout.getSelectedObjects().toArray(new ObjectInfo[0]);
    if(selection.length != 1) return false;

    return selection[0].getObject() instanceof Cloth;
  }

}
