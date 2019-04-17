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

import artofillusion.ui.EditingTool;
import artofillusion.ui.EditingWindow;
import artofillusion.ui.Translate;

/**
 * Tool used within the ClothSimEditorWindow to select vertices.
 * It is necessary to have this tool because the other tool buttons
 * set it in a state where the camera is moved but we are unable to
 * select anything in that mode.
 * 
 * @author Timothy Fish
 *
 */
public class SelectTool extends EditingTool {
  /**
   * Constructor
   * @param win
   */
  public SelectTool(EditingWindow win) {
    super(win);

    initButton("ClothMaker:selectTool");
  }

  @Override
  public void activate()
  {
    super.activate();
    theWindow.setHelpText(Translate.text("Select Vertices"));
  }

  @Override
  public int whichClicks()
  {
    return HANDLE_CLICKS;
  }

  @Override
  public boolean allowSelectionChanges()
  {
    return true;
  }

  @Override
  public String getToolTipText()
  {
    return Translate.text("Select Vertices");
  }

}
