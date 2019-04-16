/**
 * 
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
