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

import java.awt.Point;

import buoy.event.WidgetMouseEvent;
import artofillusion.Camera;
import artofillusion.LayoutWindow;
import artofillusion.Scene;
import artofillusion.UndoRecord;
import artofillusion.ViewerCanvas;
import artofillusion.animation.PositionTrack;
import artofillusion.animation.RotationTrack;
import artofillusion.math.CoordinateSystem;
import artofillusion.math.Vec3;
import artofillusion.object.ObjectInfo;
import artofillusion.ui.EditingTool;
import artofillusion.ui.EditingWindow;
import artofillusion.ui.Translate;

/**
 * Class to create a tracker, which is a special null object
 * that can be set to track the movement of another object in
 * the scene.
 * @author Timothy Fish
 *
 */
public class CreateFanTool extends EditingTool {

  static int counter = 1;
  private Point clickPoint;
  private ObjectInfo objInfo;


  /**
   * constructor
   * @param win
   */
  public CreateFanTool(EditingWindow win) {
    super(win);

    // Name is the icon name without an extension.
    initButton("ClothMaker:fanTool"); 
  }

  @Override
  public void activate()
  {
    super.activate();
    theWindow.setHelpText(Translate.text("Create Fan: Click to place fan."));
  }

  @Override
  public int whichClicks()
  {
    return ALL_CLICKS;
  }

  @Override
  public String getToolTipText()
  {
    return Translate.text("Create Fan");
  }

  /**
   * Adds the Fan object to the scene where the user clicks.
   */
  public void mousePressed(WidgetMouseEvent e, ViewerCanvas view)
  {
    clickPoint = e.getPoint();

    UndoRecord undo = new UndoRecord(theWindow, false);
    // Create the Tracker
    Scene theScene = ((LayoutWindow) theWindow).getScene();
    objInfo = new ObjectInfo(new Fan(theScene, Fan.DEFAULT_WIND_FORCE, Fan.DEFAULT_WIND_FALLOFF), new CoordinateSystem(), "Fan "+counter);

    counter++;
    objInfo.addTrack(new PositionTrack(objInfo), 0);
    objInfo.addTrack(new RotationTrack(objInfo), 1);

    int sel[] = ((LayoutWindow) theWindow).getSelectedIndices();
    ((LayoutWindow) theWindow).addObject(objInfo, undo);
    undo.addCommand(UndoRecord.SET_SCENE_SELECTION, new Object [] {sel});
    theWindow.setUndoRecord(undo);
    ((LayoutWindow) theWindow).setSelection(theScene.getNumObjects()-1);

    // Determine position of the Tracker.
    Camera cam = view.getCamera();
    Vec3 orig;

    orig = cam.convertScreenToWorld(clickPoint, Camera.DEFAULT_DISTANCE_TO_SCREEN);

    // Update the position, and redraw the display.    
    objInfo.getCoords().setOrigin(orig);
    objInfo.clearCachedMeshes();
    theWindow.updateImage();

  }

  @Override
  public void mouseReleased(WidgetMouseEvent e, ViewerCanvas view)
  {
    objInfo = null;
    theWindow.setModified();
  }

}
