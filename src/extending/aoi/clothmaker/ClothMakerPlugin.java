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

import java.util.Collection;

import artofillusion.LayoutWindow;
import artofillusion.Plugin;
import artofillusion.UndoRecord;
import artofillusion.animation.Score;
import artofillusion.math.CoordinateSystem;
import artofillusion.object.Object3D;
import artofillusion.object.ObjectInfo;
import artofillusion.ui.ComponentsDialog;
import artofillusion.ui.EditingWindow;
import artofillusion.ui.ToolPalette;
import artofillusion.ui.Translate;
import artofillusion.ui.ValueField;
import buoy.widget.BMenu;
import buoy.widget.BMenuItem;
import buoy.widget.Widget;


/**
 * This class is the entry point for Cloth Simulation. It sets up the
 * Art of Illusion environment with the classes, windows, and menu items
 * that are needed to simulate cloth.
 * 
 * @author Timothy Fish
 *
 */
public class ClothMakerPlugin implements Plugin {
  public static final double DEFAULT_MESH_TOLERANCE = 0.1;
  public static final double DEFAULT_MASS_DISTANCE = 0.2;
  public static final double DEFAULT_SPRING_CONST = 7.1;
  public static final double DEFAULT_DAMPING_CONST = 2.8;
  public static final double DEFAULT_COLLISION_DISTANCE = 0.025;
  public static final double DEFAULT_VERTEX_MASS = 0.5;
  public static final double DEFAULT_GRAVITY = -0.08;
  public static final double DEFAULT_START_TIME = 0.0;
  public static final int DEFAULT_GRAVITY_AXIS = ClothDistortion.Y_AXIS;
  public static final int DEFAULT_DRAPE_FRAMES = 0;
  public static final int DEFAULT_SIM_FRAMES = 30;
  public static final int DEFAULT_FRAMES_PER_SECOND = 30;
  public static final double DEFAULT_FLOOR = 0;
  public static final boolean DEFAULT_SELF_COLLISION = true;
  public static final boolean DEFAULT_FLOOR_COLLISION = false;
  public static final double DEFAULT_WIND_MAGNITUDE = 0.02;
  public static final int DEFAULT_SUBFRAMES = 10;
  private LayoutWindow layout;
  private ToolPalette toolPalette;
  private CreateFanTool theFanTool;
  private double meshTolerance;
  private double massDistance;
  private double springConstant;
  private double dampingConstant;
  private double collisionDistance;
  private int counter;

  /**
   * Constructor
   */
  public ClothMakerPlugin() {
    layout = null;
    meshTolerance = DEFAULT_MESH_TOLERANCE;
    massDistance = DEFAULT_MASS_DISTANCE;
    springConstant = DEFAULT_SPRING_CONST;
    dampingConstant = DEFAULT_DAMPING_CONST;
    collisionDistance = DEFAULT_COLLISION_DISTANCE;
    counter = 1;
  }

  @Override
  public void processMessage(int message, Object[] args) {
    switch (message) {
    case Plugin.SCENE_WINDOW_CREATED:
      layout = (LayoutWindow) args[0];
      toolPalette = layout.getToolPalette();
      theFanTool = new CreateFanTool(layout);

      // Add the Fan button to the Tool Palette
      toolPalette.addTool(theFanTool);

      BMenu objectMenu = layout.getObjectMenu();
      // Add menu item to the Object menu
      // Locate position after "Convert To Actor..."
      int posConvertToActor = 0; 
      for(int i = 0; i < objectMenu.getChildCount(); i++){
        BMenuItem menuIter = null;
        menuIter = (BMenuItem) objectMenu.getChild(i);
        if(menuIter.getText().equalsIgnoreCase("Convert To Actor...")){
          posConvertToActor = i;
          break;
        }
      }
      posConvertToActor++; // select position after Convert To Actor or top if not found


      BMenuItem menuItem1 = Translate.menuItem("Convert to Cloth...", this, "convertToClothMenuAction");
      objectMenu.add(menuItem1, posConvertToActor);

      BMenuItem menuItem2 = Translate.menuItem("Generate Cloth Simulation...", this, "generateClothSimulationMenuAction");
      objectMenu.add(menuItem2, posConvertToActor+1);

      new ClothMenuItemActivator(layout, menuItem1, menuItem2).start();
    }

  }


  @SuppressWarnings("unused")
  private void convertToClothMenuAction(){
    Collection<ObjectInfo> sel = layout.getSelectedObjects();
    Object3D obj;
    Object3D mesh = null;
    ObjectInfo info;

    if (sel.size() != 1)
      return;
    info = (ObjectInfo) sel.toArray()[0];
    obj = info.getObject();
    if (obj.canConvertToTriangleMesh() == Object3D.CANT_CONVERT)
      return;

    ValueField errorField = new ValueField(DEFAULT_MESH_TOLERANCE, ValueField.POSITIVE);
    ValueField massDistField = new ValueField(DEFAULT_MASS_DISTANCE, ValueField.POSITIVE & ValueField.NONZERO);
    ComponentsDialog dlg = new ComponentsDialog(layout, Translate.text("selectToleranceForMesh"),
        new Widget [] {errorField,massDistField}, new String [] {Translate.text("maxError"),Translate.text("Mass Distance")});
    if (!dlg.clickedOk())
      return;
    meshTolerance = errorField.getValue();
    massDistance = massDistField.getValue();

    Cloth cloth = new Cloth(info, meshTolerance, massDistance, springConstant, dampingConstant, collisionDistance);
    ObjectInfo C = new ObjectInfo(cloth, new CoordinateSystem(), info.getName()+" (Cloth "+counter+")");
    C.coords.setOrigin(info.coords.getOrigin());
    C.coords.setOrientation(info.coords.getZDirection(), info.coords.getUpDirection());
    counter++;

    layout.setUndoRecord(new UndoRecord(layout, false, UndoRecord.COPY_OBJECT_INFO, new Object [] {info, info.duplicate()}));

    layout.addObject(C, new UndoRecord(layout, false));

    layout.removeObject(layout.getScene().indexOf(info), new UndoRecord(layout, false));
    layout.setSelection(layout.getScene().indexOf(C));
    Score theScore = layout.getScore();
    theScore.addTrack(layout.getSelectedObjects().toArray(), ClothTrack.class, new Object [] {layout.getScene()}, true);
    layout.updateImage();
    layout.updateMenus();

  }

  @SuppressWarnings("unused")
  private void generateClothSimulationMenuAction(){
    Collection<ObjectInfo> sel = layout.getSelectedObjects();
    Object3D mesh = null;
    ObjectInfo info;

    if (sel.size() != 1)
      return;
    info = (ObjectInfo) sel.toArray()[0];
    if (info.getObject() instanceof Cloth) {
      Cloth obj = (Cloth)info.getObject();
      if(obj.getEditor() == null) {
        obj.setEditor(new ClothSimEditorWindow((EditingWindow)layout, "Simulate Cloth", info)); 
      }

      obj.getEditor().setVisible(true);
    }

  }

}
