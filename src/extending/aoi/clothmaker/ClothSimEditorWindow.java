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

import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.GraphicsEnvironment;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Vector;

import artofillusion.MeshEditorWindow;
import artofillusion.MoveViewTool;
import artofillusion.RotateViewTool;
import artofillusion.UndoRecord;
import artofillusion.ViewerCanvas;
import artofillusion.math.Vec3;
import artofillusion.object.Mesh;
import artofillusion.object.Object3D;
import artofillusion.object.ObjectInfo;
import artofillusion.object.TriangleMesh;
import artofillusion.object.TriangleMesh.Edge;
import artofillusion.object.TriangleMesh.Face;
import artofillusion.object.TriangleMesh.Vertex;
import artofillusion.ui.EditingTool;
import artofillusion.ui.EditingWindow;
import artofillusion.ui.ToolPalette;
import artofillusion.ui.Translate;
import artofillusion.ui.UIUtilities;
import artofillusion.ui.ValueField;
import buoy.widget.*;

/**
 * Window for defining a Cloth Simulation
 * @author Timothy Fish
 *
 */
public class ClothSimEditorWindow extends MeshEditorWindow implements EditingWindow {

  static int subFrames = ClothMakerPlugin.DEFAULT_SUBFRAMES;
  private BMenuItem meshMenuItem[];
  boolean hideVert[], hideFace[], hideEdge[], selected[], showQuads, tolerant;
  boolean lockedPoints[];
  boolean dynamicPoints[];
  private int selectionDistance[], maxDistance, selectMode;
  private int drapeFrames = ClothMakerPlugin.DEFAULT_DRAPE_FRAMES;
  private int simFrames= ClothMakerPlugin.DEFAULT_SIM_FRAMES;
  private int gravityAxis = ClothMakerPlugin.DEFAULT_GRAVITY_AXIS;
  private double startTime = ClothMakerPlugin.DEFAULT_START_TIME;
  private double collision_distance = ClothMakerPlugin.DEFAULT_COLLISION_DISTANCE;
  private double gravity = ClothMakerPlugin.DEFAULT_GRAVITY;
  private double spring_constant = ClothMakerPlugin.DEFAULT_SPRING_CONST;
  private double damping_constant = ClothMakerPlugin.DEFAULT_DAMPING_CONST;
  private double vertex_mass = ClothMakerPlugin.DEFAULT_VERTEX_MASS;
  private ValueField timeField;
  private ValueField gravityField;
  private BComboBox axisChoice;
  private ValueField kField;
  private ValueField cField;
  private ValueField mField;
  private ValueField collisionField;
  private ValueField drapeFramesField;
  private ValueField simFramesField;
  private ValueField fpsField;
  private ValueField subFramesField;
  static double fps = ClothMakerPlugin.DEFAULT_FRAMES_PER_SECOND;
  private BCheckBox selfCheck;
  private BCheckBox floorCheck;
  private boolean selfCollision;
  private boolean floorCollision;
  protected static boolean lastProjectOntoSurface, lastTolerant, lastShowQuads;

  /**
   * Constructor
   * @param parent
   * @param title
   * @param obj
   */
  public ClothSimEditorWindow(EditingWindow parent, String title, ObjectInfo obj) {
    super(parent, title, obj);
    objInfo = obj;
    Cloth mesh = (Cloth) objInfo.getObject();
    hideVert = new boolean [mesh.getVertices().length];

    mesh.getSpringConst();

    FormContainer content = new FormContainer(new double [] {0, 1, 0}, new double [] {1, 0, 0});
    setContent(content);
    content.setDefaultLayout(new LayoutInfo(LayoutInfo.CENTER, LayoutInfo.BOTH, null, null));
    content.add(helpText = new BLabel(), 0, 1, 2, 1);
    content.add(viewsContainer, 1, 0);

    axisChoice = new BComboBox(new String [] {"X", "Y", "Z"});
    axisChoice.setSelectedIndex(gravityAxis);
    timeField = new ValueField(startTime, ValueField.NONE, 5);
    fpsField = new ValueField(fps, ValueField.NONZERO&ValueField.NONNEGATIVE & ValueField.INTEGER, 5);
    drapeFramesField = new ValueField(drapeFrames, ValueField.NONNEGATIVE & ValueField.INTEGER, 5);
    simFramesField = new ValueField(simFrames, ValueField.NONNEGATIVE & ValueField.INTEGER, 5);
    subFramesField = new ValueField(subFrames, ValueField.NONNEGATIVE & ValueField.NONZERO & ValueField.INTEGER);
    collisionField = new ValueField(collision_distance, ValueField.NONE, 5);
    gravityField = new ValueField(gravity, ValueField.NONE, 5);
    kField = new ValueField(spring_constant, ValueField.NONE, 5);
    cField = new ValueField(damping_constant, ValueField.NONE, 5);
    mField = new ValueField(vertex_mass, ValueField.NONZERO&ValueField.NONNEGATIVE, 5);
    selfCheck = new BCheckBox("Self Collision Detection", true);
    floorCheck = new BCheckBox("Floor Collision Detection", false);

    ColumnContainer parameters = new ColumnContainer();
    parameters.add(new BLabel("Gravity Axis"));
    parameters.add(axisChoice);
    parameters.add(new BLabel("Start Time"));
    parameters.add(timeField);
    parameters.add(new BLabel("Frames/Second"));
    parameters.add(fpsField);
    parameters.add(new BLabel("Drape Frames"));
    parameters.add(drapeFramesField);
    parameters.add(new BLabel("Sim Frames"));
    parameters.add(simFramesField);
    parameters.add(new BLabel("Sub Frames"));
    parameters.add(subFramesField);
    parameters.add(new BLabel("Collision Distance"));
    parameters.add(collisionField);
    parameters.add(new BLabel("Gravity"));
    parameters.add(gravityField);
    parameters.add(new BLabel("Spring Constant"));
    parameters.add(kField);
    parameters.add(new BLabel("Damping Constant"));
    parameters.add(cField);
    parameters.add(new BLabel("Vertex Mass"));
    parameters.add(mField);
    parameters.add(selfCheck);
    parameters.add(floorCheck);
    content.add(parameters, 2, 0);


    RowContainer buttons = new RowContainer();
    buttons.add(Translate.button("ok", this, "doOk"));
    buttons.add(Translate.button("cancel", this, "doCancel"));
    content.add(buttons, 0, 2, 2, 1, new LayoutInfo());

    FormContainer toolsContainer = new FormContainer(new double [] {1}, new double [] {1, 0});
    toolsContainer.setDefaultLayout(new LayoutInfo(LayoutInfo.NORTH, LayoutInfo.BOTH));
    content.add(toolsContainer, 0, 0);
    toolsContainer.add(tools = new ToolPalette(1, 3), 0, 0);
    EditingTool metaTool, altTool;
    tools.addTool(defaultTool = new SelectTool(this));
    tools.addTool(metaTool = new MoveViewTool(this));
    tools.addTool(altTool = new RotateViewTool(this));
    tools.setDefaultTool(defaultTool);
    tools.selectTool(defaultTool);


    for (int i = 0; i < theView.length; i++)
    {
      ClothMeshViewer view = (ClothMeshViewer) theView[i];
      view.setMetaTool(metaTool);
      view.setAltTool(altTool);
      view.setScene(parent.getScene(), obj);
      view.setFreehandSelection(lastFreehand);
    }

    createEditMenu(((Cloth) obj.getObject()).getTriangleMesh());
    createMeshMenu(((Cloth) obj.getObject()).getTriangleMesh());
    createViewMenu();
    recursivelyAddListeners(this);
    UIUtilities.applyDefaultFont(content);
    UIUtilities.applyDefaultBackground(content);
    Rectangle screenBounds = GraphicsEnvironment.getLocalGraphicsEnvironment().getMaximumWindowBounds();
    Dimension windowDim = new Dimension((screenBounds.width*3)/4, (screenBounds.height*3)/4);
    setBounds(new Rectangle((screenBounds.width-windowDim.width)/2, (screenBounds.height-windowDim.height)/2, windowDim.width, windowDim.height));
    tools.requestFocus();
    selected = new boolean [mesh.getVertices().length];
    lockedPoints = new boolean [mesh.getVertices().length];
    for (int i = 0; i < lockedPoints.length; i++) {
      lockedPoints[i] = mesh.getPinnedVertices()[i];
    }
    dynamicPoints = new boolean [mesh.getVertices().length];
    findQuads();
    findSelectionDistance();
    updateMenus();
  }

  /**
   * Add menu items to the Edit Menu
   * @param obj
   */
  void createEditMenu(TriangleMesh obj)
  {
    BMenu editMenu = Translate.menu("edit");
    menubar.add(editMenu);
    editMenu.add(undoItem = Translate.menuItem("undo", this, "undoCommand"));
    editMenu.add(redoItem = Translate.menuItem("redo", this, "redoCommand"));
    editMenu.addSeparator();

    editMenu.add(Translate.menuItem("selectAll", this, "selectAllCommand"));
    editMenu.add(Translate.menuItem("invertSelection", this, "invertSelectionCommand"));
  }

  /**
   * Add menu items to the Mesh Menu
   * @param obj
   */
  void createMeshMenu(TriangleMesh obj)
  {
    BMenu meshMenu = Translate.menu("mesh");
    menubar.add(meshMenu);
    meshMenuItem = new BMenuItem [2];

    meshMenu.add(meshMenuItem[0] = Translate.menuItem("Lock in Place", this, "lockPointsCommand"));
    meshMenu.add(meshMenuItem[1] = Translate.menuItem("Set as Dynamic", this, "dynamicPointsCommand"));

  }

  @Override
  public void setMesh(Mesh mesh) {}

  @Override
  public void updateMenus()
  {
    super.updateMenus();
    int i;

    for (i = 0; i < selected.length && !selected[i]; i++) {
      ;
    };
    if (i < selected.length)    {
      for (i = 0; i < meshMenuItem.length; i++)
        meshMenuItem[i].setEnabled(true);
    }
    else    {
      meshMenuItem[0].setEnabled(false);
      for (i = 2; i < meshMenuItem.length; i++)
        meshMenuItem[i].setEnabled(false);
    }
  }

  private void findQuads() {
    TriangleMesh mesh = ((Cloth) getObject().getObject()).getTriangleMesh();
    Vertex v[] = (Vertex []) mesh.getVertices();
    Edge e[] = mesh.getEdges();
    Face f[] = mesh.getFaces();
    if (hideEdge == null || hideEdge.length != e.length)
      hideEdge = new boolean [e.length];
    if (hideFace == null)
      for (int i = 0; i < e.length; i++)
        hideEdge[i] = false;
    else
      for (int i = 0; i < e.length; i++)
        hideEdge[i] = (hideFace[e[i].f1] && (e[i].f2 == -1 || hideFace[e[i].f2]));
    if (!showQuads)
      return;

    // An edge is a candidate for hiding if the two faces it borders are in the same plane.

    boolean candidate[] = new boolean [e.length];
    Vec3 norm[] = new Vec3 [f.length];
    for (int i = 0; i < f.length; i++)
    {
      Face fc = f[i];
      norm[i] = v[fc.v2].r.minus(v[fc.v1].r).cross(v[fc.v3].r.minus(v[fc.v1].r));
      double length = norm[i].length();
      if (length > 0.0)
        norm[i].scale(1.0/length);
    }
    for (int i = 0; i < e.length; i++)
      candidate[i] = (e[i].f2 != -1 && norm[e[i].f1].dot(norm[e[i].f2]) > 0.99);

    // Give every candidate edge a score for how close the adjoining faces are to forming
    // a rectangle.

    class EdgeScore implements Comparable<EdgeScore>      {
      public int edge;
      public double score;

      public EdgeScore(int edge, double score)
      {
        this.edge = edge;
        this.score = score;
      }

      @Override
      public int compareTo(EdgeScore o)
      {
        double diff = score-((EdgeScore) o).score;
        if (diff < 0.0)
          return -1;
        if (diff > 0.0)
          return 1;
        return 0;
      }
    }
    Vector<EdgeScore> scoreVec = new Vector<EdgeScore>(e.length);
    Vec3 temp0 = new Vec3(), temp1 = new Vec3(), temp2 = new Vec3();
    for (int i = 0; i < e.length; i++)
    {
      if (!candidate[i])
        continue;

      // Find the four vertices.

      Edge ed = e[i];
      int v1 = ed.v1, v2 = ed.v2, v3, v4;
      Face fc = f[ed.f1];
      if (fc.v1 != v1 && fc.v1 != v2)
        v3 = fc.v1;
      else if (fc.v2 != v1 && fc.v2 != v2)
        v3 = fc.v2;
      else
        v3 = fc.v3;
      fc = f[ed.f2];
      if (fc.v1 != v1 && fc.v1 != v2)
        v4 = fc.v1;
      else if (fc.v2 != v1 && fc.v2 != v2)
        v4 = fc.v2;
      else
        v4 = fc.v3;

      // Find the angles formed by them.

      temp0.set(v[v1].r.minus(v[v2].r));
      temp0.normalize();
      temp1.set(v[v1].r.minus(v[v3].r));
      temp1.normalize();
      temp2.set(v[v1].r.minus(v[v4].r));
      temp2.normalize();
      if (Math.acos(temp0.dot(temp1))+Math.acos(temp0.dot(temp2)) > Math.PI)
        continue;
      double dot = temp1.dot(temp2);
      double score = (dot > 0.0 ? dot : -dot);
      temp1.set(v[v2].r.minus(v[v3].r));
      temp1.normalize();
      temp2.set(v[v2].r.minus(v[v4].r));
      temp2.normalize();
      if (Math.acos(-temp0.dot(temp1))+Math.acos(-temp0.dot(temp2)) > Math.PI)
        continue;
      dot = temp1.dot(temp2);
      score += (dot > 0.0 ? dot : -dot);
      scoreVec.addElement(new EdgeScore(i, score));
    }
    if (scoreVec.isEmpty())
      return;

    // Sort them.

    EdgeScore score[] = new EdgeScore [scoreVec.size()];
    scoreVec.copyInto(score);
    Arrays.sort(score);

    // Mark which edges to hide.

    boolean hasHiddenEdge[] = new boolean [f.length];
    for (int i = 0; i < score.length; i++)
    {
      Edge ed = e[score[i].edge];
      if (hasHiddenEdge[ed.f1] || hasHiddenEdge[ed.f2])
        continue;
      hideEdge[score[i].edge] = true;
      hasHiddenEdge[ed.f1] = hasHiddenEdge[ed.f2] = true;
    }
  }

  @Override
  public void setSelectionMode(int mode) {}

  @Override
  public ObjectInfo getObject() {
    return objInfo;
  }

  @Override
  public int getSelectionMode() {
    return POINT_MODE;
  }

  /** Set the object being edited in this window. */

  public void setObject(Object3D obj)  {
    objInfo.setObject(obj);
    objInfo.clearCachedMeshes();
  }

  @Override
  public boolean[] getSelection() {
    return selected;
  }

  @Override
  public void setSelection(boolean[] sel) {
    selected = sel;
    updateMenus();
    for (ViewerCanvas view : theView)
      view.repaint();
  }

  @Override
  public int[] getSelectionDistance() {
    if (maxDistance != getTensionDistance())
      findSelectionDistance();
    return selectionDistance;
  }

  private void findSelectionDistance() {
    int i, j;
    TriangleMesh mesh = ((Cloth) getObject().getObject()).getTriangleMesh();
    int dist[] = new int [mesh.getVertices().length];
    Edge e[] = mesh.getEdges();

    maxDistance = getTensionDistance();

    // First, set each distance to 0 or -1, depending on whether that vertex is part of the
    // current selection.

    for (i = 0; i < dist.length; i++) {
      dist[i] = selected[i] ? 0 : -1;
    }

    // Now extend this outward up to maxDistance.
    for (i = 0; i < maxDistance; i++) {
      for (j = 0; j < e.length; j++)
      {
        if (hideEdge[j])
          continue;
        if (dist[e[j].v1] == -1 && dist[e[j].v2] == i)
          dist[e[j].v1] = i+1;
        else if (dist[e[j].v2] == -1 && dist[e[j].v1] == i)
          dist[e[j].v2] = i+1;
      }
    }
    selectionDistance = dist;
  }

  @Override
  public void deleteCommand() {}

  /** Select the entire mesh. */

  public void selectAllCommand()  {
    setUndoRecord(new UndoRecord(this, false, UndoRecord.SET_MESH_SELECTION, new Object [] {this, selectMode, selected.clone()}));
    for (int i = 0; i < selected.length; i++)
      selected[i] = true;
    setSelection(selected);
  }

  /** Invert the current selection. */

  public void invertSelectionCommand()  {
    boolean newSel[] = new boolean [selected.length];
    for (int i = 0; i < newSel.length; i++)
      newSel[i] = !selected[i];
    setUndoRecord(new UndoRecord(this, false, UndoRecord.SET_MESH_SELECTION, new Object [] {this, selectMode, selected}));
    setSelection(newSel);
  }

  public void lockPointsCommand() {
    // TODO setUndoRecord(new UndoRecord(this, false, UndoRecord.SET_MESH_SELECTION, new Object [] {this, selectMode, selected}));
    UndoRecord undo = new UndoRecord(this, true);
    undo.addCommand(UndoRecord.SET_MESH_SELECTION, new Object [] {this, lockedPoints, dynamicPoints});
    for (int i = 0; i < selected.length; i++) {
      if(selected[i]) {
        lockedPoints[i] = true;
        dynamicPoints[i] = false;
      }
    }

    Cloth mesh = (Cloth) objInfo.getObject();
    mesh.setPinnedVertices(lockedPoints);
  }

  public void dynamicPointsCommand() {
    // TODO setUndoRecord(new UndoRecord(this, false, UndoRecord.SET_MESH_SELECTION, new Object [] {this, selectMode, selected}));
    UndoRecord undo = new UndoRecord(this, true);
    undo.addCommand(UndoRecord.SET_MESH_SELECTION, new Object [] {this, lockedPoints, dynamicPoints});
    for (int i = 0; i < selected.length; i++) {
      if(selected[i]) {
        dynamicPoints[i] = true;
        lockedPoints[i] = false;
      }
    }

    Cloth mesh = (Cloth) objInfo.getObject();
    mesh.setPinnedVertices(lockedPoints);
  }

  @Override
  public void adjustDeltas(Vec3[] delta) {}


  @Override
  protected void doOk() {
    // TODO add save
    doSim();
    setVisible(false);
  }

  @Override
  protected void doCancel() {
    // TODO reset to pre-edit state
    setVisible(false);
  }

  /** 
   * Calculates the distorted simframes for the simulation.
   */
  protected void doSim() {
    ClothTrack theTrack = null;
    for(int i = 0; i < objInfo.getTracks().length; i++) {
      if(objInfo.getTracks()[i] instanceof ClothTrack) {
        theTrack = (ClothTrack) objInfo.getTracks()[i];
      }
    }

    if(theTrack == null) return;

    startTime = timeField.getValue();
    fps  = fpsField.getValue();
    drapeFrames = (int) drapeFramesField.getValue();
    simFrames = (int) simFramesField.getValue();
    subFrames = (int) subFramesField.getValue();
    gravity = gravityField.getValue();
    gravityAxis = axisChoice.getSelectedIndex();
    spring_constant = kField.getValue();
    damping_constant = cField.getValue();
    vertex_mass = mField.getValue();
    collision_distance = collisionField.getValue();
    selfCollision = selfCheck.getState();
    floorCollision = floorCheck.getState();

    theTrack.setParams(startTime, fps*subFrames, gravity, gravityAxis, spring_constant, damping_constant, vertex_mass, collision_distance, selfCollision, floorCollision, drapeFrames, simFrames);

    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    Object3D original = getObject().getObject();
    String orgTitle = new String(this.getTitle());
    String workingTitle = new String(orgTitle)+"[ ";
    int maxFrames = simFrames*subFrames*subFrames;
    for(int i = -drapeFrames*subFrames; i <= maxFrames; i++) {
      SimFrame frame = theTrack.simulateCloth(i);
      
      
      this.setTitle(workingTitle+i+" : "+maxFrames+" ]");

      this.updateImage();

      setObject(frame.M);

    }
    this.setTitle(orgTitle);
    setObject(original);
    setCursor(Cursor.getDefaultCursor());

  }

  /**
   * Get reference to the lockedPoints
   * @return
   */
  public boolean[] getLockedPoints() {
    return lockedPoints;
  }

  /**
   * Get reference to the dynamicPoints
   * @return
   */
  public boolean[] getDynamicPoints() {
    return dynamicPoints;
  }

  /**
   * Get Frames Per Second
   * @return
   */
  public static double getFPS() {
    return fps;
  }

  @Override
  public void setVisible(boolean visible)
  {

    ClothTrack theTrack = null;
    for(int i = 0; i < objInfo.getTracks().length; i++) {
      if(objInfo.getTracks()[i] instanceof ClothTrack) {
        theTrack = (ClothTrack) objInfo.getTracks()[i];
      }
    }

    if(theTrack != null) {
      this.collision_distance = theTrack.collision_distance;

      timeField.setValue(theTrack.startTime);
      fpsField.setValue(theTrack.tfps/ClothSimEditorWindow.subFrames);
      drapeFramesField.setValue(theTrack.drapeFrames);
      simFramesField.setValue(theTrack.simFrames);
      gravityField.setValue(theTrack.gravity);
      axisChoice.setSelectedIndex(theTrack.gravityAxis);
      kField.setValue(theTrack.spring_constant);
      cField.setValue(theTrack.damping_constant);
      mField.setValue(theTrack.vertex_mass);
      collisionField.setValue(theTrack.collision_distance);
      selfCheck.setState(theTrack.selfCollision);
      floorCheck.setState(theTrack.floorCollision);      
    }

    super.setVisible(visible);
  }


}
