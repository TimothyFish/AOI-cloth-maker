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

import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;

import artofillusion.MeshViewer;
import artofillusion.RenderingMesh;
import artofillusion.TextureParameter;
import artofillusion.UndoRecord;
import artofillusion.math.RGBColor;
import artofillusion.math.Vec2;
import artofillusion.math.Vec3;
import artofillusion.object.MeshVertex;
import artofillusion.object.ObjectInfo;
import artofillusion.object.TriangleMesh;
import artofillusion.object.TriangleMesh.Edge;
import artofillusion.object.TriangleMesh.Face;
import artofillusion.object.TriangleMesh.Vertex;
import artofillusion.texture.FaceParameterValue;
import artofillusion.ui.EditingTool;
import artofillusion.ui.MeshEditController;
import artofillusion.view.ConstantVertexShader;
import artofillusion.view.FlatVertexShader;
import artofillusion.view.ParameterVertexShader;
import artofillusion.view.SelectionVertexShader;
import artofillusion.view.SmoothVertexShader;
import artofillusion.view.TexturedVertexShader;
import artofillusion.view.VertexShader;
import buoy.event.WidgetMouseEvent;
import buoy.widget.RowContainer;

/**
 * Viewer for displaying cloth objects.
 * @author Timothy Fish
 *
 */
public class ClothMeshViewer extends MeshViewer {
  private boolean draggingSelectionBox, dragging;
  private int deselect;
  private Point screenVert[];
  private double screenZ[];
  private Vec2 screenVec2[];
  boolean visible[];

  /**
   * Constructor
   * @param controller
   * @param p
   */
  public ClothMeshViewer(MeshEditController window, RowContainer p) {
    super(window, p);
    TriangleMesh mesh = ((Cloth) window.getObject().getObject()).getTriangleMesh();
    visible = new boolean [mesh.getVertices().length];
  }

  @Override
  public void updateImage()
  {
    TriangleMesh mesh = ((Cloth) getController().getObject().getObject()).getTriangleMesh();
    MeshVertex v[] = mesh.getVertices();

    // Calculate the screen coordinates of every vertex.

    screenVert = new Point [v.length];
    screenZ = new double [v.length];
    if (visible.length != v.length)
      visible = new boolean [v.length];
    screenVec2 = new Vec2 [v.length];
    double clipDist = theCamera.getClipDistance();
    boolean hideVert[] = (controller instanceof ClothSimEditorWindow ? ((ClothSimEditorWindow) controller).hideVert : new boolean [v.length]);
    for (int i = 0; i < v.length; i++)
    {
      Vec3 pos = v[i].r;
      screenVec2[i] = theCamera.getObjectToScreen().timesXY(pos);
      screenVert[i] = new Point((int) screenVec2[i].x, (int) screenVec2[i].y);
      screenZ[i] = theCamera.getObjectToView().timesZ(pos);
      visible[i] = (!hideVert[i] && screenZ[i] > clipDist);
    }
    super.updateImage();
  }

  @Override
  protected void drawObject() {
    // Now draw the object.
    drawSurface();
    Color meshColor = lineColor;
    Color selectedColor = highlightColor;
    Color frozenColor = Color.blue;
    Color dynamicColor = meshColor;

    drawEdges(screenVec2, disabledColor, disabledColor);
    drawVertices(meshColor, currentTool.hilightSelection() ? selectedColor : meshColor,
        dynamicColor, frozenColor);
  }



  private void drawSurface() {
    if (!showSurface)
      return;
    boolean hide[] = null;
    int faceIndex[] = null;
    ObjectInfo objInfo = controller.getObject();
    if (controller instanceof ClothSimEditorWindow && ((ClothSimEditorWindow) controller).getFaceIndexParameter() != null)
    {
      RenderingMesh mesh = objInfo.getPreviewMesh();
      TextureParameter faceIndexParameter = ((ClothSimEditorWindow) controller).getFaceIndexParameter();
      double param[] = null;
      for (int i = 0; i < mesh.param.length; i++)
        if (objInfo.getObject().getParameters()[i] == faceIndexParameter)
          param = ((FaceParameterValue) mesh.param[i]).getValue();
      faceIndex = new int [param.length];
      for (int i = 0; i < faceIndex.length; i++)
        faceIndex[i] = (int) param[i];
      boolean hideFace[] = ((ClothSimEditorWindow) controller).hideFace;
      if (hideFace != null)
      {
        hide = new boolean [param.length];
        for (int i = 0; i < hide.length; i++)
          hide[i] = hideFace[faceIndex[i]];
      }
    }
    Vec3 viewDir = getDisplayCoordinates().toLocal().timesDirection(theCamera.getViewToWorld().timesDirection(Vec3.vz()));
    if (renderMode == RENDER_WIREFRAME)
      renderWireframe(objInfo.getWireframePreview(), theCamera, surfaceColor);
    else if (renderMode == RENDER_TRANSPARENT)
    {
      VertexShader shader = new ConstantVertexShader(transparentColor);
      if (faceIndex != null && controller.getSelectionMode() == MeshEditController.FACE_MODE)
        shader = new SelectionVertexShader(new RGBColor(1.0, 0.4, 1.0), shader, faceIndex, controller.getSelection());
      renderMeshTransparent(objInfo.getPreviewMesh(), shader, theCamera, viewDir, hide);
    }
    else
    {
      RenderingMesh mesh = objInfo.getPreviewMesh();
      VertexShader shader;
      if (renderMode == RENDER_FLAT)
        shader = new FlatVertexShader(mesh, surfaceRGBColor, viewDir);
      else if (surfaceColoringParameter != null)
      {
        shader = null;
        TextureParameter params[] = objInfo.getObject().getParameters();
        for (int i = 0; i < params.length; i++)
          if (params[i].equals(surfaceColoringParameter))
          {
            shader = new ParameterVertexShader(mesh, mesh.param[i], lowValueColor, highValueColor, surfaceColoringParameter.minVal, surfaceColoringParameter.maxVal, viewDir);
            break;
          }
      }
      else if (renderMode == RENDER_SMOOTH)
        shader = new SmoothVertexShader(mesh, surfaceRGBColor, viewDir);
      else
        shader = new TexturedVertexShader(mesh, objInfo.getObject(), 0.0, viewDir).optimize();
      if (faceIndex != null && controller.getSelectionMode() == MeshEditController.FACE_MODE)
        shader = new SelectionVertexShader(new RGBColor(1.0, 0.4, 1.0), shader, faceIndex, controller.getSelection());
      renderMesh(mesh, shader, theCamera, objInfo.getObject().isClosed(), hide);
    }

  }

  /**
   * Draws the vertices in color so the user can tell if the vertices are selected or if they have
   * been set to locked or dynamic.
   * @param unselectedColor
   * @param selectedColor
   * @param dynamicColor
   * @param frozenColor
   */
  private void drawVertices(Color unselectedColor, Color selectedColor, Color dynamicColor, Color frozenColor) {
    if (!showMesh)
      return;

    MeshVertex v[] = ((Cloth) getController().getObject().getObject()).getTriangleMesh().getVertices();

    ArrayList<Rectangle> boxes = new ArrayList<Rectangle>();
    ArrayList<Double> depths = new ArrayList<Double>();
    boolean selected[] = controller.getSelection();
    boolean frozen[] = ((ClothSimEditorWindow)controller).getLockedPoints();
    boolean dynamic[] = ((ClothSimEditorWindow)controller).getDynamicPoints();

    // First, draw any unselected and unspecified portions of the object.
    for (int i = 0; i < v.length; i++) {
      if (!selected[i] && !frozen[i] && !dynamic[i] && visible[i]) {
        boxes.add(new Rectangle(screenVert[i].x-HANDLE_SIZE/2, screenVert[i].y-HANDLE_SIZE/2, HANDLE_SIZE, HANDLE_SIZE));
        depths.add(screenZ[i]-0.02);
      }
    }
    renderBoxes(boxes, depths, unselectedColor);

    // draw the portions that are fixed in place
    boxes.clear();
    depths.clear();
    for (int i = 0; i < v.length; i++) {
      if (!selected[i] && frozen[i] && visible[i]) {
        boxes.add(new Rectangle(screenVert[i].x-HANDLE_SIZE/2, screenVert[i].y-HANDLE_SIZE/2, HANDLE_SIZE, HANDLE_SIZE));
        depths.add(screenZ[i]-0.02);
      }
    }
    renderBoxes(boxes, depths, frozenColor);

    // draw the dynamic portions
    boxes.clear();
    depths.clear();
    for (int i = 0; i < v.length; i++) {
      if (!selected[i] && dynamic[i] && visible[i]) {
        boxes.add(new Rectangle(screenVert[i].x-HANDLE_SIZE/2, screenVert[i].y-HANDLE_SIZE/2, HANDLE_SIZE, HANDLE_SIZE));
        depths.add(screenZ[i]-0.02);
      }
    }
    renderBoxes(boxes, depths, dynamicColor);

    // Now draw the selected portions.
    boxes.clear();
    depths.clear();
    for (int i = 0; i < v.length; i++) {
      if (selected[i] && visible[i]) {
        boxes.add(new Rectangle(screenVert[i].x-HANDLE_SIZE/2, screenVert[i].y-HANDLE_SIZE/2, HANDLE_SIZE, HANDLE_SIZE));
        depths.add(screenZ[i]-0.02);
      }
    }
    renderBoxes(boxes, depths, selectedColor);

  }

  /**
   * Draws lines between the vertices.
   * @param p
   * @param unselectedColor
   * @param selectedColor
   */
  private void drawEdges(Vec2 p[], Color unselectedColor, Color selectedColor) {
    if (!showMesh)
      return;
    Edge e[] = ((Cloth) getController().getObject().getObject()).getTriangleMesh().getEdges();

    // Determine which edges are selected.

    int selectMode = controller.getSelectionMode();
    boolean hideEdge[] = (controller instanceof ClothSimEditorWindow ? ((ClothSimEditorWindow) controller).hideEdge : new boolean [e.length]);
    boolean selected[] = controller.getSelection();
    boolean isSelected[];
    if (selectMode == MeshEditController.POINT_MODE)
      isSelected = new boolean [e.length];
    else if (selectMode == MeshEditController.EDGE_MODE)
      isSelected = selected;
    else
    {
      isSelected = new boolean [e.length];
      for (int i = 0; i < e.length; i++) {
        isSelected[i] = (selected[e[i].f1] || (e[i].f2 > -1 && selected[e[i].f2]));
      }
    }

    // Draw the edges of the control mesh.

    for (int j = 0; j < 2; j++)
    {
      // Draw unselected edges on the first pass, selected edges on the second pass.

      boolean showSelected = (j == 1);
      Color color = (showSelected ? selectedColor : unselectedColor);
      for (int i = 0; i < e.length; i++)
        if (showSelected == isSelected[i] && !hideEdge[i] && visible[e[i].v1] && visible[e[i].v2])
          renderLine(p[e[i].v1], screenZ[e[i].v1]-0.01, p[e[i].v2], screenZ[e[i].v2]-0.01, theCamera, color);
    }    
  }

  /** When the user presses the mouse, forward events to the current tool as appropriate.
      If this is a vertex based tool, allow them to select or deselect vertices. */

  @Override
  protected void mousePressed(WidgetMouseEvent e)
  {
    TriangleMesh mesh = ((Cloth) getController().getObject().getObject()).getTriangleMesh();
    Edge ed[] = mesh.getEdges();
    Face f[] = mesh.getFaces();
    int i, j, k;

    requestFocus();
    sentClick = false;
    deselect = -1;
    dragging = false;
    clickPoint = e.getPoint();

    // Determine which tool is active.

    if (metaTool != null && e.isMetaDown())
      activeTool = metaTool;
    else if (altTool != null && e.isAltDown())
      activeTool = altTool;
    else
      activeTool = currentTool;

    // If the current tool wants all clicks, just forward the event.

    if ((activeTool.whichClicks() & EditingTool.ALL_CLICKS) != 0)
    {
      activeTool.mousePressed(e, this);
      dragging = true;
      sentClick = true;
    }
    boolean allowSelectionChange = activeTool.allowSelectionChanges();
    boolean wantHandleClicks = ((activeTool.whichClicks() & EditingTool.HANDLE_CLICKS) != 0);
    if (!allowSelectionChange && !wantHandleClicks)
      return;

    // Determine what the click was on.

    i = findClickTarget(e.getPoint(), null);

    // If the click was not on an object, start dragging a selection box.

    if (i == -1)
    {
      if (allowSelectionChange)
      {
        draggingSelectionBox = true;
        beginDraggingSelection(e.getPoint(), false);
      }
      return;
    }

    // If we are in edge or face selection mode, find a vertex of the clicked edge or face,
    // so that it can be passed to editing tools.

    if (controller.getSelectionMode() == MeshEditController.EDGE_MODE)
    {
      if (visible[ed[i].v1])
        j = ed[i].v1;
      else
        j = ed[i].v2;
    }
    else if (controller.getSelectionMode() == MeshEditController.FACE_MODE)
    {
      if (visible[f[i].v1])
        j = f[i].v1;
      else if (visible[f[i].v2])
        j = f[i].v2;
      else
        j = f[i].v3;
    }
    else
      j = i;

    // If the click was on a selected object, forward it to the current tool.  If it was a
    // shift-click, the user may want to deselect it, so set a flag.

    boolean selected[] = controller.getSelection();
    if (selected[i])
    {
      if (e.isShiftDown() && allowSelectionChange)
        deselect = i;
      if (wantHandleClicks)
        activeTool.mousePressedOnHandle(e, this, 0, j);
      sentClick = true;
      return;
    }
    if (!allowSelectionChange)
      return;

    // The click was on an unselected object.  Select it and send an event to the current tool.

    boolean oldSelection[] = selected.clone();
    if (!e.isShiftDown())
      for (k = 0; k < selected.length; k++)
        selected[k] = false;
    selected[i] = true;

    currentTool.getWindow().setUndoRecord(new UndoRecord(currentTool.getWindow(), false, UndoRecord.SET_MESH_SELECTION, new Object [] {controller, controller.getSelectionMode(), oldSelection}));
    controller.setSelection(selected);
    currentTool.getWindow().updateMenus();
    if (!e.isShiftDown() && wantHandleClicks)
    {
      activeTool.mousePressedOnHandle(e, this, 0, j);
      sentClick = true;
    }
  }

  @Override
  protected void mouseDragged(WidgetMouseEvent e)
  {
    if (!dragging)
    {
      Point p = e.getPoint();
      if (Math.abs(p.x-clickPoint.x) < 2 && Math.abs(p.y-clickPoint.y) < 2)
        return;
    }
    dragging = true;
    deselect = -1;
    super.mouseDragged(e);
  }

  /**
   * Positions based on snap to grid value.
   * @param e
   */
  void positionToGrid(WidgetMouseEvent e)
  {
    Point pos = e.getPoint();
    Vec3 v;
    Vec2 v2;

    if (!snapToGrid || isPerspective())
      return;
    v = theCamera.convertScreenToWorld(pos, theCamera.getDistToScreen());
    v2 = theCamera.getWorldToScreen().timesXY(v);
    e.translatePoint((int) v2.x - pos.x, (int) v2.y - pos.y);
  }

  @Override
  protected void mouseReleased(WidgetMouseEvent e)
  {
    TriangleMesh mesh = ((Cloth) getController().getObject().getObject()).getTriangleMesh();
    positionToGrid(e);
    endDraggingSelection();
    boolean selected[] = controller.getSelection();
    boolean oldSelection[] = selected.clone();
    if (draggingSelectionBox && !e.isShiftDown() && !e.isControlDown())
      for (int i = 0; i < selected.length; i++)
        selected[i] = false;

    // If the user was dragging a selection box, then select or deselect anything
    // it intersects.

    boolean hideVert[] = new boolean [mesh.getVertices().length];
    if (selectBounds != null)
    {
      boolean newsel = !e.isControlDown();
      for (int i = 0; i < selected.length; i++)
        if (!hideVert[i] && selectionRegionContains(screenVert[i]))
          selected[i] = newsel;
    }
    draggingBox = draggingSelectionBox = false;

    // Send the event to the current tool, if appropriate.

    if (sentClick)
    {
      if (!dragging)
      {
        Point p = e.getPoint();
        e.translatePoint(clickPoint.x-p.x, clickPoint.y-p.y);
      }
      activeTool.mouseReleased(e, this);
    }

    // If the user shift-clicked a selected point and released the mouse without dragging,
    // then deselect the point.

    if (deselect > -1)
    {
      selected[deselect] = false;
    }
    for (int k = 0; k < selected.length; k++)
      if (selected[k] != oldSelection[k])
      {
        currentTool.getWindow().setUndoRecord(new UndoRecord(currentTool.getWindow(), false, UndoRecord.SET_MESH_SELECTION, new Object [] {controller, controller.getSelectionMode(), oldSelection}));
        break;
      }
    controller.setSelection(selected);
    currentTool.getWindow().updateMenus();
  }

  /** Determine which vertex, edge, or face (depending on the current selection mode) the
      mouse was clicked on.  If the click was on top of multiple objects, priority is given
      to ones which are currently selected, and then to ones which are in front.  If the
      click is not over any object, -1 is returned. */

  public int findClickTarget(Point pos, Vec3 uvw)
  {
    double z, closestz = Double.MAX_VALUE;
    boolean sel = false;
    int which = -1;

    boolean selected[] = controller.getSelection();
    boolean priorityToSelected = (getRenderMode() == RENDER_WIREFRAME || getRenderMode() == RENDER_TRANSPARENT);

    TriangleMesh mesh = ((Cloth) getController().getObject().getObject()).getTriangleMesh();
    Vertex vt[] = (Vertex []) mesh.getVertices();
    for (int i = 0; i < vt.length; i++)
    {
      if (!visible[i])
        continue;
      if (sel && !selected[i] && priorityToSelected)
        continue;
      Point v1 = screenVert[i];
      if (pos.x < v1.x-HANDLE_SIZE/2 || pos.x > v1.x+HANDLE_SIZE/2 ||
          pos.y < v1.y-HANDLE_SIZE/2 || pos.y > v1.y+HANDLE_SIZE/2)
        continue;
      z = theCamera.getObjectToView().timesZ(vt[i].r);
      if (z < closestz || (!sel && selected[i] && priorityToSelected))
      {
        which = i;
        closestz = z;
        sel = selected[i];
      }
    }

    return which;
  }
}
