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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InvalidObjectException;
import java.lang.ref.SoftReference;
import java.util.Vector;

import artofillusion.MeshViewer;
import artofillusion.RenderingMesh;
import artofillusion.RenderingTriangle;
import artofillusion.Scene;
import artofillusion.TextureParameter;
import artofillusion.WireframeMesh;
import artofillusion.animation.Keyframe;
import artofillusion.animation.Skeleton;
import artofillusion.material.Material;
import artofillusion.material.MaterialMapping;
import artofillusion.math.BoundingBox;
import artofillusion.math.Vec3;
import artofillusion.object.Mesh;
import artofillusion.object.MeshVertex;
import artofillusion.object.Object3D;
import artofillusion.object.ObjectInfo;
import artofillusion.object.TriangleMesh;
import artofillusion.object.TriangleMesh.Edge;
import artofillusion.object.TriangleMesh.Face;
import artofillusion.object.TriangleMesh.Vertex;
import artofillusion.texture.ParameterValue;
import artofillusion.texture.Texture;
import artofillusion.texture.TextureMapping;
import artofillusion.ui.EditingWindow;
import artofillusion.ui.MeshEditController;
import buoy.widget.RowContainer;


/**
 * An scene object with vertices that can be manipulated by a cloth simulator
 * to give the appearance of cloth in the scene. Any object that can be
 * converted to a TriangleMesh can be converted to a Cloth. Springs are
 * created between vertices with resting lengths equal to the distance
 * between vertices. 
 * 
 * @author Timothy Fish
 *
 */
public class Cloth extends Object3D implements Mesh {
  private TriangleMesh theTriangleMesh; // triangle mesh used internally for consistency
  private Mass theMasses[]; // one mass per vector of theTriangleMesh
  private Vector<Spring> theSprings; // each Spring connects two masses, but not all masses have springs
  private double springConstant; // tensile spring constant
  private double dampingConstant; // damping factor
  private double collisionDistance; // how far from cloth collisions should be detected
  private double meshTolerance; // defines the fine detail of theTriangleMesh
  private double vertMass = 0.5; // per vertex weight of cloth
  private SoftReference<RenderingMesh> cachedMesh; // stored rendering mesh
  private boolean[] pinnedVerts; // specifies which vertices are locked in place during simulation
  private ClothSimEditorWindow editor; // reference to window that edits the cloth simulation
  private Vector<SimFrame> storedFrames; // reference to frames stored from simulation
  private ObjectInfo theObject; // reference to the object that was converted into a cloth


  /**
   * constructor
   * Converts obj into a cloth.    
   * @param obj
   * @param tol
   * @param massDist
   * @param springConst
   * @param dampConst
   * @param collDist
   */
  public Cloth(ObjectInfo obj, double tol, double massDist, double springConst, double dampConst, double collDist) {
    theObject = obj;
    meshTolerance = tol;
    springConstant = springConst;
    dampingConstant = dampConst;
    collisionDistance = collDist;

    theTriangleMesh = theObject.getObject().convertToTriangleMesh(tol);
    theTriangleMesh = theTriangleMesh.subdivideToLimit(massDist);
    setTexture(theObject.getObject().getTexture(), theObject.getObject().getTextureMapping());
    setMaterial(theObject.getObject().getMaterial(), theObject.getObject().getMaterialMapping());

    theMasses = new Mass[theTriangleMesh.getVertexPositions().length];
    pinnedVerts = new boolean[theTriangleMesh.getVertexPositions().length];
    for(int i = 0; i < theMasses.length; i++) {
      theMasses[i] = new Mass(theTriangleMesh.getVertexPositions()[i], i, vertMass, new Vec3());
      pinnedVerts[i] = false;
    }

    theSprings = new Vector<Spring>();

    // Add first layer springs
    Edge[] meshEdges = theTriangleMesh.getEdges();
    for(int i = 0; i < meshEdges.length; i++) {
      Mass massA = theMasses[meshEdges[i].v1];
      Mass massB = theMasses[meshEdges[i].v2];
      Spring s = new Spring(massA, massB, massA.getPosition().distance(massB.getPosition()), springConstant);
      if (!springExists(s)) {
        theSprings.add(s);
        massA.connectToSpring(s);
        massB.connectToSpring(s);
      }
    }
  }

  /**
   * Returns true if the spring already exists in the list of springs.
   * @param s
   * @return
   */
  private boolean springExists(Spring s) {
    for(Spring existingS : theSprings) {
      if(existingS.equals(s)) {
        return true;
      }
    }

    storedFrames = null;

    return false;
  }

  /**
   * Copy Constructor
   * @param cloth
   */
  public Cloth(Cloth cloth) {
    meshTolerance = cloth.meshTolerance;
    theTriangleMesh = (TriangleMesh) cloth.theTriangleMesh.duplicate();
    theMasses = new Mass[cloth.theMasses.length];
    pinnedVerts = new boolean[cloth.theMasses.length];
    for(int i = 0; i < theMasses.length; i++) {
      theMasses[i] = cloth.theMasses[i];
      pinnedVerts[i] = cloth.pinnedVerts[i];
    }
    theSprings = new Vector<Spring>();
    for(Spring S : cloth.theSprings) {
      theSprings.add(S);
    }
    springConstant = cloth.springConstant;
    dampingConstant = cloth.dampingConstant;
    collisionDistance = cloth.collisionDistance;

    storedFrames = cloth.storedFrames;

  }

  /**
   * Returns a reference to the list of masses.
   * @return
   */
  public Mass[] getMasses() {
    return theMasses;
  }

  /**
   * Returns an array with all of the springs in it.
   * @return
   */
  public Spring[] getSprings() {
    Spring[] ret = new Spring[theSprings.size()];
    for(int i = 0; i < theSprings.size(); i++) {
      ret[i] = theSprings.elementAt(i);
    }
    return ret;
  }

  /**
   * Get the spring constant.
   * @return
   */
  public double getSpringConst() { return springConstant; }

  /**
   * Get the damping constant.
   * @return
   */
  public double getDampingConst() { return dampingConstant; }

  /**
   * Get the collision distance.
   * @return
   */
  public double getCollisionDist() { return collisionDistance; }

  @Override
  public MeshVertex[] getVertices() {
    if(theTriangleMesh == null) return null;

    return theTriangleMesh.getVertices();  
  }

  @Override
  public Vec3[] getVertexPositions() {
    if(theTriangleMesh == null) return null;

    return theTriangleMesh.getVertexPositions();
  }

  @Override
  public void setVertexPositions(Vec3[] v) {
    if(theTriangleMesh == null) return;

    theTriangleMesh.setVertexPositions(v);
    for(int i = 0; i < v.length; i++) {
      theMasses[i].setPosition(v[i]);
    }
  }

  @Override
  public BoundingBox getBounds() {
    if(theTriangleMesh == null) return null;

    return theTriangleMesh.getBounds();
  }

  @Override
  public Vec3[] getNormals() {
    if(theTriangleMesh == null) return null;

    return theTriangleMesh.getNormals();
  }

  @Override
  public TextureParameter[] getParameters() {
    return super.getParameters();
  }

  @Override
  public ParameterValue[] getParameterValues() {
    return super.getParameterValues();
  }


  @Override
  public Object3D duplicate() {
    Object3D obj = (Object3D)new Cloth(this);
    obj.copyTextureAndMaterial(this);
    return obj;
  }

  @Override
  public void copyObject(Object3D obj) {
    Cloth cloth = (Cloth) obj;
    meshTolerance = cloth.meshTolerance;
    theTriangleMesh = (TriangleMesh) cloth.theTriangleMesh.duplicate();
    theTriangleMesh.copyTextureAndMaterial(obj);
    theMasses = new Mass[cloth.theMasses.length];
    pinnedVerts = new boolean[cloth.theMasses.length];
    for(int i = 0; i < theMasses.length; i++) {
      theMasses[i] = cloth.theMasses[i];
      pinnedVerts[i] = cloth.pinnedVerts[i];
    }
    theSprings.clear();
    for(Spring spring : cloth.theSprings) {
      theSprings.add(spring);
    }
    springConstant = cloth.springConstant;
    dampingConstant = cloth.dampingConstant;
    collisionDistance = cloth.collisionDistance;

    copyTextureAndMaterial(obj);
    cachedMesh = null;

    storedFrames = cloth.storedFrames;
  }

  @Override
  public Skeleton getSkeleton() {
    if(theTriangleMesh == null) return null;

    return theTriangleMesh.getSkeleton();
  }

  @Override
  public void setSkeleton(Skeleton s) {
    if(theTriangleMesh != null) {
      theTriangleMesh.setSkeleton(s);
    }

  }

  @Override
  public MeshViewer createMeshViewer(MeshEditController controller, RowContainer options) {
    return new ClothMeshViewer(controller, options);
  }

  @Override
  public void setSize(double xsize, double ysize, double zsize) {
    if(theTriangleMesh != null) {    
      theTriangleMesh.setSize(xsize,  ysize, zsize);
    }
  }

  @Override
  public WireframeMesh getWireframeMesh() {
    if(theTriangleMesh == null) return null;

    return theTriangleMesh.getWireframeMesh();
  }

  @Override
  public Keyframe getPoseKeyframe() {
    if(theTriangleMesh == null) return null;

    return theTriangleMesh.getPoseKeyframe();
  }

  @Override
  public void applyPoseKeyframe(Keyframe k) {
    if(theTriangleMesh != null) {
      theTriangleMesh.applyPoseKeyframe(k);
    }
  }

  @Override
  public boolean isClosed()
  {
    return false;
  }

  @Override
  public RenderingMesh getRenderingMesh(double tol, boolean interactive, ObjectInfo info)
  {
    TriangleMesh mesh = this.theTriangleMesh;
    Vec3 vert[], normalArray[];
    Vector<Vec3> norm;
    Vertex v[];
    Edge e[];
    Face f[], tempFace;
    RenderingTriangle tri[];
    int i, j, k, m, first, last, normals, ed[], facenorm[];
    boolean split[];


    if (interactive && cachedMesh != null)
    {
      RenderingMesh cached = cachedMesh.get();
      if (cached != null)
        return cached;
    }
    if (mesh.getFaces().length == 0)
    {
      RenderingMesh rend = new RenderingMesh(new Vec3 [] {new Vec3()}, new Vec3 [] {Vec3.vx()}, new RenderingTriangle [0], texMapping, matMapping);
      rend.setParameters(mesh.getParameterValues());
      return rend;
    }

    // If appropriate, subdivide the mesh.

    if (mesh.getSmoothingMethod() == INTERPOLATING || mesh.getSmoothingMethod() == APPROXIMATING)
    {
      double tol2 = tol*tol;
      Vec3 diff = new Vec3();
      split = new boolean [mesh.getEdges().length];
      for (i = 0; i < split.length; i++)
      {
        Vec3 r1 = mesh.getVertex(mesh.getEdges()[i].v1).r, r2 = mesh.getVertex(mesh.getEdges()[i].v2).r;
        diff.set(r1.x-r2.x, r1.y-r2.y, r1.z-r2.z);
        split[i] = (diff.length2() > tol2);
      }
      if (mesh.getSmoothingMethod() == INTERPOLATING)
        mesh = TriangleMesh.subdivideButterfly(mesh, split, tol);
      else
        mesh = TriangleMesh.subdivideLoop(mesh, split, tol);
      v = (Vertex[]) mesh.getVertices();
      e = mesh.getEdges();
      f = mesh.getFaces();
    }
    else
    {
      v = (Vertex[]) mesh.getVertices();
      e = mesh.getEdges();
      f = mesh.getFaces();
    }

    // Create the RenderingMesh.

    vert = new Vec3 [v.length];
    norm = new Vector<Vec3>();
    tri = new RenderingTriangle [f.length];
    facenorm = new int [f.length*3];
    normals = 0;
    k = last = 0;
    if (mesh.getSmoothingMethod() != NO_SMOOTHING)
    {
      // The mesh needs to be smooth shaded, so we need to calculate the normal vectors.
      // There may be more than one normal associated with a vertex, if that vertex is
      // on a crease.  Begin by finding a "true" normal for each face.

      Vec3 trueNorm[] = new Vec3 [f.length];
      for (i = 0; i < f.length; i++)
      {
        trueNorm[i] = v[f[i].v2].r.minus(v[f[i].v1].r).cross(v[f[i].v3].r.minus(v[f[i].v1].r));
        double length = trueNorm[i].length();
        if (length > 0.0)
          trueNorm[i].scale(1.0/length);
      }

      // Now loop over each vertex.

      for (i = 0; i < v.length; i++)
      {
        vert[i] = v[i].r;
        ed = v[i].getEdges();

        // If this vertex is a corner, we can just set its normal to null.

        if (v[i].smoothness < 1.0f)
        {
          norm.addElement(null);
          for (j = 0; j < ed.length; j++)
          {
            k = e[ed[j]].f1;
            tempFace = f[k];
            if (tempFace.v1 == i)
              facenorm[k*3] = normals;
            else if (tempFace.v2 == i)
              facenorm[k*3+1] = normals;
            else
              facenorm[k*3+2] = normals;
            k = e[ed[j]].f2;
            if (k != -1)
            {
              tempFace = f[k];
              if (tempFace.v1 == i)
                facenorm[k*3] = normals;
              else if (tempFace.v2 == i)
                facenorm[k*3+1] = normals;
              else
                facenorm[k*3+2] = normals;
            }
          }
          normals++;
          continue;
        }

        // If any of the edges intersecting this vertex are creases, we need to start at
        // one of them.

        for (j = 0, k = -1; j < ed.length; j++)
        {
          Edge tempEdge = e[ed[j]];
          if (tempEdge.f2 == -1 || tempEdge.smoothness < 1.0f)
          {
            if (k != -1)
              break;
            k = j;
          }
        }

        if (j == ed.length)
        {
          // There are 0 or 1 crease edges intersecting this vertex, so we will use
          // the same normal for every face.  Find it by averaging the normals of all
          // the faces sharing this point.

          Vec3 temp = new Vec3();
          int faceIndex = -1;
          for (j = 0; j < ed.length; j++)
          {
            Edge tempEdge = e[ed[j]];
            faceIndex = (tempEdge.f1 == faceIndex ? tempEdge.f2 : tempEdge.f1);
            int otherFace = (tempEdge.f1 == faceIndex ? tempEdge.f2 : tempEdge.f1);
            tempFace = f[faceIndex];
            Vec3 edge1 = v[tempFace.v2].r.minus(v[tempFace.v1].r);
            Vec3 edge2 = v[tempFace.v3].r.minus(v[tempFace.v1].r);
            Vec3 edge3 = v[tempFace.v3].r.minus(v[tempFace.v2].r);
            if (edge1.length2() < 1e-20 || edge2.length2() < 1e-20 || edge3.length2() < 1e-20)
              continue;
            edge1.normalize();
            edge2.normalize();
            edge3.normalize();
            double dot;
            if (tempFace.v1 == i)
            {
              facenorm[faceIndex*3] = normals;
              dot = edge1.dot(edge2);
            }
            else if (tempFace.v2 == i)
            {
              facenorm[faceIndex*3+1] = normals;
              dot = -edge1.dot(edge3);
            }
            else
            {
              facenorm[faceIndex*3+2] = normals;
              dot = edge2.dot(edge3);
            }
            if (dot < -1.0)
              dot = -1.0; // This can occasionally happen due to roundoff error
            if (dot > 1.0)
              dot = 1.0;
            temp.add(trueNorm[faceIndex].times(Math.acos(dot)));
            if (otherFace != -1)
            {
              tempFace = f[otherFace];
              if (tempFace.v1 == i)
                facenorm[otherFace*3] = normals;
              else if (tempFace.v2 == i)
                facenorm[otherFace*3+1] = normals;
              else
                facenorm[otherFace*3+2] = normals;
            }
          }
          temp.normalize();
          norm.addElement(temp);
          normals++;
          continue;
        }

        // This vertex is intersected by at least two crease edges, so we need to
        // calculate a normal vector for each group of faces between two creases.

        first = j = k;
        Edge tempEdge = e[ed[j]];
        groups:     do
        {
          Vec3 temp = new Vec3();
          do
          {
            // For each group of faces, find the first and last edges.  Average
            // the normals of the faces in between, and record that these faces
            // will use this normal.

            j = (j+1) % ed.length;
            m = tempEdge.f1;
            tempFace = f[m];
            if (tempFace.e1 != ed[j] && tempFace.e2 != ed[j] && tempFace.e3 != ed[j])
            {
              m = tempEdge.f2;
              if (m == -1)
                break groups;
              tempFace = f[m];
            }
            Vec3 edge1 = v[tempFace.v2].r.minus(v[tempFace.v1].r);
            Vec3 edge2 = v[tempFace.v3].r.minus(v[tempFace.v1].r);
            Vec3 edge3 = v[tempFace.v3].r.minus(v[tempFace.v2].r);
            edge1.normalize();
            edge2.normalize();
            edge3.normalize();
            double dot;
            if (tempFace.v1 == i)
            {
              facenorm[m*3] = normals;
              dot = edge1.dot(edge2);
            }
            else if (tempFace.v2 == i)
            {
              facenorm[m*3+1] = normals;
              dot = -edge1.dot(edge3);
            }
            else
            {
              facenorm[m*3+2] = normals;
              dot = edge2.dot(edge3);
            }
            if (dot < -1.0)
              dot = -1.0; // This can occasionally happen due to roundoff error
            if (dot > 1.0)
              dot = 1.0;
            temp.add(trueNorm[m].times(Math.acos(dot)));
            tempEdge = e[ed[j]];
          } while (tempEdge.f2 != -1 && tempEdge.smoothness == 1.0f);
          last = j;
          temp.normalize();
          norm.addElement(temp);
          normals++;
          j = first = last;
          tempEdge = e[ed[first]];
        } while (last != k);
      }

      // Finally, assemble all the normals into an array and create the triangles.

      normalArray = new Vec3 [norm.size()];
      for (i = 0; i < normalArray.length; i++)
        normalArray[i] = (Vec3) norm.elementAt(i);
      for (i = 0; i < f.length; i++)
      {
        tempFace = mesh.getFaces()[i];
        tri[i] = texMapping.mapTriangle(tempFace.v1, tempFace.v2, tempFace.v3,
            facenorm[i*3], facenorm[i*3+1], facenorm[i*3+2], vert);
      }
    }
    else
    {
      // The mesh is not being smooth shaded, so all the normals can be set to null.

      normalArray = new Vec3 [] {null};
      for (i = 0; i < v.length; i++)
        vert[i] = v[i].r;
      for (i = 0; i < f.length; i++)
        tri[i] = texMapping.mapTriangle(f[i].v1, f[i].v2, f[i].v3, 0, 0, 0, vert);
    }


    RenderingMesh rend = new RenderingMesh(vert, normalArray, tri, texMapping, matMapping);

    rend.setParameters(mesh.getParameterValues());
    if (interactive)
      cachedMesh = new SoftReference<RenderingMesh>(rend);
    return rend;
  }

  public TriangleMesh getTriangleMesh() {
    return this.theTriangleMesh;
  }



  @Override
  public int canConvertToTriangleMesh()
  {
    return EXACTLY;
  }

  @Override
  public TriangleMesh convertToTriangleMesh(double tol)
  {
    return theTriangleMesh;
  }

  /**
   * Set the vertices that are locked in place
   * @param lockedPoints
   */
  public void setPinnedVertices(boolean[] lockedPoints) {
    for(int i = 0; i < pinnedVerts.length; i++) {
      pinnedVerts[i] = lockedPoints[i];
    }
  }

  /**
   * Get the vertices that are locked in place
   * @return
   */
  public boolean[] getPinnedVertices() {
    return pinnedVerts;
  }

  /**
   * Get the editor that is used for editing the simulation
   * @return
   */
  public ClothSimEditorWindow getEditor() {
    return editor;
  }

  /**
   * Set the editor used for editing the simulation
   * @param clothSimEditorWindow
   */
  public void setEditor(ClothSimEditorWindow clothSimEditorWindow) {
    editor = clothSimEditorWindow;

  }

  @Override
  public boolean isEditable()
  {
    return true;
  }

  @Override
  public void edit(final EditingWindow parent, ObjectInfo info, Runnable cb)
  {
    if(editor == null) {
      editor = new ClothSimEditorWindow(parent, "Simulate Cloth '"+ info.getName() +"'", info);
    }
    editor.setVisible(true);
  }

  /** When setting the texture, we need to clear the caches. */

  @Override
  public void setTexture(Texture tex, TextureMapping mapping)
  {
    super.setTexture(tex, mapping);
    cachedMesh = null;    
  }

  /** When setting the material, we need to clear the caches. */

  @Override
  public void setMaterial(Material mat, MaterialMapping map)
  {
    super.setMaterial(mat, map);
    cachedMesh = null;
  }

  /** Set the list of texture parameters for this object. */
  public void setParameters(TextureParameter param[])
  {
    super.setParameters(param);
    if(theTriangleMesh != null) {
      theTriangleMesh.setParameters(param);
    }
  }

  /** Set the list of objects defining the values of texture parameters. */ 
  public void setParameterValues(ParameterValue val[])
  {
    super.setParameterValues(val);
    if(theTriangleMesh != null) {
      theTriangleMesh.setParameterValues(val);
    }
  }

  /** Set the object defining the value of a particular texture parameter. */ 
  public void setParameterValue(TextureParameter param, ParameterValue val)
  {
    super.setParameterValue(param, val);
    theTriangleMesh.setParameterValue(param, val);
  }

  @Override
  public void writeToFile(DataOutputStream out, Scene theScene) throws IOException{
    super.writeToFile(out, theScene);
    theTriangleMesh.writeToFile(out, theScene);
    out.writeDouble(springConstant);
    out.writeDouble(dampingConstant);
    out.writeDouble(collisionDistance);
    out.writeDouble(meshTolerance);
    out.writeDouble(vertMass);
    for(int i = 0; i < pinnedVerts.length; i++) {
      out.writeBoolean(pinnedVerts[i]);
    }
  }

  /**
   * Constructor
   * Create cloth from a file or some other input stream
   * @param in
   * @param theScene
   * @throws IOException
   * @throws InvalidObjectException
   */
  public Cloth(DataInputStream in, Scene theScene) throws IOException, InvalidObjectException{
    super(in, theScene);
    theTriangleMesh = new TriangleMesh(in, theScene);
    springConstant = in.readDouble();
    dampingConstant = in.readDouble();
    collisionDistance = in.readDouble();
    meshTolerance = in.readDouble();
    vertMass = in.readDouble();
    theMasses = new Mass[theTriangleMesh.getVertexPositions().length];
    pinnedVerts = new boolean[theTriangleMesh.getVertexPositions().length];
    for(int i = 0; i < theMasses.length; i++) {
      theMasses[i] = new Mass(theTriangleMesh.getVertexPositions()[i], i, vertMass, new Vec3());
      pinnedVerts[i] = in.readBoolean();
    }

    theSprings = new Vector<Spring>();

    // Add first layer springs
    Edge[] meshEdges = theTriangleMesh.getEdges();
    for(int i = 0; i < meshEdges.length; i++) {
      Mass massA = theMasses[meshEdges[i].v1];
      Mass massB = theMasses[meshEdges[i].v2];
      Spring s = new Spring(massA, massB, massA.getPosition().distance(massB.getPosition()), springConstant);
      if (!springExists(s)) {
        theSprings.add(s);
        massA.connectToSpring(s);
        massB.connectToSpring(s);
      }
    }   
    theTriangleMesh.copyTextureAndMaterial(this);

  }

}
