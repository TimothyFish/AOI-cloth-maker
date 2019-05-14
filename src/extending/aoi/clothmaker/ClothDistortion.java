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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Vector;

import artofillusion.Scene;
import artofillusion.animation.distortion.Distortion;
import artofillusion.math.BoundingBox;
import artofillusion.math.Mat4;
import artofillusion.math.Vec3;
import artofillusion.object.Mesh;
import artofillusion.object.Object3D;
import artofillusion.object.ObjectInfo;

/**
 * This class manipulates the vertices of a Cloth such that movement of the cloth
 * is simulated over time. The simulation frames are stored so that they can be 
 * played back during animation rather than spending the time required to calculate
 * the next frame each time the clock ticks forward. 
 * 
 * @author Timothy Fish
 *
 */
public class ClothDistortion extends Distortion {
	// The time of the distortion
  private double time;
  // The scene the object is in
  private Scene scene;
  // The ObjectInfo wrapped Cloth Object
  private ObjectInfo info;

  // Enumerated values
  public static final int X_AXIS = 0;
  public static final int Y_AXIS = 1;
  public static final int Z_AXIS = 2;
  
  // Size of a vertex regarding cloth self-collision
  final double pointRadius = 0.1;

  // Reference to Vector containing the stored frames
  private Vector<SimFrame> myStoredFrames;
  // The previous frame. During draping, we just keep one frame.
  private SimFrame myPrevDrapeMesh;
  // Direction of gravity. Possible values are X_AXIS, Y_AXIS, Z_AXIS
  private int gravityAxis;
  // Force of gravity
  private double gravity;
  // Spring constant for each spring
  private double spring_constant;
  // Damping constant
  private double damping_constant;
  // Distance at which a collision is detected
  private double collision_distance;
  // Simulation Frames per Second
  private double fps;
  // True if selfCollision should be detected
  private boolean selfCollision;
  // True if cloth should be prevented from dropping lower than 0
  private boolean floorCollision;
  // weight of a vertex
	private double vertex_mass;

  /**
   * Constructor
   * @param info
   * @param storedMeshes
   * @param prevDrapeMesh
   * @param time
   * @param fps
   * @param gravity
   * @param gravityAxis
   * @param spring_constant
   * @param damping_constant
   * @param collision_distance
   * @param selfCollision
   * @param floorCollision
   * @param scene
   */
  public ClothDistortion(ObjectInfo info, Vector<SimFrame> storedMeshes, SimFrame prevDrapeMesh, double time, double fps, double gravity, int gravityAxis, 
      double spring_constant, double damping_constant, double collision_distance, double vertex_mass, boolean selfCollision, boolean floorCollision, 
      Scene scene)  {
    this.info = info;
    this.myStoredFrames = storedMeshes;
    this.myPrevDrapeMesh = prevDrapeMesh;
    this.time = time;
    this.fps = fps;
    this.gravity = gravity;
    this.gravityAxis = gravityAxis;
    this.spring_constant = spring_constant;
    this.damping_constant = damping_constant;
    this.collision_distance = collision_distance;
    this.vertex_mass = vertex_mass;
    this.selfCollision = selfCollision;
    this.floorCollision = floorCollision;
    this.scene = scene;
  }

  @Override
  public boolean isIdenticalTo(Distortion d) {
    if (!(d instanceof ClothDistortion))
      return false;
    ClothDistortion s = (ClothDistortion) d;
    if (previous != null && !previous.isIdenticalTo(s.previous))
      return false;
    if (previous == null && s.previous != null)
      return false;
    if (time != s.time || scene != s.scene || fps != s.fps || gravity != s.gravity || spring_constant != s.spring_constant ||
        damping_constant != s.damping_constant || collision_distance != s.collision_distance || selfCollision != s.selfCollision ||
        floorCollision != s.floorCollision)
      return false;

    return true;
  }

  @Override
  public Distortion duplicate() {
    ClothDistortion d = new ClothDistortion(info, myStoredFrames, myPrevDrapeMesh, time, fps, gravity, gravityAxis, 
        spring_constant, damping_constant, collision_distance, vertex_mass,
        selfCollision, floorCollision, scene);
    if (previous != null)
      d.previous = previous.duplicate();
    return d;
  }

  @Override
  public Mesh transform(Mesh obj) {
    Mesh retObj = (Mesh)obj.duplicate();
    int frameNum = (int)(time * fps);
    if(frameNum >= myStoredFrames.size()) frameNum = myStoredFrames.size() - 1 ;

    if(!myStoredFrames.isEmpty()) {
      if(frameNum > 0) {
        SimFrame frame = myStoredFrames.lastElement();
        for(int i = myStoredFrames.size()-1; i >= 0; i--) {
          if(myStoredFrames.elementAt(i).frameNumber >= frameNum) {
            frame = myStoredFrames.elementAt(i);
          }
        }
        retObj = frame.M;
        if(obj instanceof Cloth) {
          Cloth Cth = (Cloth)obj;
          ((Object3D) retObj).copyTextureAndMaterial(Cth);
        }
      }
    }
    return retObj;
  }


  /**
   * Calculates where the vertices should be at the current time.
   * This is the primary function of the simulation.
   * @param obj
   * @param frame
   * @return
   */
  public Cloth transform(Cloth obj, int frame) {
    CollisionDetector CD = new CollisionDetector(scene);
    Mat4 fromLocal = info.getCoords().fromLocal();
    Mat4 toLocal = info.getCoords().toLocal();

    if (previous != null)
      obj = (Cloth)previous.transform(obj);

    Cloth retObj = (Cloth) obj.duplicate();
    Vec3 positions[] = retObj.getVertexPositions().clone();
    
    
    for(int n = 0; n < positions.length; n++) {
      positions[n] = fromLocal.times(positions[n]);  
    }
    retObj.setVertexPositions(positions);

    SimFrame prevSF = load_prev_mesh(frame);
    if(prevSF == null) {
      prevSF = new SimFrame(frame, retObj);
    }

    final int POINTS_TOTAL = retObj.getVertexPositions().length;

    Vec3 P[] = new Vec3[POINTS_TOTAL];
    double W[] = new double[POINTS_TOTAL];
    Vec3 V[] = new Vec3[POINTS_TOTAL];

    for(int i = 0; i < POINTS_TOTAL; i++) {
      P[i] = prevSF.M.getMasses()[i].getPosition();
      W[i] = prevSF.M.getMasses()[i].getWeight();
      V[i] = prevSF.M.getMasses()[i].getVelocity();
    }

    Vec3 newvert[] = new Vec3[POINTS_TOTAL];

    Vec3 g; // gravity
    switch(gravityAxis) {
    case X_AXIS:
      g = new Vec3(gravity, 0.0, 0.0);
      break;
    case Y_AXIS:
    default:
      g = new Vec3(0.0, gravity, 0.0);
      break;
    case Z_AXIS:
      g = new Vec3(0.0, 0.0, gravity);
      break;
    }
    double k = spring_constant; 
    double c = damping_constant;

    ArrayList<ObjectInfo> fans = new ArrayList<ObjectInfo>();

    for(ObjectInfo candidate : scene.getAllObjects()){
      if(candidate.isVisible() && (candidate.getObject() instanceof Fan)) {
        fans.add(candidate);
      }
    }

    for(int pt = 0; pt < retObj.getMasses().length; pt++) {
      Mass curMass = retObj.getMasses()[pt];
      double t = time;                      // time step
      Vec3 p = new Vec3(curMass.getPosition()); // our point
      double m = vertex_mass;
      Vec3 u = new Vec3(V[pt]); // current velocity of the mass
      Vec3 F = g.times(m).minus(u.times(c)); // F is force on mass F = gravity * m - c * u

      // Add wind forces
      for(ObjectInfo fanInfo : fans) {
      	for(int i = 0; i < fanInfo.getTracks().length; i++) {
      		fanInfo.getTracks()[i].apply(time);
      	}
      	Fan theFan = (Fan)fanInfo.getObject();
      	Vec3 fanPt = fanInfo.coords.toLocal().times(p);

      	Vec3 A = theFan.getForce(fanPt, fanInfo.coords.toLocal().times(retObj.getNormals()[pt]));
      	Vec3 forceVector = A;
     	
        double magnitude = forceVector.length();
        
        forceVector = fanInfo.getCoords().getUpDirection().times(magnitude);
        Vec3 windDelta = theFan.getWindDelta();
        
        Vec3 D = new Vec3(windDelta.x + forceVector.x, windDelta.y + forceVector.y, windDelta.z + forceVector.z);
        F = F.plus(D);
      	
      }
      
      boolean fixed_node = true;

      for(Spring curSpring: curMass.getSprings()) {
        if(!retObj.getPinnedVertices()[pt]) {
          Vec3 q;
          if(curSpring.getMassA()!=curMass) {
            q = curSpring.getMassA().getPosition();
          }
          else{
            q = curSpring.getMassB().getPosition();
          }
          Vec3 d = q.minus(p);
          double x = d.length();
          Vec3 normalizeD = new Vec3(d);
          normalizeD.normalize();
          F = F.plus(normalizeD.times(-k * (curSpring.getRestLength() - x)));

          fixed_node = false;
        }
      }


      if (fixed_node) {
        F.set(0.0, 0.0, 0.0);
      }

      // Acceleration due to force
      Vec3 a = F.times(1.0/m);

      // Displacement
      Vec3 s = u.times(t).plus(a.times(0.5*t*t));

      // final velocity
      Vec3 v = u.plus(a.times(t));

      // Constrain the absolute value of the displacement per step
      double clamp_value = 0.0045;
      s = clamp(s, clamp_value);
      Vec3 ps = p.plus(s);

      P[pt].x = ps.x;
      P[pt].y = ps.y;
      P[pt].z = ps.z;
      W[pt] = m;
      V[pt] = v;

      // maximal location the vector can move to this simulation frame
      Vec3 maxPoint = P[pt].minus(prevSF.M.getMasses()[pt].getPosition());
      maxPoint.normalize();
      maxPoint = maxPoint.times(collision_distance).plus(P[pt]);

      Collection<ObjectInfo> candidates = CD.findCandidateObjects(info, new BoundingBox(maxPoint, prevSF.M.getMasses()[pt].getPosition()), time, collision_distance, 1.0/fps);



      Vec3 prev = prevSF.M.getTriangleMesh().getVertexPositions()[pt];

      candidates = CD.findCandidateObjects(info, new BoundingBox(maxPoint, prevSF.M.getMasses()[pt].getPosition()), time/ClothSimEditorWindow.subFrames, collision_distance, 1.0/(fps/ClothSimEditorWindow.subFrames));
      if(CD.detectObjectCollision(prev, P[pt], candidates, time/ClothSimEditorWindow.subFrames, collision_distance, collision_distance))
      {
      	ps = CD.getLastCollisionPoint();
     	
      	P[pt].x = ps.x;
      	P[pt].y = ps.y;
      	P[pt].z = ps.z;
      	V[pt] = v.times(0.0);

      }
      
      if(selfCollision) {
        boolean isSelfCollision = CD.detectSelfCollision(prevSF.M, pt, pointRadius);
        if(isSelfCollision) {
          ps.x = prevSF.M.getMasses()[pt].getPosition().x;
          ps.y = prevSF.M.getMasses()[pt].getPosition().y;
          ps.z = prevSF.M.getMasses()[pt].getPosition().z;
          double moveDist = pointRadius/10000.0; 
          ps.z += moveDist;

          P[pt].x = ps.x;
          P[pt].y = ps.y;
          P[pt].z = ps.z;
          V[pt] = v.times(0.0);
        }
      }

      if(floorCollision) {
        if(P[pt].y < 0.0) {
          ps.y = collision_distance;
          P[pt].y = ps.y;
        }
      }

      newvert[pt] = new Vec3( ps.x, ps.y, ps.z);
    }

    for(int n = 0; n < newvert.length; n++) {
      newvert[n] = toLocal.times(newvert[n]);
    }

    retObj.setVertexPositions(newvert);

    save_mesh(frame, retObj);

    return retObj;
  }

  /**
   * Stores the mesh so that it can be recalled later.
   * @param frame
   * @param mesh
   */
  private void save_mesh(int frame, Cloth mesh) {
    SimFrame simFrame = new SimFrame(frame/ClothSimEditorWindow.subFrames, mesh);

    myPrevDrapeMesh = simFrame;

    if(frame >= 0 && (frame % ClothSimEditorWindow.subFrames == 0)) {
      int foundFrame = -1;
      for(int i = 0; i < myStoredFrames.size(); i++) {
        if(myStoredFrames.elementAt(i).frameNumber == frame/ClothSimEditorWindow.subFrames) {
          foundFrame = i;
          break;
        }
      }

      if(foundFrame >= 0) {
        myStoredFrames.set(foundFrame/ClothSimEditorWindow.subFrames, simFrame);
      }
      else {
        myStoredFrames.add(simFrame);
      }
    }
  }

  /**
   * Finds the frame right before the currentFrame and returns it.
   * @param currentFrame
   * @return
   */
  private extending.aoi.clothmaker.SimFrame load_prev_mesh(int currentFrame) {
    SimFrame prevFrame = null;
    if(currentFrame <= 0 || (currentFrame % ClothSimEditorWindow.subFrames != 0)) {
      prevFrame = myPrevDrapeMesh;
    }
    else {
      if(!myStoredFrames.isEmpty()) {
        prevFrame = myStoredFrames.firstElement();
        for(int i = 0; i < myStoredFrames.size(); i++) {
          if(myStoredFrames.elementAt(i).frameNumber < currentFrame) {
            prevFrame = myStoredFrames.elementAt(i);
          }
          else {
            myStoredFrames.setSize(i+1);
            break;
          }
        }

        if(prevFrame != null && prevFrame.frameNumber != currentFrame) {
          prevFrame = new SimFrame(prevFrame.frameNumber, (Cloth)prevFrame.M.duplicate());
        }
        else {
          prevFrame = null;
        }
      }
    }

    return prevFrame;

  }
  
  /**
   * Prevents s from being too long by returning a Vec3 with any values
   * that are too long dropped to within the box.
   * @param s
   * @param clampValue
   * @return
   */
  private Vec3 clamp(Vec3 s, double clampValue) {
    Vec3 ret = new Vec3(s);
    
    ret.x = (ret.x > clampValue)?clampValue:ret.x;
    ret.y = (ret.y > clampValue)?clampValue:ret.y;
    ret.z = (ret.z > clampValue)?clampValue:ret.z;
    ret.x = (ret.x < -clampValue)?-clampValue:ret.x;
    ret.y = (ret.y < -clampValue)?-clampValue:ret.y;
    ret.z = (ret.z < -clampValue)?-clampValue:ret.z;
    
    double nearlyzero = 1.0E-10;
		ret.x = (Math.abs(ret.x)<nearlyzero )?0.0:ret.x;
		ret.y = (Math.abs(ret.y)<nearlyzero )?0.0:ret.y;
		ret.z = (Math.abs(ret.z)<nearlyzero )?0.0:ret.z;

    return ret;
  }

}
