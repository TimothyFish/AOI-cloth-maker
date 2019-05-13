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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InvalidObjectException;
import java.util.Vector;

import artofillusion.LayoutWindow;
import artofillusion.Scene;
import artofillusion.animation.Track;
import artofillusion.object.ObjectInfo;
import extending.aoi.clothmaker.ClothDistortion;
import extending.aoi.clothmaker.SimFrame;

/**
 * Track that stores various aspects of the simuation of a Cloth
 * @author Timothy Fish
 *
 */
public class ClothTrack extends Track {
  ObjectInfo clothInfo;
  Vector<SimFrame> storedFrames;
  SimFrame prevDrapeFrame;
  double startTime;
  double gravity;
  int gravityAxis;
  double spring_constant;
  double damping_constant;
  double collision_distance;
  double vertex_mass;
  Scene scene;
  double tfps;
  boolean selfCollision;
  boolean floorCollision;
  double drapeFrames;
  double simFrames;

  /**
   * Constructor
   * @param info
   */
  public ClothTrack(ObjectInfo info)
  {
    super("Cloth");
    this.clothInfo = info;
    storedFrames = new Vector<SimFrame>();

    gravityAxis = ClothMakerPlugin.DEFAULT_GRAVITY_AXIS;
    tfps = ClothMakerPlugin.DEFAULT_FRAMES_PER_SECOND;
    startTime = ClothMakerPlugin.DEFAULT_START_TIME;
    collision_distance = ClothMakerPlugin.DEFAULT_COLLISION_DISTANCE;
    gravity = ClothMakerPlugin.DEFAULT_GRAVITY;
    spring_constant = ClothMakerPlugin.DEFAULT_SPRING_CONST;
    damping_constant = ClothMakerPlugin.DEFAULT_DAMPING_CONST;
    vertex_mass = ClothMakerPlugin.DEFAULT_VERTEX_MASS;
    drapeFrames = ClothMakerPlugin.DEFAULT_DRAPE_FRAMES;
    simFrames = ClothMakerPlugin.DEFAULT_SIM_FRAMES;
    selfCollision = ClothMakerPlugin.DEFAULT_SELF_COLLISION;
    floorCollision = ClothMakerPlugin.DEFAULT_FLOOR_COLLISION;

    scene = null;
  }

  /**
   * Constructor
   * @param info
   * @param s
   */
  public ClothTrack(ObjectInfo info, Scene s)
  {
    super("Cloth");
    this.clothInfo = info;
    storedFrames = new Vector<SimFrame>();

    gravityAxis = ClothMakerPlugin.DEFAULT_GRAVITY_AXIS;
    tfps = ClothMakerPlugin.DEFAULT_FRAMES_PER_SECOND*ClothSimEditorWindow.subFrames;
    startTime = ClothMakerPlugin.DEFAULT_START_TIME;
    collision_distance = ClothMakerPlugin.DEFAULT_COLLISION_DISTANCE;
    gravity = ClothMakerPlugin.DEFAULT_GRAVITY;
    spring_constant = ClothMakerPlugin.DEFAULT_SPRING_CONST;
    damping_constant = ClothMakerPlugin.DEFAULT_DAMPING_CONST;
    vertex_mass = ClothMakerPlugin.DEFAULT_VERTEX_MASS;
    drapeFrames = ClothMakerPlugin.DEFAULT_DRAPE_FRAMES;
    simFrames = ClothMakerPlugin.DEFAULT_SIM_FRAMES;
    selfCollision = ClothMakerPlugin.DEFAULT_SELF_COLLISION;
    floorCollision = ClothMakerPlugin.DEFAULT_FLOOR_COLLISION;

    scene = s;
  }


  @Override
  public void edit(LayoutWindow win) {
    ClothSimEditorWindow editor = new ClothSimEditorWindow(win, "Simulate Cloth", clothInfo);
    editor.setVisible(true);
  }

  @Override
  public void apply(double time) {
    clothInfo.addDistortion(new ClothDistortion(clothInfo, storedFrames, prevDrapeFrame, time, tfps, gravity, gravityAxis, spring_constant, damping_constant, collision_distance, vertex_mass, selfCollision, floorCollision, scene));  

  }

  @Override
  public Track duplicate(Object parent) {
    ClothTrack t = new ClothTrack((ObjectInfo) parent);

    t.name = name;
    t.enabled = enabled;
    t.quantized = quantized;

    return t;
  }

  @Override
  public void copy(Track tr) {
    ClothTrack t = (ClothTrack) tr;

    name = t.name;
    enabled = t.enabled;
    quantized = t.quantized;
  }

  @Override
  public double[] getKeyTimes() {
    return new double [0];
  }

  @Override
  public int moveKeyframe(int which, double time) {
    return -1;
  }

  @Override
  public void deleteKeyframe(int which) {}

  @Override
  public boolean isNullTrack() {

    return false;
  }

  @Override
  public void writeToStream(DataOutputStream out, Scene scene) throws IOException {

    out.writeShort(0); // Version number
    out.writeUTF(name);
    out.writeBoolean(enabled);

    out.writeInt(storedFrames.size());
    for(int i = 0; i < storedFrames.size(); i++) {
      out.writeInt(storedFrames.elementAt(i).frameNumber);
      storedFrames.elementAt(i).M.writeToFile(out, scene);
    }

    out.writeDouble(startTime);
    out.writeDouble(gravity);
    out.writeInt(gravityAxis);
    out.writeDouble(spring_constant);
    out.writeDouble(damping_constant);
    out.writeDouble(collision_distance);
    out.writeDouble(vertex_mass);
    out.writeDouble(tfps/ClothSimEditorWindow.subFrames);
    out.writeBoolean(selfCollision);
    out.writeBoolean(floorCollision);
    out.writeDouble(drapeFrames);
    out.writeDouble(simFrames);

  }

  @Override
  public void initFromStream(DataInputStream in, Scene scene) throws IOException, InvalidObjectException {

    short version = in.readShort();
    if (version != 0)
      throw new InvalidObjectException("");
    name = in.readUTF();
    enabled = in.readBoolean();

    int size = in.readInt();
    this.storedFrames.clear();
    for(int i = 0; i < size; i++) {
      int frameNum = in.readInt();
      Cloth cloth = new Cloth(in, scene);
      SimFrame SF = new SimFrame(frameNum, cloth); 
      storedFrames.add(SF);
    }

    startTime = in.readDouble();
    gravity = in.readDouble();
    gravityAxis = in.readInt();
    spring_constant = in.readDouble();
    damping_constant = in.readDouble();
    collision_distance = in.readDouble();
    vertex_mass = in.readDouble();
    tfps = in.readDouble()*ClothSimEditorWindow.subFrames;
    selfCollision = in.readBoolean();
    floorCollision = in.readBoolean();
    drapeFrames = in.readDouble();
    simFrames = in.readDouble();
    this.scene = scene;

  }

  /**
   * Calls the functions to simulate cloth at frame i.
   * @param i
   * @return
   */
  public SimFrame simulateCloth(int i) {
    SimFrame frame = new SimFrame(i, (Cloth) clothInfo.getObject());
    ClothDistortion distort = new ClothDistortion(clothInfo, storedFrames, prevDrapeFrame, (double)(i) / tfps, tfps, gravity, gravityAxis, spring_constant, damping_constant, collision_distance, vertex_mass, selfCollision, floorCollision, scene);
    frame.M = distort.transform(frame.M, frame.frameNumber);

    return frame;
  }

  /**
   * Set the parameters that are used during simuation of the cloth.
   * @param startTime2
   * @param fps2
   * @param gravity2
   * @param gravityAxis2
   * @param spring_constant2
   * @param damping_constant2
   * @param vertex_mass2
   * @param collision_distance2
   * @param selfCollision2
   * @param floorCollision2
   * @param drapeFrames2
   * @param simFrames2
   */
  public void setParams(double startTime2, double fps2, double gravity2, int gravityAxis2, double spring_constant2,
      double damping_constant2, double vertex_mass2, double collision_distance2, boolean selfCollision2, boolean floorCollision2, 
      double drapeFrames2, double simFrames2) {

    gravityAxis = gravityAxis2;
    tfps = fps2;
    startTime = startTime2;
    gravity = gravity2;
    spring_constant = spring_constant2;
    damping_constant = damping_constant2;
    collision_distance = collision_distance2;
    vertex_mass = vertex_mass2;
    selfCollision = selfCollision2;
    floorCollision = floorCollision2;
    drapeFrames = drapeFrames2;
    simFrames = simFrames2;
  }
}
