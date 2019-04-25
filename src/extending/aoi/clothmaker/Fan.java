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
import artofillusion.Scene;
import artofillusion.WireframeMesh;
import artofillusion.animation.Keyframe;
import artofillusion.math.BoundingBox;
import artofillusion.math.Vec3;
import artofillusion.object.Object3D;
import artofillusion.object.ObjectInfo;
import artofillusion.ui.ComponentsDialog;
import artofillusion.ui.EditingWindow;
import artofillusion.ui.Translate;
import artofillusion.ui.ValueField;
import buoy.widget.Widget;
import buoy.widget.WindowWidget;

/**
 * Special wind producing object that simulates wind in a cloth simulation.
 * @author Timothy Fish
 *
 */
public class Fan extends Object3D {
  static final double DEFAULT_WIND_FORCE = 0.25;
  static final double DEFAULT_WIND_FALLOFF = 10.0;
  private static BoundingBox bounds = new BoundingBox(-0.25, 0.25, -0.25, 0.25, -0.25, 0.25);
  private static WireframeMesh mesh;
  private static final int PLUG_IN_VERSION = 1;

  private Scene theScene = null;
  private double magnitude = DEFAULT_WIND_FORCE;
  private double falloff = DEFAULT_WIND_FALLOFF;
  private java.util.Random rand;

  static
  {
    Vec3 vert[];
    double r = 0.25;
    int from[];
    int to[];

    vert = new Vec3 [12];
    from = new int [14];
    to = new int [14];
    vert[0] = new Vec3(0.0, 0.0, 0.0);
    vert[1] = new Vec3(0.0, r, 0.0);
    vert[2] = new Vec3(0.0, r*0.5, r*0.1);
    vert[3] = new Vec3(0.0, r*0.5, -r*0.1);
    vert[4] = new Vec3(r*0.1, r*0.5, 0.0);
    vert[5] = new Vec3(-r*0.1, r*0.5, 0.0);

    from[0] = 0;
    to[0] = 1;

    from[1] = 1;
    to[1] = 2;
    from[2] = 1;
    to[2] = 3;
    from[3] = 1;
    to[3] = 4;
    from[4] = 1;
    to[4] = 5;

    vert[6] = new Vec3(-r*0.37, 0.0, r);
    vert[7] = new Vec3( r*0.37, 0.0, r);
    vert[8] = new Vec3(-r*1.066, 0.0, 0.0);
    vert[9] = new Vec3(-r*0.814, 0.0, -r*0.69);
    vert[10] = new Vec3(r*1.066, 0.0, 0.0);
    vert[11] = new Vec3(r*0.814, 0.0, -r*0.69);

    from[5] = 0;
    to[5] = 7;
    from[6] = 7;
    to[6] = 6;
    from[7] = 6;
    to[7] = 0;

    from[8] = 0;
    to[8] = 8;
    from[9] = 8;
    to[9] = 9;
    from[10] = 9;
    to[10] = 0;

    from[11] = 0;
    to[11] = 11;
    from[12] = 11;
    to[12] = 10;
    from[13] = 10;
    to[13] = 0;

    mesh = new WireframeMesh(vert, from, to);
  }

  /**
   * constructor
   * @param theScene 
   * @param undo 
   * 
   */
  public Fan(Scene theScene, double magnitude, double falloff) {
    this.theScene  = theScene;
    this.magnitude = magnitude;
    this.falloff = falloff;
    this.rand = new java.util.Random();
  }

  /**
   * Returns a duplicate of the existing object. 
   * @return
   */
  public Object3D duplicate()
  {
    Fan t = new Fan(this.theScene, this.magnitude, this.falloff);
    return t;
  }

  /**
   * Makes this object like the one passed in. 
   * 
   * @param obj
   */
  public void copyObject(Object3D obj)
  {
    if(obj instanceof Fan){
      Fan t = (Fan)obj;

      this.theScene = t.theScene;
      this.magnitude = t.magnitude;
      this.falloff = t.falloff;
    }
  }

  /**
   * Constructor for reading data from the file.
   * @param in
   * @param theScene
   * @throws IOException
   * @throws InvalidObjectException
   */
  public Fan(DataInputStream in, Scene theScene) throws IOException, InvalidObjectException{
    super(in, theScene);

    this.theScene = theScene;
    short version = in.readShort();
    if (version != PLUG_IN_VERSION){
      throw new InvalidObjectException("");
    }
    magnitude = in.readDouble();
    falloff = in.readDouble();

  }

  /**
   * Writes the Fan to the file.
   * @param out
   * @param theScene
   * @throws IOException
   */
  public void writeToFile(DataOutputStream out, Scene theScene) throws IOException{
    super.writeToFile(out, theScene);

    // Write Version Info
    out.writeShort(PLUG_IN_VERSION);
    out.writeDouble(magnitude);
    out.writeDouble(falloff);
  }

  @Override
  public WireframeMesh getWireframeMesh()
  {
    return mesh;
  }

  @Override
  public BoundingBox getBounds() {
    return bounds;
  }

  @Override
  public void setSize(double xsize, double ysize, double zsize) {}

  @Override
  public Keyframe getPoseKeyframe() {
    return null;
  }

  @Override
  public void applyPoseKeyframe(Keyframe k) {}

  /**
   * Given a position relative to the fan, returns the force vector produced by the
   * wind from the fan at that point.
   * @param point
   * @param posNormal
   * @return
   */
  public Vec3 getForce(Vec3 position, Vec3 posNormal) {
    Vec3 force = new Vec3(0.0, rand.nextDouble()*10.0, rand.nextDouble()*1.0);
    force.normalize();
    double dist = position.length();
    double F = magnitude*(1-(dist/falloff));

    if(F < 0) {
    	F = 0;
    }
    force = force.times(F);
    
    return force;

  }

  /**
   * Set the amount of force that comes from the fan.
   * @param magnitude
   */
  public void setMagnitude(double magnitude) {
    this.magnitude = magnitude;
  }

  /**
   * Set how quickly the wind force dissipates. A larger number
   * results in cloth farther from the fan being affected while
   * a smaller number will cause it to be more localized.
   * @param falloff
   */
  public void setFalloff(double falloff) {
    this.falloff = falloff;
  }

  @Override
  public boolean isEditable()
  {
    return true;
  }

  /** Display a window in which the user can edit this object.
  @param parent   the window from which this command is being invoked
  @param info     the ObjectInfo corresponding to this object
  @param cb       a callback which will be executed when editing is complete.  If the user
                  cancels the operation, it will not be called.
   */
  @Override
  public void edit(EditingWindow parent, ObjectInfo info, Runnable cb)
  {
    ValueField magnitudeField = new ValueField(magnitude, ValueField.POSITIVE);
    ValueField falloffField = new ValueField(falloff , ValueField.POSITIVE);
    ComponentsDialog dlg = new ComponentsDialog((WindowWidget) parent, Translate.text("Select Wind Magnitude"),
        new Widget [] {magnitudeField,falloffField}, new String [] {Translate.text("Wind Magnitude"), Translate.text("Wind Falloff")});
    if (!dlg.clickedOk())
      return;
    magnitude = magnitudeField.getValue();
    falloff = falloffField.getValue();
  }
}
