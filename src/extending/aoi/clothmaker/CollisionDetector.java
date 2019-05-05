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

import java.util.ArrayList;
import java.util.Collection;

import artofillusion.Scene;
import artofillusion.animation.PositionTrack;
import artofillusion.animation.RotationTrack;
import artofillusion.animation.Track;
import artofillusion.animation.VisibilityTrack;
import artofillusion.math.BoundingBox;
import artofillusion.math.Mat4;
import artofillusion.math.Vec3;
import artofillusion.object.Cube;
import artofillusion.object.Cylinder;
import artofillusion.object.Light;
import artofillusion.object.NullObject;
import artofillusion.object.Object3D;
import artofillusion.object.ObjectInfo;
import artofillusion.object.ReferenceImage;
import artofillusion.object.SceneCamera;
import artofillusion.object.Sphere;
import artofillusion.object.TriangleMesh;
import artofillusion.object.TriangleMesh.Face;

/**
 * Helper class with functions that can be used to detect collisions between
 * objects in a scene.
 * @author Timothy Fish
 *
 */
public class CollisionDetector {

	private Scene scene;

	private double lastDistanceToCollision;
	private Vec3 lastCollisionPoint;
	private Triangle lastCollisionTriangle;
	private static final double MESH_TOLERANCE_NOMINEE = 0.010;

	private static final double TOL = 1e-12;

	private static final int PROP_RX = 0;
	private static final int PROP_RZ = 1;
	private static final int PROP_RATIO = 2;
	private static final int PROP_HEIGHT = 3;
	public static final int TOP = 0;
	public static final int BOTTOM = 1;
	public static final int SIDE = 2;

	private static final double FLT_EPSILON = 1e-23;

	/**
	 * Constructor
	 * @param s
	 */
	public CollisionDetector(Scene s) {
		scene = s;
	}


	/**
	 * returns true if the object is a special object, like a camera,
	 * a light, or a null object.
	 */
	public boolean isSpecial(Object obj) {
		boolean ret = false;

		if((obj instanceof ObjectInfo)) {
			// recursion because an ObjectInfo is an Object that wraps an Object
			return isSpecial(((ObjectInfo) obj).getObject());
		}
		if((obj instanceof SceneCamera)){
			ret = true;
		}
		else if(obj instanceof Light){
			ret = true;
		}
		else if(obj instanceof NullObject){
			ret = true;
		}
		else if(obj instanceof ReferenceImage){
			ret = true;
		}
		else{
			ret = false;
		}

		return ret;
	}

	/**
	 * Return the bounding box of the object.
	 * @param obj
	 * @return
	 */
	public BoundingBox getBounds(ObjectInfo obj) {
		BoundingBox B;
		if(obj.isDistorted()) {
			TriangleMesh TMesh = obj.getObject().convertToTriangleMesh(ClothMakerPlugin.DEFAULT_MESH_TOLERANCE);
			if(TMesh != null) {
				B = new BoundingBox(TMesh.getVertex(0).r.x, TMesh.getVertex(0).r.y, TMesh.getVertex(0).r.z, 
						TMesh.getVertex(0).r.x,  TMesh.getVertex(0).r.y, TMesh.getVertex(0).r.z);

				for(int i = 0; i < TMesh.getVertexPositions().length; i++) {
					B.maxx = Math.max(B.maxx, TMesh.getVertex(i).r.x);
					B.minx = Math.min(B.minx, TMesh.getVertex(i).r.x);
					B.maxy = Math.max(B.maxy, TMesh.getVertex(i).r.y);
					B.miny = Math.min(B.miny, TMesh.getVertex(i).r.y);
					B.maxz = Math.max(B.maxz, TMesh.getVertex(i).r.z);
					B.minz = Math.max(B.minz, TMesh.getVertex(i).r.z);
				}
			}
			else {
				B = new BoundingBox(0,0,0,0,0,0);
			}
		}
		else {
			B = obj.getBounds().transformAndOutset(obj.getCoords().fromLocal());
		}

		return B;
	}

	/**
	 * Finds objects that are close enough they could collide.
	 * @param obj
	 * @param deformationBox
	 * @param time
	 * @param collisionDistance
	 * @param timeIncrement
	 * @return
	 */
	public Collection<ObjectInfo> findCandidateObjects(ObjectInfo obj, BoundingBox deformationBox, double time, double collisionDistance, double timeIncrement) {
		// Determine if the boundary boxes are close enough to
		// what we're looking for.
		BoundingBox OB = deformationBox;

		ArrayList<ObjectInfo> objects = new ArrayList<ObjectInfo>();

		for(ObjectInfo candidate : scene.getAllObjects()){
			if(!candidate.isVisible() || isSpecial(candidate)) {
			} // skip this object
			else if(candidate.isDistorted()) { // TODO figure out how to deal with distorted recursion
			}
			else if(candidate == obj) {
			} // no self collisions
			else
			{
				BoundingBox CB = addCollisionDistance(getBounds(candidate), collisionDistance);
				if (!OB.intersects(CB) && !objectMoved(candidate, time-timeIncrement, time)){
				}
				else 
				{
					// bounding boxes intersect
					objects.add(candidate);
				}
			}

		}
		return objects;
	}

	/**
	 * Increases the size of the boundingBox to include the collisionDistance.
	 * @param bounds
	 * @param collisionDistance
	 * @return
	 */
	private BoundingBox addCollisionDistance(BoundingBox bounds, double collisionDistance) {
		bounds.maxx += collisionDistance;
		bounds.maxy += collisionDistance;
		bounds.maxz += collisionDistance;
		bounds.minx -= collisionDistance;
		bounds.miny -= collisionDistance;
		bounds.minz -= collisionDistance;
		return bounds;
	}

	/**
	 * Returns true if the object has moved between the previous time and now.
	 * @param obj
	 * @param timePrev
	 * @param timeNow
	 * @return
	 */
	public boolean objectMoved(ObjectInfo obj, double timePrev, double timeNow) {
		// see if any movement has occurred
		boolean ret = false;
		if(timePrev <= 0 ) timePrev = 0.001;
		if(timePrev >= timeNow) return ret;

		ObjectInfo prevObj = obj.duplicate();
		ObjectInfo nowObj = obj.duplicate();
		Vec3 prevOrigin = new Vec3(prevObj.coords.getOrigin());
		Vec3 nowOrigin = new Vec3(nowObj.coords.getOrigin());


		for(Track T : prevObj.getTracks()) {
			if(T instanceof PositionTrack) {
				T.apply(timePrev);
				prevOrigin = new Vec3(prevObj.coords.getOrigin());
			}
		}
		for(Track T : nowObj.getTracks()) {
			if(T instanceof RotationTrack) {
				T.apply(timeNow);
				nowOrigin = new Vec3(nowObj.coords.getOrigin());
			}
		}

		if(!nowOrigin.equals(prevOrigin)) {
			ret = true;
		}

		return ret;
	}

	/**
	 * Returns the vector from the previous location to the new location.
	 * @param obj
	 * @param timePrev
	 * @param timeNow
	 * @return
	 */
	public Vec3 objectMovement(ObjectInfo obj, double timePrev, double timeNow) {
		// get the distance it has moved between time1 and time2
		Vec3 ret = new Vec3();
		if(timePrev <= 0 ) timePrev = 0.001;
		if(timePrev >= timeNow) return ret;

		ObjectInfo prevObj = obj.duplicate();
		ObjectInfo nowObj = obj.duplicate();
		Vec3 prevOrigin = new Vec3(prevObj.coords.getOrigin());
		Vec3 nowOrigin = new Vec3(nowObj.coords.getOrigin());

		for(Track T : prevObj.getTracks()) {
			if(T instanceof PositionTrack) {
				T.apply(timePrev);
				prevOrigin = new Vec3(prevObj.coords.getOrigin());
			}
		}

		for(Track T : nowObj.getTracks()) {
			if(T instanceof PositionTrack) {
				T.apply(timeNow);
				nowOrigin = new Vec3(nowObj.coords.getOrigin());
			}
		}

		ret = nowOrigin.minus(prevOrigin);
		return ret;

	}

	/**
	 * Finds the distance a point will travel before colliding with an object.
	 * @param point
	 * @param nominee
	 * @param direction
	 * @param collisionDistance
	 * @return
	 */
	public double findDistanceToCollisionPoint(Vec3 point, ObjectInfo nominee, Vec3 direction,  double collisionDistance, boolean isInMotion) {  
		double ret = Double.MAX_VALUE;

		if(nominee.getObject() instanceof Sphere) {
			ret = findDistanceToEllipsoid(point, nominee, direction, collisionDistance);
		}
		else if(nominee.getObject() instanceof Cube) {
			ret = findDistanceToCube(point, nominee, direction, collisionDistance);
		}
		else if(nominee.getObject() instanceof Cylinder && !isInMotion) {
			// TODO Figure out why findDistanceToCylinder doesn't work when the cylinder is in motion, so we can handle all with special case.
			ret = findDistanceToCylinder(point, nominee, direction, collisionDistance);
		}
		else if(nominee.getDistortedObject(MESH_TOLERANCE_NOMINEE).canConvertToTriangleMesh() == Object3D.CANT_CONVERT){
			// if we can't do anything, just quit
		}
		else {
			TriangleMesh meshB = nominee.getDistortedObject(MESH_TOLERANCE_NOMINEE).convertToTriangleMesh(MESH_TOLERANCE_NOMINEE);

			// check for collision with each face
			for(Face faceB : meshB.getFaces()){
				Mat4 fromB = nominee.getCoords().fromLocal();
				Triangle triangleB = convertFaceToTriangle(meshB, faceB, fromB);

				double currentCollisionDistance = findPointTriangleCollisionDistance(point, triangleB, direction);
				if(currentCollisionDistance < ret){
					lastCollisionTriangle = triangleB;

					// check to see if this is the first collision
					ret = currentCollisionDistance;
				}  
			}
		}

		return ret;

	}

	/**
	 * Finds distance a point will travel before colliding with the cylinder
	 * @param point
	 * @param nominee
	 * @param direction
	 * @param collisionDistance
	 * @return
	 */
	private double findDistanceToCylinder(Vec3 point, ObjectInfo nominee, Vec3 direction, double collisionDistance) {
		double ret = Double.MAX_VALUE;

		Cylinder dup = (Cylinder) nominee.getObject().duplicate();
		Mat4 fromLocal = nominee.getCoords().fromLocal();
		Mat4 toLocal = nominee.getCoords().toLocal();
		double cx = fromLocal.m14/fromLocal.m44;
		double cy = fromLocal.m24/fromLocal.m44;
		double cz = fromLocal.m34/fromLocal.m44;
		double height = (double) dup.getPropertyValue(PROP_HEIGHT)+collisionDistance*2.0; 
		double halfh = height/2.0;
		double ratio = (double) dup.getPropertyValue(PROP_RATIO);
		double rx = (double) dup.getPropertyValue(PROP_RX)+collisionDistance;
		double rz = (double) dup.getPropertyValue(PROP_RZ)+collisionDistance;
		double rx2 = rx*rx;
		double rz2 = rz*rz;
		double toprx2 = rx2*ratio*ratio;
		double sy = rx*(ratio-1.0)/height;
		double sz = rx2/rz2;
		Vec3 orig = point;
		Vec3 rdir = direction;
		Vec3 v1 = new Vec3(point);
		Vec3 v2 = new Vec3(point);
		Vec3 dir = new Vec3(rdir);
		double a, b, c, d, e, temp1, temp2, mint;
		double dist1 = Double.MAX_VALUE;
		double dist2 = Double.MAX_VALUE;

		Boolean cone = (ratio == 0.0);

		int intersections;
		int hit = -1;

		v1.set(cx-orig.x, cy-orig.y, cz-orig.z);
		toLocal.transformDirection(v1);
		v1.y -= halfh;
		dir.set(rdir);
		toLocal.transformDirection(dir);


		mint = Double.MAX_VALUE;
		if (dir.y != 0.0) {
			// See if the ray hits the top or bottom face of the cylinder.

			temp1 = v1.y/dir.y;
			if (temp1 > TOL) {
				a = temp1*dir.x - v1.x;
				b = temp1*dir.z - v1.z;
				if (a*a+sz*b*b < rx2) {
					hit = BOTTOM;
					mint = temp1;
				}
			}
			if (!cone) {
				temp1 = (v1.y+height)/dir.y;
				if (temp1 > TOL)  {
					a = temp1*dir.x - v1.x;
					b = temp1*dir.z - v1.z;
					if (a*a+sz*b*b < toprx2) {
						if (mint < Double.MAX_VALUE) {
							// The ray hit both the top and bottom faces, so we know it
							// didn't hit the sides.

							intersections = 2;
							if (temp1 < mint) {
								hit = TOP;
								dist1 = temp1;
								dist2 = mint;
							}
							else {
								dist1 = mint;
								dist2 = temp1;
							}
							v1.set(orig.x+dist1*rdir.x, orig.y+dist1*rdir.y, orig.z+dist1*rdir.z);
							v2.set(orig.x+dist2*rdir.x, orig.y+dist2*rdir.y, orig.z+dist2*rdir.z);

							if(dist1 < dist2) {
								lastCollisionPoint = v1;
								lastDistanceToCollision = dist1;
							}
							else {
								lastCollisionPoint = v2;
								lastDistanceToCollision = dist2;
							}

							ret = lastDistanceToCollision;
							return ret;
						}
						else {
							hit = TOP;
							mint = temp1;
						}
					}
				}
			}
		}

		// Now see if it hits the sides of the cylinder.

		if (sy == 0.0) {
			// A simple cylinder

			temp1 = sz*dir.z;
			temp2 = 0.0;
			d = rx;
			b = dir.x*v1.x + temp1*v1.z;
			c = v1.x*v1.x + sz*v1.z*v1.z - d*d;
		}
		else {
			temp1 = sz*dir.z;
			temp2 = sy*dir.y;
			d = rx - sy*v1.y;
			b = dir.x*v1.x + d*sy*dir.y + temp1*v1.z;
			c = v1.x*v1.x + sz*v1.z*v1.z - d*d;
		}
		dist1 = Double.MAX_VALUE;
		dist2 = mint;
		if (c > TOL){ // Ray origin is outside cylinder.

			if (b > 0.0){  // Ray points toward cylinder.

				a = dir.x*dir.x + temp1*dir.z - temp2*temp2;
				e = b*b - a*c;
				if (e >= 0.0)
				{
					temp1 = Math.sqrt(e);
					dist1 = (b - temp1)/a;
					if (dist2 == Double.MAX_VALUE)
						dist2 = (b + temp1)/a;
				}
			}
		}
		else if (c < -TOL) { // Ray origin is inside cylinder.

			a = dir.x*dir.x + temp1*dir.z - temp2*temp2;
			e = b*b - a*c;
			if (e >= 0.0)
				dist1 = (b + Math.sqrt(e))/a;
		}
		else{  // Ray origin is on the surface of the cylinder.

			if (b > 0.0){  // Ray points into cylinder.

				a = dir.x*dir.x + temp1*dir.z - temp2*temp2;
				e = b*b - a*c;
				if (e >= 0.0)
					dist1 = (b + Math.sqrt(e))/a;
			}
		}
		if (dist1 < mint){

			a = dist1*dir.y-v1.y;
			if (a > 0.0 && a < height){

				hit = SIDE;
				mint = dist1;
			}
		}
		if (mint == Double.MAX_VALUE)
			return Double.MAX_VALUE;
		if (dist2 < mint){

			temp1 = dist2;
			dist2 = mint;
			mint = temp1;
		}
		dist1 = mint;
		v1.set(orig.x+dist1*rdir.x, orig.y+dist1*rdir.y, orig.z+dist1*rdir.z);
		if (hit == SIDE) {
			double dx = v1.x-cx, dz = v1.z-cz;
			double r = rx + sy*(v1.y-cy+halfh);
			double scale = r/Math.sqrt(dx*dx+sz*dz*dz);
			v1.set(cx+dx*scale, v1.y, cz+dz*scale);
		}
		if (dist2 == Double.MAX_VALUE) {
			intersections = 1;
		}
		else
		{
			intersections = 2;
			v2.set(orig.x+dist2*rdir.x, orig.y+dist2*rdir.y, orig.z+dist2*rdir.z);

		}

		if((intersections == 1) || (dist1 < dist2)) {
			lastCollisionPoint = v1;
			lastDistanceToCollision = TOL;
		}
		else {
			lastCollisionPoint = v2;
			lastDistanceToCollision = TOL;
		}

		ret = lastDistanceToCollision;
		return ret;
	}

	/**
	 * Find distance point will travel before colliding with the ellipsoid
	 * @param point
	 * @param nominee
	 * @param direction
	 * @param collisionDistance
	 * @return
	 */
	private double findDistanceToEllipsoid(Vec3 point, ObjectInfo nominee, Vec3 direction, double collisionDistance) {
		double ret = Double.MAX_VALUE;
		Sphere localSphere = (Sphere) nominee.getObject().duplicate();
		Mat4 fromLocal = nominee.getCoords().fromLocal();
		double cx = fromLocal.m14/fromLocal.m44;
		double cy = fromLocal.m24/fromLocal.m44;
		double cz = fromLocal.m34/fromLocal.m44;
		Vec3 rayDir = direction;

		Vec3 ellipsoidCenter = new Vec3(cx, cy, cz);
		Vec3 origin = point.minus(ellipsoidCenter);  

		origin = fromLocal.times(origin);
		rayDir = fromLocal.times(rayDir);

		// raySphere is the vector from the ray to the sphere center
		Vec3 raySphere = origin.times(-1);

		lastCollisionPoint = ellipsoidCenter.plus(raySphere);

		Vec3 surfacePoint = new Vec3(raySphere);
		surfacePoint.normalize();
		surfacePoint.multiply(localSphere.getRadii());
		lastDistanceToCollision = ellipsoidCenter.distance(point)-surfacePoint.length();

		ret = lastDistanceToCollision;
		return ret;
	}

	/**
	 * Find distance point will travel before colliding with the cube
	 * @param point
	 * @param nominee
	 * @param direction
	 * @param collisionDistance
	 * @return
	 */
	private double findDistanceToCube(Vec3 point, ObjectInfo nominee, Vec3 direction, double collisionDistance) {
		// TODO Auto-generated method stub
		Mat4 toLocal = nominee.getCoords().toLocal();

		Cube localCube = (Cube)nominee.getObject().duplicate();
		BoundingBox B = localCube.getBounds();
		double retval = Double.MAX_VALUE;

		Vec3 localPt = new Vec3(toLocal.times(point));
		Vec3 localDir = new Vec3(toLocal.times(direction));
		localDir.normalize();

		Triangle T[] = new Triangle[12];
		T[0] = new Triangle(new Vec3(B.minx, B.miny, B.minz), new Vec3(B.maxx, B.miny, B.minz), new Vec3(B.minx, B.miny, B.maxz));
		T[1] = new Triangle(new Vec3(B.maxx, B.miny, B.minz), new Vec3(B.minx, B.miny, B.maxz), new Vec3(B.maxx, B.miny, B.maxz));
		T[2] = new Triangle(new Vec3(B.minx, B.maxy, B.minz), new Vec3(B.maxx, B.maxy, B.minz), new Vec3(B.minx, B.maxy, B.maxz));
		T[3] = new Triangle(new Vec3(B.maxx, B.maxy, B.minz), new Vec3(B.minx, B.maxy, B.maxz), new Vec3(B.maxx, B.maxy, B.maxz));
		T[4] = new Triangle(new Vec3(B.minx, B.miny, B.minz), new Vec3(B.maxx, B.miny, B.minz), new Vec3(B.minx, B.maxy, B.minz));
		T[5] = new Triangle(new Vec3(B.maxx, B.miny, B.minz), new Vec3(B.minx, B.maxy, B.minz), new Vec3(B.maxx, B.maxy, B.minz));
		T[6] = new Triangle(new Vec3(B.minx, B.miny, B.maxz), new Vec3(B.maxx, B.miny, B.maxz), new Vec3(B.maxx, B.maxy, B.maxz));
		T[7] = new Triangle(new Vec3(B.minx, B.miny, B.maxz), new Vec3(B.maxx, B.maxy, B.maxz), new Vec3(B.minx, B.maxy, B.maxz));
		T[8] = new Triangle(new Vec3(B.maxx, B.miny, B.minz), new Vec3(B.maxx, B.miny, B.maxz), new Vec3(B.maxx, B.maxy, B.maxz));
		T[9] = new Triangle(new Vec3(B.maxx, B.miny, B.minz), new Vec3(B.maxx, B.maxy, B.minz), new Vec3(B.maxx, B.maxy, B.maxz));
		T[10] = new Triangle(new Vec3(B.minx, B.miny, B.minz), new Vec3(B.minx, B.miny, B.maxz), new Vec3(B.minx, B.maxy, B.maxz));
		T[11] = new Triangle(new Vec3(B.minx, B.miny, B.minz), new Vec3(B.minx, B.maxy, B.minz), new Vec3(B.minx, B.maxy, B.maxz));

		for(int i = 0; i < 12; i++) {

			double dist = findPointTriangleCollisionDistance(localPt, T[i], localDir);
			if(dist < retval) {
				lastCollisionPoint = nominee.getCoords().fromLocal().times(localDir.times(dist-collisionDistance).plus(localPt));
				retval = lastCollisionPoint.distance(point);
			}
		}

		return retval;
	}


	/**
	 * Return true if a collision will occur with an object.
	 * @param dir
	 * @param newV
	 * @param candidate_objects
	 * @param distance
	 * @param collisionDistance
	 * @return
	 */
	public boolean detectObjectCollision(Vec3 dir, Vec3 newV, Collection<ObjectInfo> candidate_objects, double distance, double collisionDistance) {
		Vec3 direction = new Vec3(dir);
		direction.normalize();
		for(ObjectInfo I : candidate_objects) {
			double distanceToCollision = findDistanceToCollisionPoint(newV, I, direction, collisionDistance, false);
			if( distanceToCollision < distance) {
				lastDistanceToCollision = TOL;
				lastCollisionPoint = direction.times(lastDistanceToCollision).plus(newV);

				return true;
			}

		}
		return false;
	}

	/**
	 * Return true if a collision will occur with an object.
	 * @param dir
	 * @param newV
	 * @param candidate_objects
	 * @param time
	 * @param distance
	 * @param collisionDistance
	 * @return
	 */
	public boolean detectObjectCollision(Vec3 dir, Vec3 newV, Collection<ObjectInfo> candidate_objects, double time, double distance, double collisionDistance) {
		Vec3 direction = new Vec3(dir);
		direction.normalize();
		double prevTime = time-(1.0/3000.0); // TODO base on proper sim inputs
		Vec3 moveVec = new Vec3();
		for(ObjectInfo I : candidate_objects) {
			moveVec = objectMovement(I, prevTime, time);
			for(Track T : I.getTracks()) {
				if((T instanceof PositionTrack) || (T instanceof RotationTrack) || (T instanceof VisibilityTrack)) {
					T.apply(time);        
				}
			}
			if(I.isVisible()) {

				Vec3 dynamicDirection = new Vec3(moveVec.times(-1.0));
				dynamicDirection.normalize();

				double distanceToCollision = findDistanceToCollisionPoint(newV.plus(moveVec), I, direction, collisionDistance, objectMoved(I, prevTime, time));
				if( distanceToCollision < distance) {
					lastCollisionPoint = direction.times(distanceToCollision).plus(newV).plus(moveVec);
					lastDistanceToCollision = newV.distance(lastCollisionPoint);

					return true;
				}      
			}
		}
		return false;
	}

	/**
	 * Returns the distance that was found the last time a distance to collision was calculated.
	 * @return
	 */
	public double getLastDistanceToCollision() { return lastDistanceToCollision; }

	/**
	 * Returns the point of collision that was calculated last.
	 * @return
	 */
	public Vec3 getLastCollisionPoint() {
		Vec3 ps = lastCollisionPoint;
		ps.x = Math.round(ps.x*100.0)/100.0;
		ps.y = Math.round(ps.y*100.0)/100.0;
		ps.z = Math.round(ps.z*100.0)/100.0;
		return ps;
	}
	/**
	 * Returns the triangle that was found the last time a collision was found.
	 * @return
	 */
	public Triangle getLastCollisionTriangle() { return lastCollisionTriangle; }

	/**
	 * Finds the nearest triangle to point V.
	 * @param candidate
	 * @param V
	 * @param time
	 * @return
	 */
	public Triangle findNearestTriangle(ObjectInfo candidate, Vec3 V, double time) {
		Triangle T = null;

		// if we can't do anything, just quit
		if(candidate.getDistortedObject(MESH_TOLERANCE_NOMINEE).canConvertToTriangleMesh() == Object3D.CANT_CONVERT){ return T; }

		TriangleMesh meshB = candidate.getDistortedObject(MESH_TOLERANCE_NOMINEE).convertToTriangleMesh(MESH_TOLERANCE_NOMINEE);

		double leastDistance = Double.MAX_VALUE;
		// check distance of each face
		for(Face faceB : meshB.getFaces()){
			Mat4 fromB = candidate.getCoords().fromLocal();
			Triangle triangleB = convertFaceToTriangle(meshB, faceB, fromB);

			double triangleLeast = distanceToNearestTriangleVertex(V, triangleB);
			if(triangleLeast < leastDistance) {
				T = triangleB;
				leastDistance = triangleLeast;
			}  
		}

		return T;
	}

	/**
	 * Find the vertex of the triangle that is nearest to V.
	 * @param V
	 * @param triangleB
	 * @return
	 */
	public double distanceToNearestTriangleVertex(Vec3 V, Triangle triangleB) {
		double triangleLeast = Math.min(V.distance(triangleB.getP0()), V.distance(triangleB.getP1()));
		triangleLeast = Math.min(triangleLeast, V.distance(triangleB.getP2()));
		return triangleLeast;
	}


	/**
	 * Find the distance the point is from the plane moving in this direction.
	 * @param ptA
	 * @param plane
	 * @param direction
	 * @return
	 */
	public double findPointPlaneCollisionDistance(Vec3 ptA,
			Triangle plane, Vec3 direction) {


		Vec3 b = plane.getP0();
		Vec3 normal = plane.getNormal();

		double dist = 0;
		double denominator = direction.dot(normal);
		double numerator = 0;

		// point A0
		Vec3 a0 = ptA;
		if(denominator != 0) {
			numerator = a0.minus(b).dot(normal);
			dist = Math.abs(numerator/denominator);
		}

		return dist;
	}


	/**
	 * Returns true if the cloth collides with itself. Neighboring vertices are ignored.
	 * @param mesh
	 * @param pt
	 * @param ptRadius
	 * @return
	 */
	public boolean detectSelfCollision(Cloth mesh, int pt, double ptRadius) {
		Mass point = mesh.getMasses()[pt];
		for(int i = 0; i < mesh.getMasses().length; i++) {
			// i is the point or a vertex connected to the point then ignore how close they are
			Mass point2 = mesh.getMasses()[i];
			if(i == pt) continue;

			boolean matched = false;
			for(int j = 0; j < point.getSprings().size(); j++) {
				Spring S = point.getSprings().elementAt(j);

				if(S.getMassA() == point2 || S.getMassB() == point2) {
					matched = true;
				}
			}

			if(!matched && Math.abs(point2.getPosition().distance(point.getPosition())) <= (ptRadius*2.0)) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Given a face defined in a Triangle Mesh, creates a triangle
	 * with point values transformed from local to by the matrix.
	 * @param mesh
	 * @param face
	 * @param matrix
	 * @return
	 */
	public Triangle convertFaceToTriangle(TriangleMesh mesh, Face face, Mat4 matrix) {
		Vec3 p0 = new Vec3(mesh.getVertexPositions()[face.v1]);
		matrix.transform(p0);
		Vec3 p1 = new Vec3(mesh.getVertexPositions()[face.v2]);
		matrix.transform(p1);
		Vec3 p2 = new Vec3(mesh.getVertexPositions()[face.v3]);
		matrix.transform(p2);

		Triangle triangle = new Triangle(p0, p1, p2);
		return triangle;
	}

	/**
	 * Given a point and a triangle in the same coordinate system, finds the
	 * distance between them along the direction vector. If ptA
	 * will not collide with triangleB in that direction, this method
	 * will return Double.MAX_VALUE. 
	 * @param ptA
	 * @param triangleB
	 * @param direction
	 * @return
	 */
	public double findPointTriangleCollisionDistance(Vec3 ptA,
			Triangle triangleB, Vec3 direction) {


		Vec3 b = triangleB.getP0();
		Vec3 normal = triangleB.getNormal();

		double dist = 0;
		double denominator = direction.dot(normal);
		double numerator = 0;

		if(denominator == 0) { return Double.MAX_VALUE; } // exit now, direction parallel to plane

		// point A0
		Vec3 a0 = ptA;
		numerator = a0.minus(b).dot(normal);
		dist = Math.abs(numerator/denominator);

		Vec3 P0 = a0.plus(direction.times(dist));

		if(!pointInTriangle(triangleB, P0, direction)){
			dist = Double.MAX_VALUE;
		}
		else {
			lastCollisionTriangle = triangleB;
		} 

		return dist;
	}

	/**
	 * Returns true of a line given by P and direction passes through the triangle.
	 * @param triangleB
	 * @param P
	 * @param direction
	 * @return
	 */
	public boolean pointInTriangle(Triangle triangleB, Vec3 P, Vec3 direction) {
		// Transpose triangleB onto the plane perpendicular to direction.
		Mat4 xRotMatrix = Mat4.xrotation(calculateXRotate(direction));
		Mat4 yRotMatrix = Mat4.yrotation(calculateYRotate(direction));
		Mat4 zRotMatrix = Mat4.zrotation(calculateZRotate(direction));
		Mat4 rotationMatrix = xRotMatrix.times(yRotMatrix).times(zRotMatrix);

		Triangle projB = new Triangle(
				rotationMatrix.times(triangleB.getP0()), 
				rotationMatrix.times(triangleB.getP1()), 
				rotationMatrix.times(triangleB.getP2())); 

		// Transpose P onto the plane perpendicular to direction.
		Vec3 projP = rotationMatrix.times(P);

		// Determine if P is on the same side of each edge as the opposite point.
		if(SameSide(projP, projB.getP0(), projB.getP1(), projB.getP2())
				&& SameSide(projP, projB.getP1(), projB.getP2(), projB.getP0())
				&& SameSide(projP, projB.getP2(), projB.getP0(), projB.getP1())){
			return true;
		}
		// else
		return false;
	}

	/**
	 * Returns true if the two points are on the same side of the line
	 * given by linePoint1 and linePoint2.
	 * @param point1
	 * @param point2
	 * @param linePoint1
	 * @param linePoint2
	 * @return
	 */
	public boolean SameSide(Vec3 point1, Vec3 point2, Vec3 linePoint1, Vec3 linePoint2) {
		// Two points are on the same side of a line if the dot product of the
		// cross product of each point minus a point on the line and a second
		// point on the line minus the first point is >= 0.
		Vec3 crossPoint1 = (linePoint2.minus(linePoint1)).cross(point1.minus(linePoint1));
		Vec3 crossPoint2 = (linePoint2.minus(linePoint1)).cross(point2.minus(linePoint1));
		if(crossPoint1.dot(crossPoint2) >= 0){
			return true;
		}
		//else

		return false;
	}

	/**
	 * Calculates the second rotation, which will align the direction
	 * vector with the y axis.
	 * @param x
	 * @param y
	 * @param z
	 * @return
	 */
	private double calculateZRotate(Vec3 direction) {
		// Note that this is essentially the same code as we used
		// for the Point At Plugin
		double x = direction.x;
		double y = direction.y;
		double z = direction.z;
		double sqrtYYZZ = Math.sqrt(y*y + z*z);
		double ret = 0; // return value;
		if(sqrtYYZZ != 0){

			ret = Math.toDegrees(Math.atan(x/sqrtYYZZ));

			if(sqrtYYZZ < 0){
				// angle is on the other size
				ret += 180.0; 
			}
		}
		// Handle differently because the angle gives two possible directions.
		else if(x > 0){
			ret = 90.0;
		}
		else if(x < 0){
			ret = -90.0;
		}
		return ret;
	}

	/**
	 * Returns a valid angle for Y.
	 * @return
	 */
	private double calculateYRotate(Vec3 direction) {
		// Note that this is essentially the same code as we used
		// for the Point At Plugin

		// Zero is as good as any value.
		return 0.0;
	}

	/**
	 * Calculates rotation that will put direction vector in the XY plane. 
	 * @param y
	 * @param z
	 * @return
	 */
	private double calculateXRotate(Vec3 direction) {
		// Note that this is essentially the same code as we used
		// for the Point At Plugin
		double y = direction.y;
		double z = direction.z;
		double ret = 0;

		if(y != 0){
			ret = Math.toDegrees(Math.atan(-z/y));

			if(y < 0){
				ret += 180.0;
			}
		}
		// Handle differently because the angle gives two possible directions.
		else if(z > 0){
			ret = -90.0;
		}
		else if(z < 0){
			ret = 90.0;
		}
		return  ret;
	}

	// Begin Blender stuff------------------------

	/**
	 * Collision Detection code based on Blender's Collision Detection
	 * @param dir
	 * @param newV
	 * @param candidate_objects
	 * @param time
	 * @param distance
	 * @param collisionDistance
	 * @return
	 */

	public class Compute_collision_point_return{
		double dist = 0.0;
		Vec3 r_a = new Vec3();
		Vec3 r_b = new Vec3();
		Vec3 r_vec = new Vec3();
	}
	public Compute_collision_point_return compute_collision_point(Triangle A, Triangle B, Boolean culling, Boolean use_normal) {
		Compute_collision_point_return ret = new Compute_collision_point_return();
				
		Vec3[] AVecs = new Vec3[4];
		AVecs[0] = A.getP0();
		AVecs[1] = A.getP1();
		AVecs[2] = A.getP2();
		
		Vec3[] BVecs = new Vec3[4];
		BVecs[0] = B.getP0();
		BVecs[1] = B.getP1();
		BVecs[2] = B.getP2();

		//	public boolean compute_collision_point(Vec3 dir, 
		//			                                   Vec3 newV, 
		//			                                   Collection<ObjectInfo> candidate_objects, 
		//			                                   double time, 
		//			                                   double distance, 
		//			                                   double collisionDistance) {

		Triangle a = new Triangle(A);
		Triangle b = new Triangle(B);
		double FLT_MAX = Double.MAX_VALUE;
		double dist = FLT_MAX;
		Vec3 tmp_co1 = new Vec3();
		Vec3 tmp_co2 = new Vec3();
		Vec3 isect_a = new Vec3();
		Vec3 isect_b = new Vec3();
		int isect_count = 0;
		Isect_line_segment_tri_v3_return tmp = new Isect_line_segment_tri_v3_return();
		Vec3 tmp_vec = new Vec3();
		Vec3 normal = new Vec3();
		Vec3 cent = new Vec3();
		Boolean backside = false;

		/* Find intersections. */
		for(int i = 0; i < 3; i++) {
			tmp = isect_line_segment_tri_v3(AVecs[i], AVecs[next_ind(i)], b);
			if (tmp.values_set) {
				double t = tmp.r_lambda;
				double s = 1.0f - t;
				isect_a = new Vec3( 
						s * AVecs[i].x + t * AVecs[next_ind(i)].x, 
						s * AVecs[i].y + t * AVecs[next_ind(i)].y, 
						s * AVecs[i].z + t * AVecs[next_ind(i)].z 
						);
				isect_count++;
			}
		}

		if (isect_count == 0) {
			for(int i = 0; i < 3; i++) {
			tmp = isect_line_segment_tri_v3(BVecs[i], BVecs[next_ind(i)], a);
			if (tmp.values_set) {
				double t = tmp.r_lambda;
				double s = 1.0f - t;
				isect_a = new Vec3( 
						s * BVecs[i].x + t * BVecs[next_ind(i)].x, 
						s * BVecs[i].y + t * BVecs[next_ind(i)].y, 
						s * BVecs[i].z + t * BVecs[next_ind(i)].z 
				);
				isect_count++;
			}
			}
		}
		else if (isect_count == 1) {
			int i = 0;
			do {
			tmp = isect_line_segment_tri_v3(BVecs[i], BVecs[next_ind(i)], a);
			if (tmp.values_set) {
				double t = tmp.r_lambda;
				double s = 1.0f - t;
				isect_b = new Vec3( 
						s * BVecs[i].x + t * BVecs[next_ind(i)].x, 
						s * BVecs[i].y + t * BVecs[next_ind(i)].y, 
						s * BVecs[i].z + t * BVecs[next_ind(i)].z 
				);
				isect_count++;
			}
			i++;
			}while((i < 3) && (isect_count < 1));
		}


		/* Determine collision side. */
		if (culling) {
			normal = b.getP0().minus(b.getP1()).cross(b.getP1().minus(b.getP2()));
			normal.normalize();

			cent = new Vec3((b.getP0().x+b.getP1().x+b.getP2().x)/3.0,
					(b.getP0().y+b.getP1().y+b.getP2().y)/3.0,
					(b.getP0().z+b.getP1().z+b.getP2().z)/3.0);

			if (isect_count == 2) {
				backside = true;
			}
			else if (isect_count == 0) {
				tmp_vec = a.getP0().minus(cent);
				if(tmp_vec.dot(normal) < 0.0) {
					backside = true;
				}
				else {
					tmp_vec = a.getP1().minus(cent);
					if(tmp_vec.dot(normal) < 0.0) {
						backside = true;
					}
					else {
						tmp_vec = a.getP2().minus(cent);
						if(tmp_vec.dot(normal) < 0.0) {
							backside = true;
						}	
					}
				}
			}
		}
		else if (use_normal) {
			normal = b.getP0().minus(b.getP1()).cross(b.getP1().minus(b.getP2()));
			normal.normalize();
		}

		if (isect_count == 1) {
			/* Edge intersection. */
			ret.r_a  = new Vec3(isect_a);
			ret.r_b  = new Vec3(isect_b);


			if (use_normal) {
				ret.r_vec = new Vec3(normal);
			}
			else {
				ret.r_vec = ret.r_b.minus(ret.r_a);
			}

			ret.dist = 0.0;
			return ret;
		}

		if (backside) {
			double maxdist = 0.0f;
			Boolean found = false;
			Isect_ray_tri_v3_return tmp2 = new Isect_ray_tri_v3_return();
			
			/* Point projections. */
			for(int i = 0; i < 3; i++) {
			tmp2 = isect_ray_tri_v3(AVecs[i], normal, b);
			if (tmp2.data_valid) {
				if (tmp2.r_lambda > maxdist) {
					maxdist = tmp2.r_lambda;
					ret.r_a = new Vec3( AVecs[i]);
					ret.r_b = AVecs[i].plus(normal).times(tmp2.r_lambda);
					found = true;
				}
			}
			}

			normal = normal.times(-1.0);

			for(int i = 0; i < 3; i++) {
			tmp2 = isect_ray_tri_v3(BVecs[i], normal, a);
			if (tmp2.data_valid) {
				if (tmp2.r_lambda > maxdist) {
					maxdist = tmp2.r_lambda;
					ret.r_a = new Vec3( BVecs[i]);
					ret.r_b = BVecs[i].plus(normal).times(tmp2.r_lambda);
					found = true;
				}
			}
			}

			normal = normal.times(-1.0);
		
			/* Edge projections. */
			for (int i = 0; i < 3; i++) {
				Vec3 dir = new Vec3();

				tmp_vec = BVecs[next_ind(i)].minus(BVecs[i]);
				dir = tmp_vec.cross(normal);

				for (int j = 0; j < 3; j++) {
					if (isect_line_plane_v3(tmp_co1, AVecs[j], AVecs[next_ind(j)], BVecs[i], dir) &&
							point_in_slice_seg(tmp_co1, AVecs[j], AVecs[next_ind(j)]) &&
							point_in_slice_seg(tmp_co1, BVecs[i], BVecs[next_ind(i)])) {
						tmp_co2 = closest_to_line_v3(tmp_co2, tmp_co1, BVecs[i], BVecs[next_ind(i)]).r_close;
						tmp_vec = tmp_co1.minus(tmp_co2);
						tmp.r_lambda = tmp_vec.length();

						if ((tmp.r_lambda > maxdist) && (tmp_vec.dot(normal) < 0.0f)) {
							maxdist = tmp.r_lambda;
							ret.r_a = new Vec3( tmp_co1);
							ret.r_b = new Vec3( tmp_co2);
							found = true;
						}
					}
				}
			}

			/* If no point is found, will fallback onto regular proximity test below. */
			if (found) {
				ret.r_vec = ret.r_b.minus(ret.r_a);

				if (use_normal) {
					if (normal.dot(ret.r_vec) >= 0.0f) {
						ret.r_vec = new Vec3( normal);
					}
					else {
						ret.r_vec = normal.times(-1.0);
					}
				}

				ret.dist = 0.0;
				return ret;
			}
		}

				
		/* Closest point. */
		for (int i = 0; i < 3; i++) {
			tmp_co1 = closest_on_tri_to_point_v3(AVecs[i], B);
			tmp.r_lambda = AVecs[i].minus(tmp_co1).length2();

			if (tmp.r_lambda < dist) {
				dist = tmp.r_lambda;
				ret.r_a = new Vec3( AVecs[i]);
				ret.r_b = new Vec3( tmp_co1);
			}
		}

		for (int i = 0; i < 3; i++) {
			tmp_co1 = closest_on_tri_to_point_v3(BVecs[i], A);
			tmp.r_lambda = BVecs[i].minus(tmp_co1).length2();

			if (tmp.r_lambda < dist) {
				dist = tmp.r_lambda;
				ret.r_a = new Vec3( tmp_co1);
				ret.r_b = new Vec3( BVecs[i]);
			}
		}

		/* Closest edge. */
		if (isect_count == 0) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					Vec3[] tmp_vecs = isect_seg_seg_v3(AVecs[i], AVecs[next_ind(i)], BVecs[j], BVecs[next_ind(j)]);
					tmp_co1 = tmp_vecs[0];
					tmp_co2 = tmp_vecs[1];
					
					tmp.r_lambda = tmp_co2.minus(tmp_co1).length2();

					if (tmp.r_lambda < dist) {
						dist = tmp.r_lambda;
						ret.r_a = new Vec3(tmp_co1);
						ret.r_b = new Vec3(tmp_co2);
					}
				}
			}
		}

		if (isect_count == 0) {
			ret.r_vec = ret.r_a.minus(ret.r_b);
			dist = Math.sqrt(dist);
		}
		else {
			ret.r_vec = ret.r_b.minus(ret.r_a);
			dist = 0.0f;
		}

		if (culling && use_normal) {
			ret.r_vec = new Vec3(normal);
		}
		else if (use_normal) {
			if (normal.dot(ret.r_vec) >= 0.0f) {
				ret.r_vec = new Vec3(normal);
			}
			else {
				ret.r_vec = normal.times(-1.0);
			}
		}
		else if (culling && (ret.r_vec.dot(normal) < 0.0f)) {
			ret.dist = Double.MAX_VALUE;
			return ret;
		}

		ret.dist = dist;
		return ret;
	}


	private boolean point_in_slice_seg(Vec3 p, Vec3 l1, Vec3 l2) {
		
		Vec3 normal = l2.minus(l1);
		
		double h = 0.0;
		Vec3 rp = p.minus(l1);
		h = normal.dot(rp) / normal.dot(normal);
		
		if(h < 0.0 || h > 1.0) {
			return false;
		}
		
		return true;
	}

	private class Isect_ray_tri_v3_return{
		public Boolean data_valid = false;
		public double r_lambda = 0.0;
		public double r_uv[] = new double[2];
	}
	private Isect_ray_tri_v3_return isect_ray_tri_v3(Vec3 ray_origin, Vec3 ray_direction, Triangle A) {
		Isect_ray_tri_v3_return ret = new Isect_ray_tri_v3_return();
		double epsilon = 0.00000001f;
		Vec3 p = new Vec3();
		Vec3 s = new Vec3();
		Vec3 e1 = new Vec3();
		Vec3 e2 = new Vec3();
		Vec3 q = new Vec3();
		double a, f, u, v;

		e1 = A.getP1().minus(A.getP0());
		e2 = A.getP2().minus(A.getP0());

		p = ray_direction.cross(e2);
		a = e1.dot(p);
		if ((a > -epsilon) && (a < epsilon)) {
			ret.data_valid = false;
			return ret;
		}
		f = 1.0f / a;

		s = ray_origin.minus(A.getP0());

		u = f * s.dot(p);
		if ((u < 0.0f) || (u > 1.0f)) {
			ret.data_valid = false;
			return ret;
		}

		q = s.cross(e1);

		v = f * ray_direction.dot(q);
		if ((v < 0.0f) || ((u + v) > 1.0f)) {
			ret.data_valid = false;
			return ret;
		}

		ret.r_lambda = f * e2.dot(q);
		if ((ret.r_lambda < 0.0f)) {
			ret.data_valid = false;
			return ret;
		}

		ret.r_uv[0] = u;
		ret.r_uv[1] = v;
		ret.data_valid = true;

		return ret;
	}

	Boolean isect_line_plane_v3(Vec3 r_isect_co, Vec3 l1, Vec3 l2, Vec3 plane_co, Vec3 plane_no)
	{
		Vec3 u;
		Vec3 h;
		double dot;

		u = l2.minus(l1);
		h = l1.minus(plane_co);
		dot = plane_no.dot(u);

		if (Math.abs(dot) > FLT_EPSILON) {
			double lambda = -(plane_no.dot(h)) / dot;
			r_isect_co = l1.plus(u).times(lambda);
			return true;
		}
		else {
			/* The segment is parallel to plane */
			return false;
		}
	}
	
	/* find closest point to p on line through (l1, l2) and return lambda,
	 * where (0 <= lambda <= 1) when cp is in the line segment (l1, l2)
	 */
	private class Closest_to_line_v3_return{
		public double lambda = 0.0;
		public Vec3 r_close = new Vec3();
	}
	Closest_to_line_v3_return closest_to_line_v3(Vec3 r_close, Vec3 p, Vec3 l1, Vec3 l2)
	{
		Closest_to_line_v3_return ret = new Closest_to_line_v3_return();
	  Vec3 h;
	  Vec3 u;

	  u = l2.minus(l1);
	  h = p.minus(l1);
	  ret.lambda = u.dot(h) / u.dot(u);
	  r_close = l1.plus(u).times(ret.lambda);
	  return ret;
	}
	
	/* Adapted from "Real-Time Collision Detection" by Christer Ericson,
	 * published by Morgan Kaufmann Publishers, copyright 2005 Elsevier Inc.
	 *
	 * Set 'r' to the point in triangle (a, b, c) closest to point 'p' */
	Vec3 closest_on_tri_to_point_v3(Vec3 p, Triangle T)
	{
		Vec3 a = T.getP0();
		Vec3 b = T.getP1();
		Vec3 c = T.getP2();
	  Vec3 ab = new Vec3();
	  Vec3 ac = new Vec3();
	  Vec3 ap = new Vec3();
	  double d1;
	  double d2;
	  Vec3 bp = new Vec3();
	  double d3;
	  double d4;
	  double vc;
	  Vec3 cp = new Vec3();
	  double d5;
	  double d6;
	  double vb;
	  double va;
	  double denom;
	  double v;
	  double w;
	  Vec3 r;

	  /* Check if P in vertex region outside A */
	  ab = b.minus(a);
	  ac = c.minus(a);
	  ap = p.minus(a);
	  d1 = ab.dot(ap);
	  d2 = ac.dot(ap);
	  if (d1 <= 0.0f && d2 <= 0.0f) {
	    /* barycentric coordinates (1,0,0) */
	    r = new Vec3(a);
	    return r;
	  }

	  /* Check if P in vertex region outside B */
	  bp = p.minus(b);
	  d3 = ab.dot(bp);
	  d4 = ac.dot(bp);
	  if (d3 >= 0.0f && d4 <= d3) {
	    /* barycentric coordinates (0,1,0) */
	    r = new Vec3(b);
	    return r;
	  }
	  /* Check if P in edge region of AB, if so return projection of P onto AB */
	  vc = d1 * d4 - d3 * d2;
	  if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
	    v = d1 / (d1 - d3);
	    /* barycentric coordinates (1-v,v,0) */
	    r = a.cross(ab).times(v);
	    return r;
	  }
	  /* Check if P in vertex region outside C */
	  cp = p.minus(c);
	  d5 = ab.dot(cp);
	  d6 = ac.dot(cp);
	  if (d6 >= 0.0f && d5 <= d6) {
	    /* barycentric coordinates (0,0,1) */
	    r = new Vec3(c);
	    return r;
	  }
	  /* Check if P in edge region of AC, if so return projection of P onto AC */
	  vb = d5 * d2 - d1 * d6;
	  if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
	    w = d2 / (d2 - d6);
	    /* barycentric coordinates (1-w,0,w) */
	    r = a.cross(ac).times(w);
	    return r;
	  }
	  /* Check if P in edge region of BC, if so return projection of P onto BC */
	  va = d3 * d6 - d5 * d4;
	  if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
	    w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
	    /* barycentric coordinates (0,1-w,w) */
	    r = c.minus(b).times(w).plus(b);
	    return r;
	  }

	  /* P inside face region. Compute Q through its barycentric coordinates (u,v,w) */
	  denom = 1.0f / (va + vb + vc);
	  v = vb * denom;
	  w = vc * denom;

	  r = a.plus(ab).times(v).plus(ac.times(w));
	  return r;
	}
	
	/* Returns a point on each segment that is closest to the other. */
	Vec3[] isect_seg_seg_v3(Vec3 a0, Vec3 a1, Vec3 b0, Vec3 b1)
	{
		Vec3 [] ret = new Vec3[2];
		Vec3 r_a = new Vec3();
		Vec3 r_b = new Vec3();
	  double fac_a;
	  double fac_b;
	  Vec3 a_dir = new Vec3();
	  Vec3 b_dir = new Vec3();
    Vec3 a0b0 = new Vec3();
    Vec3 crs_ab = new Vec3();
    
	  a_dir = a1.minus(a0);
	  b_dir = b1.minus(b0);
	  a0b0 = b0.minus(a0);
	  crs_ab = b_dir.cross(a_dir);
	  double nlen = crs_ab.length2();

	  if (nlen == 0.0f) {
	    /* Parallel Lines */
	    /* In this case return any point that
	     * is between the closest segments. */
	    Vec3 a0b1 = new Vec3();
	    Vec3 a1b0 = new Vec3();
	    double len_a;
	    double len_b;
	    double fac1;
	    double fac2;
	    
	    a0b1 = b1.minus(a0);
	    a1b0 = b0.minus(a1);
	    len_a = a_dir.length2();
	    len_b = b_dir.length2();

	    if (len_a > 0) {
	      fac1 = a0b0.dot(a_dir);
	      fac2 = a0b1.dot(a_dir);
	      fac1 = CLAMP(fac1, 0.0f, len_a);
	      fac2 = CLAMP(fac2, 0.0f, len_a);
	      fac_a = (fac1 + fac2) / (2 * len_a);
	    }
	    else {
	      fac_a = 0.0f;
	    }

	    if (len_b > 0) {
	      fac1 = -a0b0.dot(b_dir);
	      fac2 = -a1b0.dot(b_dir);
	      fac1 = CLAMP(fac1, 0.0f, len_b);
	      fac2 = CLAMP(fac2, 0.0f, len_b);
	      fac_b = (fac1 + fac2) / (2 * len_b);
	    }
	    else {
	      fac_b = 0.0f;
	    }
	  }
	  else {
	    Vec3 c;
	    Vec3 cray;
	    c = crs_ab.minus(a0b0);

	    cray = c.cross(b_dir);
	    fac_a = cray.dot(crs_ab) / nlen;

	    cray = c.cross(a_dir);
	    fac_b = cray.dot(crs_ab) / nlen;

	    fac_a = CLAMP(fac_a, 0.0f, 1.0f);
	    fac_b = CLAMP(fac_b, 0.0f, 1.0f);
	  }

	  r_a = a0.plus(a_dir).times(fac_a);
	  r_b = b0.plus(b_dir).times(fac_b);
	  
	  ret[0] = r_a;
	  ret[1] = r_b;
	  return ret;
	}
	
	double CLAMP(double a, double b, double c) {
		double ret = a;
		if(a < b) {
			ret = b;
		}
		else if(a > c) {
			ret = c;
		}
		
		return ret;
	}
	
	private int next_ind(int i)
	{
		return (++i < 3) ? i : 0;
	}

	/* 
	 * test if the line starting at p1 ending at p2 intersects the triangle v0..v2
	 * return non zero if it does
	 */
	private class Isect_line_segment_tri_v3_return{
		public Boolean values_set = false;
		public double r_lambda = 0.0;
		public double r_uv[] = new double[2];
	}
	Isect_line_segment_tri_v3_return isect_line_segment_tri_v3(Vec3 p1,
			Vec3 p2,
			Triangle V)
	{
		Isect_line_segment_tri_v3_return ret = new Isect_line_segment_tri_v3_return();

		Vec3 p = new Vec3();
		Vec3 s = new Vec3();
		Vec3 d = new Vec3(p2);
		Vec3 e1 = new Vec3(V.getP1());
		Vec3 e2 = new Vec3(V.getP2());
		Vec3 q = new Vec3();
		double a, f, u, v;

		e1 = V.getP1().minus(V.getP0());
		e2 = V.getP2().minus(V.getP0());
		d = p2.minus(p1);

		p = d.cross(e2);

		a = e1.dot(p);
		if (a == 0.0f) {
			ret.values_set = false;
			return ret;
		}
		f = 1.0f / a;

		s = p1.minus(V.getP0());

		u = f * s.dot(p);
		if ((u < 0.0f) || (u > 1.0f)) {
			ret.values_set = false;
			return ret;
		}

		q = s.cross(e1);

		v = f * d.dot(q);
		if ((v < 0.0f) || ((u + v) > 1.0f)) {
			ret.values_set = false;
			return ret;
		}

		ret.r_lambda = f * e2.dot(q);
		if ((ret.r_lambda < 0.0f) || (ret.r_lambda > 1.0f)) {
			ret.values_set = false;
			return ret;
		}

		ret.r_uv[0] = u;
		ret.r_uv[1] = v;
		ret.values_set = true;
		return ret;
	}

}
