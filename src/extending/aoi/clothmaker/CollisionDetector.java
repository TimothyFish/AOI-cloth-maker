/**
 * 
 */
package extending.aoi.clothmaker;

import java.util.ArrayList;
import java.util.Collection;

import artofillusion.Scene;
import artofillusion.animation.PositionTrack;
import artofillusion.animation.RotationTrack;
import artofillusion.animation.Track;
import artofillusion.math.BoundingBox;
import artofillusion.math.Mat4;
import artofillusion.math.Vec3;
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
			if(T instanceof PositionTrack) {
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
				if((T instanceof PositionTrack) || (T instanceof RotationTrack)) {
					T.apply(time);        
				}
			}

			Vec3 dynamicDirection = new Vec3(moveVec.times(-1.0));
			dynamicDirection.normalize();

			double distanceToCollision = findDistanceToCollisionPoint(newV.plus(moveVec), I, direction, collisionDistance, objectMoved(I, prevTime, time));
			if( distanceToCollision < distance) {
				lastCollisionPoint = direction.times(distanceToCollision).plus(newV).plus(moveVec);
				lastDistanceToCollision = newV.distance(lastCollisionPoint);

				return true;
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

}
