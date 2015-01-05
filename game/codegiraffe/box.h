#pragma once

#include "aabb.h"
#include "circle.h"

/**
* A roated Bounding-Box, which includes methods to help with 2D math, including 2D collision
* Intended to be used with OpenGL, and more specifically, freeglut http://freeglut.sourceforge.net/docs/api.php
* Written to be useful, readable, and informative for students studying 2D Vector math.
*
* @author mvaganov@hotmail.com October 2014
*/
template<typename TYPE>
struct Box {
	V2<TYPE> center;
	V2<TYPE> size;
	V2<TYPE> rotation;

	Box(V2<TYPE> center, V2<TYPE> size, TYPE rotationInRadians) : center(center), size(size), rotation(rotationInRadians){}
	Box(V2<TYPE> center, V2<TYPE> size, V2<TYPE> rotationVector) : center(center), size(size), rotation(rotationVector){}
	Box(AABB<TYPE> aabb) :center(aabb.getCenter()), size(aabb.diagonal()), rotation((TYPE)0){}
	Box() :rotation(V2<TYPE>::ZERO_DEGREES()){}

	void set(V2<TYPE> const center, V2<TYPE> const size, TYPE const rotation) {
		this->center = center;
		this->size = size;
		this->rotation = V2<TYPE>(rotation);
	}

	void set(V2<TYPE> const center, V2<TYPE> const size, V2<TYPE> const rotation) {
		this->center = center;
		this->size = size;
		this->rotation = rotation;
	}

	void rotate(V2<TYPE> cos_sin){}
	void rotateAround(V2<TYPE> cos_sin, V2<TYPE> position){}

	void calculateAxis(V2<TYPE> & xAxis, V2<TYPE> & yAxis) const {
		yAxis.set(0, 1);
		xAxis.set(1, 0);
		yAxis.rotate(rotation);
		xAxis.rotate(rotation);
	}

	/**
	 * writes {(minx,miny), (minx,maxy), (maxx,maxy), (maxx,miny)}
	 */
	bool writePoints(V2<TYPE> * points, const int pointCount) const {
		if (pointCount != 4) return false;
		V2<TYPE> rad = size / 2, xAxis, yAxis;
		calculateAxis(xAxis, yAxis);
		points[0].set(center - xAxis * rad.x - yAxis * rad.y); // min
		points[1].set(center - xAxis * rad.x + yAxis * rad.y);
		points[2].set(center + xAxis * rad.x + yAxis * rad.y); // max
		points[3].set(center + xAxis * rad.x - yAxis * rad.y);
		return true;
	}

	bool writePointsRelative(V2<TYPE> * points, const int pointCount) {
		if (pointCount != 4) return false;
		V2<TYPE> rad = size / 2;
		points[0].set(-rad.x, -rad.y); // min
		points[1].set(-rad.x, +rad.y);
		points[2].set(+rad.x, +rad.y); // max
		points[3].set(+rad.x, -rad.y);
		return true;
	}

	AABB<TYPE> getLocalSpaceAABB() const {
		V2<TYPE> halfDiagonal = size / 2;
		return AABB<TYPE>(-halfDiagonal, halfDiagonal);
	}

	bool intersectsCircle(Circle<TYPE> const c) const {
		return intersectsCircle(c.center, c.radius);
	}
	bool intersectsCircle(V2<TYPE> const & a_center, const TYPE a_radius) const  {
		AABB<TYPE> aabb = getLocalSpaceAABB();
		V2<TYPE> yAxis, xAxis;
		calculateAxis(xAxis, yAxis);
		V2<TYPE> point;
		point = a_center.toLocalSpace(xAxis, yAxis, center);
		return Circle<TYPE>::intersects(point, a_radius, aabb.min, aabb.max);
	}


	Circle<TYPE> getCollisionCircle() const {
		return Circle<TYPE>(center, size.magnitude() / 2);
	}

	/**
	* @param point
	* @param out_rh will be assigned to details about the closest raycast hit to that point
	* @return true if on linearedge, false if on corner
	*/
	bool getClosestRaycastHit(V2<TYPE> const & point, RaycastHit_<TYPE> & out_rh) const {
		AABB<TYPE> aabb = getLocalSpaceAABB();
		V2<TYPE> yAxis, xAxis;
		calculateAxis(xAxis, yAxis);
		// do math in local space
		V2f rotatedPoint = point.toLocalSpace(xAxis, yAxis, center);
		bool onLine = aabb.getClosestRaycastHit(rotatedPoint, out_rh);
		// bring results to world space
		out_rh.normal = out_rh.normal.toWorldSpace(xAxis, yAxis);
		out_rh.point = out_rh.point.toWorldSpace(xAxis, yAxis, center);
		return onLine;
	}

	bool contains(V2<TYPE> point) const {
		AABB<TYPE> aabb = getLocalSpaceAABB();
		V2<TYPE> yAxis, xAxis;
		calculateAxis(xAxis, yAxis);
		V2f rotatedPoint = point.toLocalSpace(xAxis, yAxis, center);
		return aabb.contains(rotatedPoint);
	}

	bool contains(V2<TYPE> * points, int pointsCount) const {
		AABB<TYPE> aabb = getLocalSpaceAABB();
		V2<TYPE> yAxis, xAxis;
		calculateAxis(xAxis, yAxis);
		for (int i = 0; i < pointsCount; ++i) {
			V2f rotatedPoint = points[i].toLocalSpace(xAxis, yAxis, center);
			if (aabb.contains(rotatedPoint))
				return true;
		}
		return false;
	}

	// TODO make this faster by checking against AABB
	bool intersects(Box<TYPE> const b) const {
		Circle<TYPE> c = getCollisionCircle(), bc = b.getCollisionCircle();
		// if they are even close to the same region
		if (c.intersectsCircle(bc)) {
			const int NUM_POINTS = 4;
			// check if any points are inside the other
			V2<TYPE> ap[NUM_POINTS], bp[NUM_POINTS];
			b.writePoints(bp, NUM_POINTS);
			if (contains(bp, NUM_POINTS))
				return true;
			writePoints(ap, NUM_POINTS);
			if (b.contains(ap, NUM_POINTS))
				return true;
			// check if any lines cross
			V2<TYPE> point;
			TYPE dist;
			for (int a = 0; a < NUM_POINTS; ++a) {
				for (int b = 0; b < NUM_POINTS; ++b) {
					if (V2<TYPE>::lineIntersection(
						ap[a], ap[(a + 1) % NUM_POINTS],
						bp[b], bp[(b + 1) % NUM_POINTS],
						dist, point)) return true;
				}
			}
		}
		return false;
	}

	bool intersectsAABB(AABB<TYPE> aabb) const {
		Box<TYPE> b(aabb.getCenter(), aabb.diagonal(), 0);
		return intersects(b);
	}

	bool intersectsAABB(V2<TYPE> const & min, V2<TYPE> const & max) const {
		V2f delta = max - min;
		Box<TYPE> b(min + delta/2 , delta, 0);
		return intersects(b);
	}

	/**
	 * Does ray-casting against one polygon by detecting line collision
	 *
	 * @param rayStart
	 * @param rayDirection must be a unit vector
	 * @param out_rh information about the eventual hit (if any)
	 * @return true if raycast collision happened
	 */
	bool raycast(Ray_<TYPE> const & ray, RaycastHit & out_rh) const {
		V2<TYPE> points[4];
		writePoints(points, 4);
		return V2<TYPE>::raycastPolygon(ray.start, ray.direction, points, 4, out_rh.distance, out_rh.point, out_rh.normal) != -1;
	}

#ifdef __GL_H__
	bool glDraw(bool filled) const {
		AABB<TYPE> aabb = getLocalSpaceAABB();
		glPushMatrix();
		center.glTranslate();
		rotation.glRotate();
		aabb.glDraw(filled);
		glPopMatrix();
		//V2<TYPE> points[4];
		//writePoints(points, 4);
		//glBegin(filled ? GL_POLYGON : GL_LINE_LOOP);
		//for (int i = 0; i < 4; ++i)
		//	points[i].glVertex();
		//glEnd();
		return true;
	}
	bool glDraw() const { return glDraw(false); }
#endif
};

typedef Box<float> Boxf;
