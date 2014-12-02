#pragma once

#include "aabb.h"
#include "circle.h"

//#define SHOW_LOCALSPACE_COLLISION

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

	Box(V2<TYPE> center, V2<TYPE> size, float rotationInRadians) : center(center), size(size), rotation(rotationInRadians){}
	Box(V2<TYPE> center, V2<TYPE> size, V2<TYPE> rotationVector) : center(center), size(size), rotation(rotationVector){}
	Box() :rotation(V2<TYPE>::ZERO_DEGREES()){}

	void set(V2<TYPE> const center, V2<TYPE> const size, float const rotation) {
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

#ifdef SHOW_LOCALSPACE_COLLISION
	V2<TYPE> point;
#endif
	AABB<TYPE> getLocalSpaceAABB() const {
		V2<TYPE> halfDiagonal = size / 2;
		return AABB<TYPE>(-halfDiagonal, halfDiagonal);
	}

	bool intersects(Circle<TYPE> const c)
#ifndef SHOW_LOCALSPACE_COLLISION
		const
#endif
	{
		AABB<TYPE> aabb = getLocalSpaceAABB();
		V2<TYPE> yAxis, xAxis;
		calculateAxis(xAxis, yAxis);
#ifndef SHOW_LOCALSPACE_COLLISION
		V2<TYPE> point;
#endif
		point = c.center.toLocalSpace(xAxis, yAxis, center);
		return Circle<TYPE>::intersects(point, c.radius, aabb.min, aabb.max);
	}


	Circle<TYPE> getCollisionCircle() const {
		return Circle<TYPE>(center, size.magnitude() / 2);
	}

	V2<TYPE> getClosestPointOnEdge(V2<TYPE> const point, V2<TYPE> & out_normal) const {
		AABB<TYPE> aabb = getLocalSpaceAABB();
		V2<TYPE> yAxis, xAxis;
		calculateAxis(xAxis, yAxis);
		// do math in local space
		V2<TYPE> close;
		V2f rotatedPoint = point.toLocalSpace(xAxis, yAxis, center);
		close = aabb.getClosestPointOnEdge(rotatedPoint, out_normal);
		// bring results to world space
		out_normal = out_normal.toWorldSpace(xAxis, yAxis);
		close = close.toWorldSpace(xAxis, yAxis, center);
		return close;
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
		if (c.intersects(bc)) {
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
			float dist;
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

	bool intersects(AABB<TYPE> aabb) const {
		Box<TYPE> b(aabb.getCenter(), aabb.diagonal(), 0);
		return intersects(b);
	}

	/**
	 * Does ray-casting against one polygon by detecting line collision
	 *
	 * @param rayStart
	 * @param rayDirection must be a unit vector
	 * @param out_dist how far along the ray the collision happened
	 * @param out_point where the collision happened
	 * @param out_normal the normal of the collision line
	 * @return true if raycast collision happened
	 */
	bool raycast(V2<TYPE> const & rayStart, V2<TYPE> const & rayDirection,
		float & out_dist, V2<TYPE> & out_point, V2<TYPE> & out_normal) const {
		V2<TYPE> points[4];
		writePoints(points, 4);
		return V2<TYPE>::raycastPolygon(rayStart, rayDirection, points, 4, out_dist, out_point, out_normal) != -1;
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
#ifdef SHOW_LOCALSPACE_COLLISION
		AABB<TYPE> aabb = getLocalSpaceAABB();
		aabb.writePoints(points, 4);
		glBegin(GL_LINE_LOOP);
		for (int i = 0; i < 4; ++i)
			points[i].glVertex();
		glEnd();
		glDrawCircle(point, .1f, false);
#endif
		return true;
	}
	bool glDraw() const { return glDraw(false); }
#endif
};

typedef Box<float> BoxF;

#ifdef SHOW_LOCALSPACE_COLLISION
#undef SHOW_LOCALSPACE_COLLISION
#endif
