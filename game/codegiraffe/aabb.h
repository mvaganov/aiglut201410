#pragma once
#include "v2.h"

// to get around the windows API definition of min and max macros
#define NOMINMAX
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif max

/**
 * An Axis-Aligned Bounding-Box, which includes methods to help with 2D math, including 2D collision
 * Comments assume standard graphics axis alignment: y increases as y goes down (0,0 in upper-left)
 * Intended to be used with OpenGL, and more specifically, freeglut http://freeglut.sourceforge.net/docs/api.php
 * Written to be useful, readable, and informative for students studying 2D Vector math.
 *
 * @author mvaganov@hotmail.com October 2014
 */
template<typename TYPE>
struct AABB {
	/**
	 * corners of box are stored as opposed to x/y w/h.
	 * m_ (name by scope) used to avoid collision with possible min/max macros
	 */
	V2<TYPE> min, max;

	/** @return top-left */
	const V2<TYPE> getMin() const { return min; }

	/** @return bottom-right */
	const V2<TYPE> getMax() const { return max; }

	/** @return y value of top edge */
	const TYPE getMinY() const { return min.y; }
	/** @return x value of left edge */
	const TYPE getMinX() const { return min.x; }
	/** @return y value of bottom edge */
	const TYPE getMaxY() const { return max.y; }
	/** @return x value of right edge */
	const TYPE getMaxX() const { return max.x; }

	/** set the left-edge location (resizes the rectangle) */
	void setMinX(const TYPE a_value) {	min.x = a_value; }
	/** set the right-edge location (resizes the rectangle) */
	void setMaxX(const TYPE a_value) {	max.x = a_value; }
	/** set the top-edge location (resizes the rectangle) */
	void setMinY(const TYPE a_value) {	min.y = a_value; }
	/** set the bottom-edge location (resizes the rectangle) */
	void setMaxY(const TYPE a_value) { max.y = a_value; }

	/** add to the left-edge location (resizes the rectangle) */
	void addMinX(const TYPE a_value) { min.x += a_value; }
	/** add to the right-edge location (resizes the rectangle) */
	void addMaxX(const TYPE a_value) { max.x += a_value; }
	/** add to the top-edge location (resizes the rectangle) */
	void addMinY(const TYPE a_value) { min.y += a_value; }
	/** add to the bottom-edge location (resizes the rectangle) */
	void addMaxY(const TYPE a_value) { max.y += a_value; }
	/** add values to each side */
	void add(const TYPE minX, const TYPE minY, const TYPE maxX, const TYPE maxY) {
		min.add(minX, minY);
		max.add(maxX, maxY);
	}
	void addMin(const V2<TYPE> delta) { min.add(delta); }
	void addMax(const V2<TYPE> delta) { max.add(delta); }

	bool isEqual(const AABB<TYPE> & b) const {
		return min == b.min && max == b.max;
	}

	bool isValid() const { return min.x <= max.x && min.y <= max.y; }

	bool operator==(const AABB<TYPE> & b){ return isEqual(b); }

	/**
	* writes {(minx,miny), (minx,maxy), (maxx,maxy), (maxx,miny)}
	*/
	bool writePoints(V2<TYPE> * a_array, const int pointCount) const {
		if (pointCount != 4) return false;
		a_array[0].set(min.x, min.y);
		a_array[1].set(min.x, max.y);
		a_array[2].set(max.x, max.y);
		a_array[3].set(max.x, min.y);
		return true;
	}

	/** set the upper-left corner location (resizes the rectangle) */
	void setMin(const V2<TYPE> & p) { min = p; }
	/** set the lower-right corner location (resizes the rectangle) */
	void setMax(const V2<TYPE> & p) { max = p; }

	TYPE getWidth() const { return max.x - min.x; }
	TYPE getHeight() const { return max.y - min.y; }

	/** set the width value (moves the max value, keeping min value stationary) */
	void setWidth(const TYPE a_value) {	setMaxX(getMinX()+a_value);	}
	/** set the height value (moves the max value, keeping min value stationary) */
	void setHeight(const TYPE a_value) {	setMaxY(getMinY()+a_value);	}

	/** @return x position of rectangle (left edge location) */
	const TYPE getX() const { return min.x; }
	/** @return y position of rectangle (top edge location) */
	const TYPE getY() const { return min.y; }

	/** set the x position value (moves the rectangle) */
	void setX(const TYPE a_value) {
		TYPE d = getMaxX() - getMinX();
		setMinX(a_value);	
		setMaxX(a_value+d);
	}

	/** set the y position value (moves the rectangle) */
	void setY(const TYPE a_value) {
		TYPE d = getMaxY() - getMinY();
		setMinY(a_value);
		setMaxY(a_value+d);
	}
	/** the distance vector between the two AABB points */
	V2<TYPE> diagonal() const { return getMax() - getMin(); }
	
	/** width/height */
	V2<TYPE> getDimension() const { return diagonal(); }

	void setDimension(const V2<TYPE> & a_dim) {setWidth(a_dim.getX());	setHeight(a_dim.getY());}

	TYPE getArea() const { return getWidth()*getHeight(); }

	/** re-constructor */
	void set(const V2<TYPE> & a_min, const V2<TYPE> & a_max) {	setMin(a_min);	setMax(a_max);	}

	void set(const AABB<TYPE> & r) { setMin(r.getMin());	setMax(r.getMax()); }

	/** given points in arbitrary order, makes this an inscribed rectangle */
	void setFromPoints(const V2<TYPE> & a, const V2<TYPE> & b) {
		V2<TYPE> min(((a.x<b.x)?a.x:b.x),((a.y<b.y)?a.y:b.y));
		V2<TYPE> max(((a.x>b.x)?a.x:b.x),((a.y>b.y)?a.y:b.y));
		set(min, max);
	}

	/** re-constructor */
	void set(const V2<TYPE> & a_center, const TYPE a_radius) {
		V2<TYPE> corner(a_radius, a_radius);
		set(a_center - corner, a_center + corner);
	}

	/** default constructor, everything is zero */
	AABB(){}

	/** */
	AABB(const V2<TYPE> & a_min, const V2<TYPE> & a_max) { setFromPoints(a_min, a_max); }

	/** @param where to center the box, a square a_radius*2 tall and wide */
	AABB(const V2<TYPE> & a_center, const TYPE a_radius) {
		V2<TYPE> corner(a_radius, a_radius);
		set(a_center - corner, a_center + corner);
	}

	/** @return the radius of the circle that this rectangle is circumscribed by (circle outside) */
	TYPE getCircumscribedRadius() const { return diagonal().quotient(2).magnitude(); }

	/** @return the radius of the circle that this rectangle is circumscribed by (circle inside) */
	TYPE getInscribedRadius() const { 
		TYPE w = getWidth(), h = getHeight();
		return ((w<=h)?w:h) / 2;
	}

	/** resizes the rectangle so it could be circumscribed by (circle outside) a circle with the given radius */
	void setCircumscribedRadius(const TYPE a_radius) {
		V2<TYPE> diffRad = diagonal();
		if(diffRad.isZero()) {
			diffRad.set(1,1);
		}
		TYPE currentRad = (diffRad / 2).magnitude();
		TYPE ratio = a_radius / currentRad;
		V2<TYPE> center(getCenter());
		diffRad *= ratio;
		setMin(center - diffRad);
		setMax(center + diffRad);
	}

	bool intersectsCircle(V2<TYPE> const & center, const TYPE radius) const {
		V2<TYPE> rad(radius, radius);
		V2<TYPE> outerMin = min - rad, outerMax = max + rad;
		// check if collision is possible at the most extreme estimate
		if (center.x >= outerMin.x && center.y >= outerMin.y && center.x < outerMax.x && center.y < outerMax.y) {
			// if it's within range horizontally or vertically, that's enough for line intersection.
			if ((center.y >= min.y && center.y < max.y) || (center.x >= min.x && center.x < max.x))
				return true;
			// otherwise, find the closest corner
			V2<TYPE> closeCorner;
			closeCorner.x = (center.x < min.x) ? min.x : max.x;
			closeCorner.y = (center.y < min.y) ? min.y : max.y;
			// and test distance against radius
			TYPE dist = (center - closeCorner).magnitude();
			return dist < radius;
		}
		return false;
	}

	/** resizes the rectangle so it could be inscribed (circle inside) with a circle with the given radius */
	void setInscribedRadius(const TYPE a_radius) {
		TYPE w = getWidth(), h = getHeight();
		if (w == 0) setWidth((w = 1));
		if (h == 0) setWidth((h = 1));
		float constrainingDimension = ((w <= h) ? w : h);
		TYPE ratio = a_radius / constrainingDimension;
		V2<TYPE> center(getCenter());
		V2<TYPE> diffRad = diagonal() / 2;
		diffRad *= ratio;
		setMin(center - diffRad);
		setMax(center + diffRad);
	}

	V2<TYPE> getCenter() const { return V2<TYPE>::between(getMin(),getMax()); }

	void setCenter(const V2<TYPE> & c) {
		V2<TYPE> rad = diagonal() / 2;
		setMin(c - rad);
		setMax(c + rad);
	}

	/** @return a new rectangle that is beveled by the given delta x and delta y */
	AABB expandBorder(const V2<TYPE> & a_bevel) {
		return AABB<TYPE>(getMin() - a_bevel, getMax() + a_bevel);
	}

	bool intersectsAABB(const V2<TYPE> & a_min, const V2<TYPE> a_max) const {
		return!(a_max.x < min.x || a_max.y < min.y || a_min.x > max.x || a_min.y > max.y);
	}

	void clampToInt() { min.clampToInt(); max.clampToInt(); }

	bool intersectsAABB(const AABB<TYPE> & b) const {
		return!(b.max.x < min.x || b.max.y < min.y || b.min.x > max.x || b.min.y > max.y);
	}

	/** @return true if the given AABB is totally contained in this AABB */
	bool contains(const AABB<TYPE> & b) const {
		return min.x <= b.min.x && min.y <= b.min.y
			&& max.x >= b.max.x && max.y >= b.max.y;
	}	
	/** @return true if the given point is in this AABB */
	bool contains(const V2<TYPE> & p) const {
		return min.x <= p.x && min.y <= p.y
			&& max.x >= p.x && max.y >= p.y;
	}

	/**
	 * Does ray-casting against one polygon by detecting line collision
	 *
	 * @param rayStart
	 * @param rayDirection must be a unit vector
	 * @param out_rh information about the eventual hit (if any)
	 * @return true if raycast collision happened
	 */
	bool raycast(Ray_<TYPE> const & ray, RaycastHit_<TYPE> & out_rh) const {
		V2<TYPE> points[4];
		writePoints(points, 4);
		return V2<TYPE>::raycastPolygon(ray.start, ray.direction, points, 4, out_rh.distance, out_rh.point, out_rh.normal) != -1;
	}

	static const int BAD_VALUE = -1;
	/** used to determine X and Y. should not be used in setField() or getField() */
	static const int X = 0, Y = 1;
	/** MINX/MINY */
	static const int MIN = 0;	//V2::NUM_DIMENSIONS*0
	/** used by getField() and setField() */
	static const int MINX = MIN+X, MINY = MIN+Y;
	/** MAXX/MAXY */
	static const int MAX = 2;	//V2::NUM_DIMENSIONS*1
	/** used by getField() and setField() */
	static const int MAXX = MAX+X, MAXY = MAX+Y;
	/** X Position/Y Position */
	static const int POSITION = 4;	//V2::NUM_DIMENSIONS*2
	/** used by getField() and setField() */
	static const int POSITIONX = POSITION+X, POSITIONY = POSITION+Y;
	/** WIDTH/HEIGHT */
	static const int DIMENSION = 6;	//V2::NUM_DIMENSIONS*3
	/** used by getField() and setField() */
	static const int WIDTH = DIMENSION+X, HEIGHT = DIMENSION+Y;

	/** 
	 * @param point 
	 * @param out_normal will be set to the normal of the-point-at-the-edge (returned value)
	 * @return a point on the edge of this AABB that is as close as possible to point
	 */
	V2<TYPE> getClosestPointOnEdge(const V2<TYPE> point, V2<TYPE> & out_normal) const {
		V2<TYPE> p = point;
		out_normal.setZero();
		// check if the point is far out of range
		bool insideXrange = false, insideYrange = false;
		int justOutsideOfEdge = BAD_VALUE;
		if (p.x < min.x) { p.x = min.x; justOutsideOfEdge = MINX; }
		else if (p.x > max.x) { p.x = max.x; justOutsideOfEdge = MAXX; }
		else insideXrange = true;
		if (p.y < min.y) { p.y = min.y; justOutsideOfEdge = MINY; }
		else if (p.y > max.y) { p.y = max.y; justOutsideOfEdge = MAXY; }
		else insideYrange = true;
		// if is in neither the x or y range (it's closest to a corner)
		if (!insideXrange && !insideYrange) {
			out_normal = (point - p).normal();
			return p;
		}
		// if point is only one of the ranges, and outside another, no problem
		if ((insideXrange || insideYrange) && insideXrange ^ insideYrange) {
			out_normal = normalOfEdge(justOutsideOfEdge);
			return p;
		}
		// otherwise, it's inside the aabb. find out which range should be popped to the edge
		const int sides = V2<TYPE>::NUM_DIMENSIONS * 2;
		float distance, smallestDistance;
		int closestSide = BAD_VALUE;
		// check each side
		for (int i = 0; i < sides; ++i) {
			distance = abs(getField(i) - point.getField(i % V2<TYPE>::NUM_DIMENSIONS));
			// if the distance to this side is the smallest yet, that is the correct side to snap to
			if (closestSide == BAD_VALUE || distance < smallestDistance) {
				closestSide = i;
				smallestDistance = distance;
			}
		}
		out_normal = normalOfEdge(closestSide);
		p.setField(closestSide % V2<TYPE>::NUM_DIMENSIONS, getField(closestSide));
		return p;
	}

	/**
	* @param point
	* @param out_rh will be assigned to details about the closest raycast hit to that point
	* @return true if on linearedge, false if on corner
	*/
	bool getClosestRaycastHit(const V2<TYPE> point, RaycastHit_<TYPE> & out_rh) {
		V2<TYPE> p = point;
		out_rh.normal.setZero();
		// check if the point is far out of range
		bool insideXrange = false, insideYrange = false;
		int justOutsideOfEdge = BAD_VALUE;
		if (p.x < min.x) { p.x = min.x; justOutsideOfEdge = MINX; }
		else if (p.x > max.x) { p.x = max.x; justOutsideOfEdge = MAXX; }
		else insideXrange = true;
		if (p.y < min.y) { p.y = min.y; justOutsideOfEdge = MINY; }
		else if (p.y > max.y) { p.y = max.y; justOutsideOfEdge = MAXY; }
		else insideYrange = true;
		// if is in neither the x or y range (it's closest to a corner)
		if (!insideXrange && !insideYrange) {
			out_rh.point = p;
			out_rh.normal = (point - p);
			out_rh.distance = out_rh.normal.magnitude();
			out_rh.normal /= out_rh.distance;
			return false;
		}
		// if point is only one of the ranges, and outside another, no problem
		if (insideXrange ^ insideYrange) {
			out_rh.point = p;
			out_rh.normal = normalOfEdge(justOutsideOfEdge);
			out_rh.distance = (insideXrange ? abs(p.y - point.y) : abs(p.x - point.x));
			return true;
		}
		// otherwise, it's inside the aabb. find out which range should be popped to the edge
		const int sides = V2<TYPE>::NUM_DIMENSIONS * 2;
		float distance, smallestDistance;
		justOutsideOfEdge = BAD_VALUE;
		// check each side
		for (int i = 0; i < sides; ++i) {
			distance = abs(getField(i) - point.getField(i % V2<TYPE>::NUM_DIMENSIONS));
			// if the distance to this side is the smallest yet, that is the correct side to snap to
			if (justOutsideOfEdge == BAD_VALUE || distance < smallestDistance) {
				justOutsideOfEdge = i;
				smallestDistance = distance;
			}
		}
		p.setField(justOutsideOfEdge % V2<TYPE>::NUM_DIMENSIONS, getField(justOutsideOfEdge));
		out_rh.point = p;
		out_rh.distance = smallestDistance;
		out_rh.normal = normalOfEdge(justOutsideOfEdge);
		return true;
	}

	static V2<TYPE> normalOfEdge(const int edgeIndex) {
		switch (edgeIndex) {
		case MINX:	return V2<TYPE>(-1,0);
		case MINY:	return V2<TYPE>(0,-1);
		case MAXX:	return V2<TYPE>(1, 0);
		case MAXY:	return V2<TYPE>(0, 1);
		}
		return V2<TYPE>();
	}

	/**
	 * used to grab the opposite side of the range of the given dimension. 
	 * eg: MINX being passed in returns MAXX. MAXY passed in returns MINY.
	 * @return MINX<->MAXX, MINY<->MAXY
	 */
	static int getOppositeSide(int const a_side) {
		switch (a_side) {
		case MINX:	return MAXX;
		case MINY:	return MAXY;
		case MAXX:	return MINX;
		case MAXY:	return MINY;
		}
		return BAD_VALUE;
	}

	/** @param a_value {MINX, MINY, MAXX, MAXY} */
	TYPE getField(int const a_value) const {
		switch(a_value) {
		case MINX:	return getMinX();
		case MINY:	return getMinY();
		case MAXX:	return getMaxX();
		case MAXY:	return getMaxY();
		case POSITIONX:	return getX();
		case POSITIONY:	return getY();
		case WIDTH:		return getWidth();
		case HEIGHT:	return getHeight();
		}
		return 0;
	}
	/**
	 * a more data-driven way to call setMinX, setMaxX, setMinY, or setMaxY
	 * @param a_dimension {MINX,MAXX,MINY,MAXY}
	 */
	void setField(const int a_dimension, const TYPE a_value) {
		switch(a_dimension) {
		case MINX:	return setMinX(a_value);
		case MINY:	return setMinY(a_value);
		case MAXX:	return setMaxX(a_value);
		case MAXY:	return setMaxY(a_value);
		case POSITIONX:	return setX(a_value);
		case POSITIONY:	return setY(a_value);
		case WIDTH:		return setWidth(a_value);
		case HEIGHT:	return setHeight(a_value);
		}
	}

	/**
	 * @param a_range {MINX,MAXX,MINY,MAXY}
	 * @return true if this box and the given box share the passed in range 
	 * (x or y), meaing their width/height would overlap
	 */
	bool commonRange(const AABB<TYPE> & a_him, int a_range) const {
		switch(a_range) {
		case MINX:
		case MAXX:	return (getMaxY() >= a_him.getMinY() && getMinY() <= a_him.getMaxY());
		case MINY:
		case MAXY:	return (getMaxX() >= a_him.getMinX() && getMinX() <= a_him.getMaxX());
		}
		return false;
	}
	/**
	 * @param a_range {MINX,MAXX,MINY,MAXY}
	 * @return true if this box and the given box share the passed in range, excluding equality
	 * (x or y), meaing their width/height would overlap
	 */
	bool nonOrthogonalCommonRange(const AABB<TYPE> * a_him, int a_range) const {
		switch(a_range) {
		case MINX:
		case MAXX:	return (getMaxY() > a_him->getMinY() && getMinY() < a_him->getMaxY());
		case MINY:
		case MAXY:	return (getMaxX() > a_him->getMinX() && getMinX() < a_him->getMaxX());
		}
		return false;
	}

	/**
	 * @return true if this box and the given box share the passed in range, excluding equality
	 * (x or y), meaing their width/height would overlap
	 * this is different from intersects because intersects returns true if the bounds are adjacent
	 */
	bool nonOrthogonalIntersect(const AABB<TYPE> & b) const {
		return nonOrthogonalCommonRange(b, X) && nonOrthogonalCommonRange(b, Y);
	}

	/** @return true if this rectangle's a_myEdge is touching the given rectangles a_hisEdge */
	bool isEdgeTouchingEdge(const int & a_myEdge, const AABB<TYPE> * a_him,
		const int & a_hisEdge, TYPE a_touchRange = 1.0/1024) const
	{
		if(a_touchRange == 0.0)
			return (nonOrthogonalCommonRange(a_him, a_myEdge) 
			&& getField(a_myEdge) == a_him->getField(a_hisEdge));
		return (commonRange(a_him, a_myEdge)
			&& (abs(getField(a_myEdge) - a_him->getField(a_hisEdge)) < a_touchRange));
	}

	/**
	 * @return {MINX, MINY, MAXX, MAXY} if this Rectangle touches a_him. 
	 * BAD_VALUE (-1) if they don't touch.
	 */
	int getEdgeTouching(const AABB<TYPE> * a_him, TYPE const & a_neighborRange = 1.0/1024) {
		int myEdge, hisEdge;
		for(myEdge = 0; myEdge <= MAXY; ++myEdge) {
			hisEdge = (myEdge + MAX) % POSITION; // gets the opposite edge
			if(isEdgeTouchingEdge(myEdge, a_him, hisEdge, a_neighborRange))
				return myEdge;
		}
		return BAD_VALUE;
	}

	static const int NUM_DIMENSIONS = V2<TYPE>::NUM_DIMENSIONS;

#ifdef __GL_H__
	void glDraw(bool filled) const {
		glBegin(filled ? GL_POLYGON : GL_LINE_LOOP);
		V2<TYPE> Xy(getMaxX(), getMinY());
		V2<TYPE> xY(getMinX(), getMaxY());
		getMin().glVertex();
		Xy.glVertex();
		getMax().glVertex();
		xY.glVertex();
		glEnd();
	}
	void glDraw() const { glDraw(false); }
#endif
};

typedef AABB<float> RectF;
