#pragma once

#include "v2.h"
#include "aabb.h"

/**
 * Intended to be used with OpenGL, and more specifically, freeglut http://freeglut.sourceforge.net/docs/api.php
 * Written to be useful, readable, and informative for students studying 2D Vector math.
 *
 * @author mvaganov@hotmail.com October 2014
 */
template<typename TYPE>
struct Circle {
	V2<TYPE> center;
	TYPE radius;

	Circle():radius(0){}
	Circle(const V2<TYPE> center, const TYPE radius) :center(center), radius(radius){}
	Circle(const TYPE x, const TYPE y, const TYPE radius) :center(x,y), radius(radius){}

	void set(const V2<TYPE> center, const TYPE radius){ this->center.set(center); this->radius = radius; }
	void set(const TYPE x, const TYPE y, const TYPE radius){ center.set(x, y); this->radius = radius; }

	AABB<TYPE> getBounds() const {
		V2<TYPE> d(radius, radius);
		return AABB(center-d, center+d);
	}

	bool intersects(const Circle<TYPE> c) const {
		V2<TYPE> delta = center - c.center;
		TYPE minDist = radius + c.radius;
		return delta.magnitude() < minDist;
	}

	bool contains(V2<TYPE> const & p) const {
		V2<TYPE> delta = center - p;
		return delta.magnitude() < radius;
	}

	static bool intersects(const V2<TYPE> center, const TYPE radius, const V2<TYPE> & a_min, const V2<TYPE> & a_max) {
		// check if circle is inside the AABB
		bool xRange = center.x >= a_min.x && center.x <= a_max.x;
		bool yRange = center.y >= a_min.y && center.y <= a_max.y;
		if (xRange && yRange) return true; // inside the AABB
		// be ready to check if circle is within the extended range, which is radius units away from the AABB
		bool xRange2 = (center.x >= (a_min.x - radius) && center.x <= (a_max.x + radius)),
			yRange2 = (center.y >= (a_min.y - radius) && center.y <= (a_max.y + radius));
		if ((yRange && xRange2)// in extended range left or right			
		|| (xRange && yRange2)) {// in extended range above or below
			return true;
		}
		// if it is possibly in the extended range...
		if (xRange2 && yRange2) {
			// one of the 4 corners
			return(((center - a_max).magnitude() <= radius)
				|| ((center - a_min).magnitude() <= radius)
				|| ((center - V2<TYPE>(a_min.x, a_max.y)).magnitude() <= radius)
				|| ((center - V2<TYPE>(a_max.x, a_min.y)).magnitude() <= radius));
		}
		return false;
	}

	bool intersects(const AABB<TYPE> r) const {
		return intersects(center, radius, r.getMin(), r.getMax());
	}

	bool isEqual(const Circle<TYPE> & c) const {
		return center == c.center && radius == c.radius;
	}

	TYPE distanceBetween(const Circle<TYPE> & c) const {
		V2<TYPE> delta = center - c.center;
		return delta.magnitude() - c.radius - radius;
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
		float & out_dist, V2<TYPE> & out_point, V2<TYPE> & out_normal) {
		return V2<TYPE>::raycastCircle(
			rayStart, rayDirection, center, radius, out_dist, out_point, out_normal);
	}

	/**
	* @param point
	* @param out_normal will be set to the normal of the-point-at-the-edge (returned value)
	* @return a point on the edge of this Circle that is as close as possible to point
	*/
	V2<TYPE> getClosestPointOnEdge(const V2<TYPE> point, V2<TYPE> & out_normal) const {
		out_normal = point - center;
		out_normal.normalize();
		return center + (out_normal * radius);
	}

	bool operator==(const Circle<TYPE> & c) const { return isEqual(c); }

	TYPE getArea() const { return V_PI * radius * radius; }

#ifdef __GL_H__
	void glDraw(bool filled) const { glDrawCircle(center, radius, filled); }

	void glDraw() const { glDraw(false); }
#endif
};

typedef Circle<float> CircF;

#ifdef __GL_H__

bool glDrawCircle(const float a_x, const float a_y, const float a_radius, bool filled);

// inlined to ignore redefinition errors
inline bool glDrawCircle(const V2f & a_center, const float a_radius, bool filled) {
	return glDrawCircle(a_center.x, a_center.y, a_radius, filled);
}

#endif