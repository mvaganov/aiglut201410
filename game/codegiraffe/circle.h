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

	bool intersectsCircle(const Circle<TYPE> c) const {
		V2<TYPE> delta = center - c.center;
		return delta.magnitude() < radius + c.radius;
	}

	bool intersectsCircle(V2<TYPE> const & a_center, const TYPE a_radius) const {
		V2<TYPE> delta = center - a_center;
		return delta.magnitude() < radius + a_radius;
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

	bool intersectsAABB(const AABB<TYPE> r) const {
		return intersects(center, radius, r.getMin(), r.getMax());
	}

	bool intersectsAABB(V2<TYPE> const & min, V2<TYPE> const & max) const {
		return intersects(center, radius, min, max);
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
	 * @param out_rh information about the eventual hit (if any)
	 * @return true if raycast collision happened
	 */
	bool raycast(Ray_<TYPE> const & ray, RaycastHit_<TYPE> & out_rh) const {
		V2f p1, p2;
		if (raycast(ray, p1, p2)) {
			V2<TYPE> dp1 = (p1 - ray.start), dp2 = (p2 - ray.start);
			// check if the points are in the correct direction
			bool p1good = V2<TYPE>::dot(ray.direction, dp1) > 0;
			bool p2good = V2<TYPE>::dot(ray.direction, dp2) > 0;
			float dist1 = p1good ? dp1.magnitude() : -1, dist2 = p2good ? dp2.magnitude() : -1;
			// if both are good, invalidate the further one
			if (p1good && p2good) {
				if (dist1 < dist2) {
					p2good = false;
				} else {
					p1good = false;
				}
			}
			// if at least one is good...
			if (p1good || p2good) {
				// give the one that is good
				if (p1good) {
					out_rh.distance = dist1;
					out_rh.point = p1;
				} else {
					out_rh.distance = dist2;
					out_rh.point = p2;
				}
				out_rh.normal = (out_rh.point - center).normal();
				return true;
			}
		}
		return false;
	}

	/**
	 * @param ray
	 * @param out_p1 a point on the circle that the given ray interects
	 * @param out_p2 a point on the circle that the given ray interects, may not be different than out_p1
	 */
	bool raycast(Ray_<TYPE> const & ray, V2<TYPE> & out_p1, V2<TYPE> & out_p2) const {
		V2f rayNorm = ray.direction.normal();
		if (ray.start == center) {
			out_p1 = out_p2 = center + rayNorm * radius;
			return true;
		}
		TYPE radHit;
		Ray_<TYPE> radiusLine(center, rayNorm.perp());
		if (radiusLine.intersection(ray, radHit)) {
			if (abs(radHit) <= radius) {
				V2<TYPE> intersection = center + radiusLine.direction * radHit;
				// pythagorean theorum to get missing triangle side, where radius is the hypontinuse
				float side = sqrt(radius*radius - radHit*radHit);
				V2<TYPE> hitDelta = radiusLine.direction * side;
				// find which of the two solutions is closer to the ray start later
				out_p1 = intersection + hitDelta;
				out_p2 = intersection - hitDelta;
				return true;
			}
		}
		return false;
	}

	/**
	* @param point
	* @param out_normal will be set to the normal of the-point-at-the-edge (returned value)
	* @return a point on the edge of this Circle that is as close as possible to point
	*/
	V2<TYPE> getClosestPointOnEdge(const V2<TYPE> point, V2<TYPE> & out_normal) const {
		out_normal = point - center;
		if (out_normal.isZero()) out_normal = V2<TYPE>::ZERO_DEGREES();
		else out_normal.normalize();
		return center + (out_normal * radius);
	}

	/**
	 * @param point
	 * @param out_rh will be assigned to details about the closest raycast hit to that point
	 */
	void getClosestRaycastHit(const V2<TYPE> point, RaycastHit_<TYPE> & out_rh) {
		out_rh.normal = point - center;
		if (out_rh.normal.isZero()) {
			out_rh.set(point, V2<TYPE>::ZERO_DEGREES(), 0);
		} else {
			out_rh.distance = out_rh.normal.magnitude();
			out_rh.normal /= out_rh.distance;
			out_rh.point = center + (out_rh.normal * radius);
			out_rh.distance -= radius;
		}
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
