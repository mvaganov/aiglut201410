#pragma once

#include "v2.h"

/** for line calculations */
template<typename TYPE>
class Line {
public:
	const V2<TYPE> * start;
	const V2<TYPE> * end;

	Line(const V2<TYPE> * s, const V2<TYPE> * e) :start(s), end(e) {}
	Line(){}

	const V2<TYPE> & getStart() const { return *start; }
	const V2<TYPE> & getEnd() const { return *end; }
	V2<TYPE> getDelta() const { return getEnd() - getStart(); }

	/**
	 * @param point
	 * @param out_rh will be assigned to details about the closest point
	 * @return true if on the line
	 */
	bool getClosestRaycastHit(const V2<TYPE> point, RaycastHit_<TYPE> & out_rh) {
		TYPE lineDistance;
		Ray_<TYPE> line(getStart(), getDelta()), ray(point, line.direction.perp());
		bool hit = line.intersection(ray, lineDistance);
		if (!hit) {
			out_rh.set(point, V2<TYPE>::ZERO_DEGREES(), 0);
		} else if (lineDistance < 0) {
			out_rh.setFromLine(point, getStart());
			hit = false;
		} else if (lineDistance >= 1) {
			out_rh.setFromLine(point, getEnd());
			hit = false;
		} else {
			out_rh.point = getStart() + line.direction * lineDistance;
			out_rh.normal = ray.direction.normal();
			out_rh.distance = (point - out_rh.point).magnitude();
		}
		return hit;
	}

	bool raycast(Ray_<TYPE> const & ray, RaycastHit_<TYPE> & out_rh) const {
		TYPE lineDistance;
		Ray_<TYPE> line(getStart(), getDelta());
		bool hit = line.intersection(ray, lineDistance);
		if (!hit || lineDistance < 0 || lineDistance >= 1) {
			hit = false;
		} else {
			out_rh.point = getStart() + line.direction * lineDistance;
			out_rh.normal = line.direction.perp().normal();
			out_rh.distance = (out_rh.point - ray.start).magnitude();
		}
		return hit;
	}

	void gatherLineIntoRaycastHit(V2<TYPE> const & from, V2<TYPE> const & to, RaycastHit_<TYPE> & out_rh) const {
		out_rh.point = to;
		out_rh.normal = to - from;
		out_rh.distance = out_rh.normal.magnitude();
		out_rh.normal /= out_rh.distance;
	}
};

typedef Line<float> LineF;
