#pragma once

#include "v2.h"

template<typename TYPE>
class _LineRef {
public:
	const V2<TYPE> * start;
	const V2<TYPE> * end;
	const V2<TYPE> & getStart() const { return *start; }
	const V2<TYPE> & getEnd() const { return *end; }
	void setStart(V2<TYPE> const & v) { *start = v; }
	void setEnd(V2<TYPE> const & v) { *end = v; }
	_LineRef(const V2<TYPE> * s, const V2<TYPE> * e) :start(s), end(e) {}
	_LineRef(){}
};

template<typename TYPE>
class _LinePoints {
public:
	V2<TYPE> start, end;
	const V2<TYPE> & getStart() const { return start; }
	const V2<TYPE> & getEnd() const { return end; }
	void setStart(V2<TYPE> const & v) { start = v; }
	void setEnd(V2<TYPE> const & v) { end = v; }
	_LinePoints(const V2<TYPE> * s, const V2<TYPE> * e) :start(*s), end(*e) {}
	_LinePoints(){}
};

/** for line calculations */
template<typename TYPE, typename REFTYPE>
class Line : public REFTYPE {
public:
	Line(){}
	Line(const V2<TYPE> * s, const V2<TYPE> * e) : REFTYPE(s, e){}

	V2<TYPE> getDelta() const { return getEnd() - getStart(); }

	TYPE length() const { return getDelta().magnitude(); }

	void translate(V2f const & moveDelta) { setStart(getStart() + moveDelta); setEnd(getEnd() + moveDelta); }

	/**
	 * @param point
	 * @param out_rh will be assigned to details about the closest point
	 * @return true if on the line
	 */
	bool getClosestRaycastHit(const V2<TYPE> point, RaycastHit_<TYPE> & out_rh) const {
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
		TYPE lineDistance, rayDistance;
		Ray_<TYPE> line(getStart(), getDelta());
		bool hit = line.intersection(ray, lineDistance, rayDistance);
		if (!hit || lineDistance < 0 || lineDistance >= 1 || rayDistance < 0) {
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

	bool rotate(const TYPE radians) {
		V2<TYPE> d = getDelta();
		d.rotate(radians);
		setEnd(getStart() + d);
		return true;
	}

	bool setRotation(const TYPE radians) {
		V2<TYPE> d(radians);
		d *= length();
		setEnd(getStart() + d);
		return true;
	}

	bool intersectsAABB(V2<TYPE> const & min, V2<TYPE> const & max) const {
		AABB<TYPE> r(min, max);
		RaycastHit_<TYPE> rh;
		V2<TYPE> d = getDelta();
		TYPE len = d.magnitude();
		if (len <= 0) return false;
		bool hit = r.raycast(Ray(start, d / len), rh);
		return (hit && rh.distance >= 0 && rh.distance < len);
	}

	bool intersectsCircle(V2<TYPE> const & center, const TYPE radius) const {
		V2<TYPE> p;
		return V2<TYPE>::lineCrossesCircle(getStart(), getEnd(), center, radius, p);
	}

};

typedef Line<float, _LineRef<float> > LinePF;
typedef Line<float, _LinePoints<float> > Linef;

template<typename TYPE>
class LineP : public Line<TYPE, _LineRef<TYPE> >{
public:
	LineP(){}
	LineP(const V2<TYPE> * s, const V2<TYPE> * e) : Line<TYPE, _LineRef<TYPE> >(s, e){}
};