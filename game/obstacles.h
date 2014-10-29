#pragma once

#include "codegiraffe/v2.h"
#include "codegiraffe/circle.h"
#include "codegiraffe/box.h"

// pure virtual class -- an Interface
class Obstacle {
public:
	virtual bool intersects(const Obstacle * o) const = 0;

	virtual bool contains(V2f const & p) const = 0;
	virtual bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const = 0;
	virtual V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const = 0;
	virtual void glDraw(bool filled) const = 0;
};

class CircleObject : public Obstacle, public CircF {
public:
	CircleObject(CircF c):CircF(c){}
	bool intersects(const Obstacle * o) const {
		return false; // TODO
	}
	bool contains(V2f const & p) const { return CircF::contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
			return CircF::raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return CircF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { CircF::glDraw(filled); }
};
class BoxObject : public Obstacle, public BoxF {
public:
	BoxObject(BoxF b):BoxF(b){}
	BoxObject(RectF r):BoxF(r.getCenter(), r.getDimension(), 0){}
		bool intersects(const Obstacle * o) const {
		return false; // TODO
	}
	bool contains(V2f const & p) const {return BoxF::contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
			return BoxF::raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return BoxF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { BoxF::glDraw(filled); }
};
