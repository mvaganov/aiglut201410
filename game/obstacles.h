#pragma once

#include <GL/freeglut.h>
#include "codegiraffe/v2.h"
#include "codegiraffe/circle.h"
#include "codegiraffe/box.h"
#include "codegiraffe/cone.h"
#include "codegiraffe/polygon.h"

// pure virtual class -- an Interface
class Obstacle {
public:
	virtual bool intersects(const Obstacle * o) const = 0;

	virtual V2f getCenter() const = 0;
	virtual float getRadius() const = 0;
	virtual bool contains(V2f const & p) const = 0;
	virtual bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const = 0;
	virtual V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const = 0;
	virtual void glDraw(bool filled) const = 0;
	/** @return a data structure that identifies how this object should resolve collision with the given otherObject */
	virtual void * calculateCollisionResolution(Obstacle * otherObject) = 0;
	/** @param collisionData the result of calculateCollisionResolution. If memory is allocated, this method should de-allocate it */
	virtual void resolveCollision(Obstacle * otherObject, void * collisionData) = 0;
	virtual ~Obstacle(){};
};

class CircleObject : public CircF, public Obstacle {
public:
	CircleObject(CircF c):CircF(c){}
	V2f getCenter() const {return center;}
	float getRadius() const { return radius; }
	bool intersects(const Obstacle * o) const;
	// this method needs BoxObject to be defined before the method can be defined
	bool contains(V2f const & p) const { return CircF::contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
			return CircF::raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return CircF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { CircF::glDraw(filled); }
	void * calculateCollisionResolution(Obstacle * otherObject){ return 0; }
	void resolveCollision(Obstacle * o, void * collisionData){}
	~CircleObject(){}
};
class BoxObject : public BoxF, public Obstacle {
public:
	BoxObject(BoxF b):BoxF(b){}
	BoxObject(RectF r):BoxF(r.getCenter(), r.getDimension(), 0){}
	V2f getCenter() const {return center;}
	float getRadius() const { return size.magnitude()/2; }
	bool intersects(const Obstacle * o) const;
	bool contains(V2f const & p) const {return BoxF::contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
			return BoxF::raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return BoxF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { BoxF::glDraw(filled); }
	void * calculateCollisionResolution(Obstacle * otherObject){ return 0; }
	void resolveCollision(Obstacle * o, void * collisionData){}
	~BoxObject(){}
};
class ConeObject : public ConeF, public Obstacle {
public:
	ConeObject(ConeF c) :ConeF(c){}
	V2f getCenter() const { return ConeF::getCenter(); }
	float getRadius() const { return radius; }
	bool intersects(const Obstacle * o) const;
	bool contains(V2f const & p) const { return ConeF::contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
		return ConeF::raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return ConeF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { ConeF::glDraw(filled); }
	void * calculateCollisionResolution(Obstacle * otherObject){ return 0; }
	void resolveCollision(Obstacle * o, void * collisionData){}
	~ConeObject(){}
};
class PolygonObject : public Polygon2f, public Obstacle{
public:
	PolygonObject(Polygon2f c) :Polygon2f(c){}
	V2f getCenter() const { return Polygon2f::getCenter(); }
	float getRadius() const { return Polygon2f::getRadius(); }
	bool intersects(const Obstacle * o) const;
	bool contains(V2f const & p) const { return Polygon2f::contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
		return Polygon2f::raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return Polygon2f::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { Polygon2f::glDraw(filled); }
	void * calculateCollisionResolution(Obstacle * otherObject){ return 0; }
	void resolveCollision(Obstacle * o, void * collisionData){}
	~PolygonObject(){}
};
