#pragma once

#include <GL/freeglut.h>
#include "codegiraffe/v2.h"
#include "codegiraffe/circle.h"
#include "codegiraffe/box.h"
#include "codegiraffe/cone.h"
#include "codegiraffe/polygon.h"

// pure virtual class -- an Interface. TODO rename this class to be something that equates to a general-boundary, which may belong to a solid object or not, which may be in 2 or 3 dimensions. Shape?
class Shape {
public:
	virtual const char* getTypeName() const = 0;
	virtual bool intersects(const Shape * o) const = 0;
	virtual bool intersectsAABB(V2f const & min, V2f const & max) const = 0;
	virtual bool intersectsCircle(V2f const & center, const float radius) const = 0;

	virtual V2f getCenter() const = 0;
	virtual float getRadius() const = 0;
	virtual bool contains(V2f const & p) const = 0;
	virtual bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const = 0;
	virtual V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const = 0;
	virtual void glDraw(bool filled) const = 0;
	virtual ~Shape(){};
};

class Collidable {
public:
	/** @return a data structure that identifies how this object should resolve collision with the given otherObject. NULL for no collision. */
	virtual void * calculateCollisionResolution(Collidable * otherObject) = 0;
	/** @param collisionData the result of calculateCollisionResolution. If memory is allocated, this method should de-allocate it */
	virtual void resolveCollision(Collidable * otherObject, void * & collisionData) = 0;
	virtual ~Collidable(){}
};

class Obstacle : public Shape, public Collidable { };

class AABBObject : public RectF, public Obstacle {
public:
	const char* getTypeName() const { return "AABB"; }
	AABBObject(RectF r) :RectF(r){}
	V2f getCenter() const { return RectF::getCenter(); }
	float getRadius() const { return RectF::getCircumscribedRadius(); }
	bool intersects(const Shape * o) const { return o->intersectsAABB(min, max); }
	bool intersectsAABB(V2f const & min, V2f const & max) const { return RectF::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return RectF::intersectsCircle(center, radius); }
	// this method needs BoxObject to be defined before the method can be defined
	bool contains(V2f const & p) const { return RectF::contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
		return RectF::raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return RectF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { RectF::glDraw(filled); }
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & collisionData){}
	~AABBObject(){}
};

class CircleObject : public CircF, public Obstacle {
public:
	const char* getTypeName() const { return "Circle"; }
	CircleObject(CircF c) :CircF(c){}
	V2f getCenter() const {return center;}
	float getRadius() const { return radius; }
	bool intersects(const Shape * o) const { return o->intersectsCircle(center, radius); }
	bool intersectsAABB(V2f const & min, V2f const & max) const { return CircF::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return CircF::intersectsCircle(center, radius); }
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
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & collisionData){}
	~CircleObject(){}
};
class BoxObject : public BoxF, public Obstacle {
public:
	const char* getTypeName() const { return "Box"; }
	BoxObject(BoxF b) :BoxF(b){}
	BoxObject(RectF r):BoxF(r.getCenter(), r.getDimension(), 0){}
	V2f getCenter() const {return center;}
	float getRadius() const { return size.magnitude()/2; }
	bool intersects(const Shape * o) const;
	bool intersectsAABB(V2f const & min, V2f const & max) const { return BoxF::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return BoxF::intersectsCircle(center, radius); }
	bool contains(V2f const & p) const { return BoxF::contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
			return BoxF::raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return BoxF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { BoxF::glDraw(filled); }
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & collisionData){}
	~BoxObject(){}
};
class ConeObject : public ConeF, public Obstacle {
public:
	const char* getTypeName() const { return "Cone"; }
	ConeObject(ConeF c) :ConeF(c){}
	V2f getCenter() const { return ConeF::getCenter(); }
	float getRadius() const { return radius; }
	bool intersects(const Shape * o) const;
	bool intersectsAABB(V2f const & min, V2f const & max) const { return ConeF::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return ConeF::intersectsCircle(center, radius); }
	bool contains(V2f const & p) const { return ConeF::contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
		return ConeF::raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return ConeF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { ConeF::glDraw(filled); }
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & collisionData){}
	~ConeObject(){}
};
class PolygonObject : public Polygon2f, public Obstacle{
public:
	const char* getTypeName() const { return "Poly"; }
	PolygonObject(Polygon2f c) :Polygon2f(c){}
	V2f getCenter() const { return Polygon2f::getCenter(); }
	float getRadius() const { return Polygon2f::getRadius(); }
	bool intersects(const Shape * o) const;
	bool intersectsAABB(V2f const & min, V2f const & max) const { return Polygon2f::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return Polygon2f::intersectsCircle(center, radius); }
	bool contains(V2f const & p) const { return Polygon2f::contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
		return Polygon2f::raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return Polygon2f::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { Polygon2f::glDraw(filled); }
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & collisionData){}
	~PolygonObject(){}
};
