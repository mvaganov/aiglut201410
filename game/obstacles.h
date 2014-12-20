#pragma once

#include <GL/freeglut.h>
#include "codegiraffe/v2.h"
#include "codegiraffe/circle.h"
#include "codegiraffe/box.h"
#include "codegiraffe/cone.h"
#include "codegiraffe/polygon.h"

/** pure virtual class -- an Interface. if it is a shape, then all of this is also knowable. */
class Shape {
public:
	virtual const char* getShapeName() const = 0;
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

/** if something has a shape, but should not be considered a shape (composition vs inheritance relationship), then it should inherit this class */
class Shaped {
public:
	virtual Shape * getShape() = 0;
	virtual const Shape * getShape() const = 0;
	virtual ~Shaped(){}
};

/** anything that can be collided with has a shape, but is not neccessarily a shape. */
class Collidable : public Shaped {
public:
	/** @return a data structure that identifies how this object should resolve collision with the given otherObject. NULL for no collision. */
	virtual void * calculateCollisionResolution(Collidable * otherObject) = 0;
	/** @param collisionData the result of calculateCollisionResolution. If memory is allocated, this method should de-allocate it */
	virtual void resolveCollision(Collidable * otherObject, void * & collisionData) = 0;
	virtual ~Collidable(){}
};
/** an obstacle is Collidable, though a basic obstacle doesn't move or change based on collision */
class Obstacle : public Collidable {
public:
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & colledisionData){}
	void glDraw(bool filled) const { getShape()->glDraw(filled); }
	bool intersects(const Shaped * o) const { return getShape()->intersects(o->getShape()); }
};

class ShapeAABB : public RectF, public Shape {
public:
	const char* getShapeName() const { return "AABB"; }
	ShapeAABB(RectF r) :RectF(r){}
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
	~ShapeAABB(){}
};

class ShapeCircle : public CircF, public Shape {
public:
	const char* getShapeName() const { return "Circle"; }
	ShapeCircle(CircF c) :CircF(c){}
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
	~ShapeCircle(){}
};

class ShapeBox : public BoxF, public Shape {
public:
	const char* getShapeName() const { return "Box"; }
	ShapeBox(BoxF b) :BoxF(b){}
	ShapeBox(RectF r) :BoxF(r.getCenter(), r.getDimension(), 0){}
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
	~ShapeBox(){}
};
class ShapeCone : public ConeF, public Shape {
public:
	const char* getShapeName() const { return "Cone"; }
	ShapeCone(ConeF c) :ConeF(c){}
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
	~ShapeCone(){}
};
class ShapePolygon : public Polygon2f, public Shape {
public:
	const char* getShapeName() const { return "Poly"; }
	ShapePolygon(Polygon2f c) :Polygon2f(c){}
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
	~ShapePolygon(){}
};

class AABBObject : public Obstacle {
	ShapeAABB shape; public: Shape * getShape() { return &shape; } const Shape * getShape() const { return &shape; }
	AABBObject(RectF r) : shape(r){}
	~AABBObject(){}
	ShapeAABB * getAABB() { return &shape; }
	const ShapeAABB * getAABB() const { return &shape; }
};
class CircleObject : public Obstacle {
	ShapeCircle shape; public: Shape * getShape() { return &shape; } const Shape * getShape() const { return &shape; }
	CircleObject(CircF c) : shape(c){}
	~CircleObject(){}
	ShapeCircle * getCircle() { return &shape; }
	const ShapeCircle * getCircle() const { return &shape; }
};
class BoxObject : public Obstacle {
	ShapeBox shape;	 public: Shape * getShape() { return &shape; } const Shape * getShape() const { return &shape; }
	BoxObject(BoxF b) :shape(b){}
	BoxObject(RectF r) :shape(r){}
	~BoxObject(){}
	ShapeBox * getBox() { return &shape; }
	const ShapeBox * getBox() const { return &shape; }
};
class ConeObject : public Obstacle {
	ShapeCone shape; public: Shape * getShape() { return &shape; } const Shape * getShape() const { return &shape; }
	ConeObject(ConeF c) : shape(c){}
	~ConeObject(){}
	ShapeCone * getCone() { return &shape; }
	const ShapeCone * getCone() const { return &shape; }
};
class PolygonObject : public Obstacle {
	ShapePolygon shape; public: Shape * getShape() { return &shape; } const Shape * getShape() const { return &shape; }
	PolygonObject(Polygon2f c) : shape(c){}
	~PolygonObject(){}
	ShapePolygon * getPolygon() { return &shape; }
	const ShapePolygon * getPolygon() const { return &shape; }
};
