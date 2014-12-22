#pragma once

#include <GL/freeglut.h>
#include "codegiraffe/v2.h"
#include "codegiraffe/circle.h"
#include "codegiraffe/box.h"
#include "codegiraffe/cone.h"
#include "codegiraffe/polygon.h"

class Shape;

/** if something has a shape, but should not be considered a shape (composition vs inheritance relationship), then it should inherit this class */
class Shaped {
public:
	/** how is this thing shaped again? */
	virtual Shape * getShape() = 0;
	/** for const correctness */
	virtual const Shape * getShape() const = 0;
	/** interfaces without virtual destructors allow potential memory leaks. */
	virtual ~Shaped(){}
};

/** if it is a shape, then all of this must also be knowable. */
class Shape : public Shaped {
public:
	/** name of the shape */
	virtual const char* getShapeName() const = 0;
	/** general-case collision detection with other shapes */
	virtual bool intersects(const Shape * o) const = 0;
	/** AABB collision is common and important enough to deserve it's own access */
	virtual bool intersectsAABB(V2f const & min, V2f const & max) const = 0;
	/** Spherical collision is common and important enough to deserve it's own access */
	virtual bool intersectsCircle(V2f const & center, const float radius) const = 0;
	/** check if the given point is inside the shape */
	virtual bool contains(V2f const & p) const = 0;
	/** must return a point inside of the shape */
	virtual V2f getCenter() const = 0;
	/** if the shape were generalized into a sphere, centered on getCenter, what is the minimum radius required to make sure all of the points of this shape are within the sphere */
	virtual float getRadius() const = 0;
	/** general-case collision detection */
	virtual bool raycast(Ray const & ray, RaycastHit & out_rh) const = 0;
	/** used for general-case collision, and spatially aware logic */
	virtual V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const = 0;
	/** TODO rename to be more general, use a rendering context rather than expecting the context is globally accessible */
	virtual void glDraw(bool filled) const = 0;
	virtual Shape * getShape() { return this; }
	virtual const Shape * getShape() const { return this; }
	virtual ~Shape(){};
};

/** anything that can be collided with has a shape, but is not neccessarily a shape. */
class Collidable : public Shaped {
public:
	/** @return a data structure that identifies how this object should resolve collision with the given otherObject. NULL for no collision. */
	virtual void * calculateCollisionResolution(Collidable * otherObject) = 0;
	/** @param collisionData the result of calculateCollisionResolution. If memory is allocated, this method should de-allocate it */
	virtual void resolveCollision(Collidable * otherObject, void * & collisionData) = 0;
	virtual ~Collidable(){}

	/** allows a shape to ignore shapes that do not share bits in the mask */
	virtual long getMask() const = 0;
	/** allows a shape to ignore shapes that do not share bits in the mask */
	virtual void setMask(const long a_mask) = 0;
	/** should be called before collision checks, to make sure collision is even possible between the given shapes */
	bool masksCollide(const Collidable * s) const { return (getMask() & s->getMask()) != 0; }
	bool masksCollide(const long otherMask) const { return (getMask() & otherMask) != 0; }
};

/** an obstacle is Collidable, though a basic obstacle doesn't move or change based on collision */
class Obstacle : public Collidable {
	long mask;
public:
	static const int EVERYTHING = -1;
	static const int VORONOI_NODE = 1 << 0;
	static const int DELAUNY_TRIANGULATION = 1 << 1;
	static const int VORONOI_EDGE = 1 << 2;
	static const int VORONOI_POLYHEDRON = 1 << 3;

	/** collide with everything by default (all the bits are set) */
	Obstacle() : mask(-1){}
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & colledisionData){}
	void glDraw(bool filled) const { getShape()->glDraw(filled); }
	bool intersects(const Shaped * o) const {
		const Collidable * c = dynamic_cast<const Collidable*>(o);
		if (c == NULL || masksCollide(c))
			return getShape()->intersects(o->getShape());
		return false;
	}
	/** allows a shape to ignore shapes that do not share bits in the mask */
	long getMask() const { return mask; }
	/** allows a shape to ignore shapes that do not share bits in the mask */
	void setMask(const long a_mask) { mask = a_mask; }
};

class ShapeAABB : public RectF, public Shape {
public:
	const char* getShapeName() const { return "AABB"; }
	ShapeAABB() {}
	ShapeAABB(RectF r) :RectF(r){}
	V2f getCenter() const { return RectF::getCenter(); }
	float getRadius() const { return RectF::getCircumscribedRadius(); }
	bool intersects(const Shape * o) const { return o->intersectsAABB(min, max); }
	bool intersectsAABB(V2f const & min, V2f const & max) const { return RectF::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return RectF::intersectsCircle(center, radius); }
	// this method needs BoxObject to be defined before the method can be defined
	bool contains(V2f const & p) const { return RectF::contains(p); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return RectF::raycast(ray, out_rh); }
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
	ShapeCircle() {}
	ShapeCircle(CircF c) :CircF(c){}
	V2f getCenter() const {return center;}
	float getRadius() const { return radius; }
	bool intersects(const Shape * o) const { return o->intersectsCircle(center, radius); }
	bool intersectsAABB(V2f const & min, V2f const & max) const { return CircF::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return CircF::intersectsCircle(center, radius); }
	// this method needs BoxObject to be defined before the method can be defined
	bool contains(V2f const & p) const { return CircF::contains(p); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return CircF::raycast(ray, out_rh); }
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return CircF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { CircF::glDraw(filled); }
	~ShapeCircle(){}
};

class ShapeBox : public BoxF, public Shape {
public:
	const char* getShapeName() const { return "Box"; }
	ShapeBox() {}
	ShapeBox(BoxF b) :BoxF(b){}
	ShapeBox(RectF r) :BoxF(r.getCenter(), r.getDimension(), 0){}
	V2f getCenter() const {return center;}
	float getRadius() const { return size.magnitude()/2; }
	bool intersects(const Shape * o) const;
	bool intersectsAABB(V2f const & min, V2f const & max) const { return BoxF::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return BoxF::intersectsCircle(center, radius); }
	bool contains(V2f const & p) const { return BoxF::contains(p); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return BoxF::raycast(ray, out_rh); }
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return BoxF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { BoxF::glDraw(filled); }
	~ShapeBox(){}
};
class ShapeCone : public ConeF, public Shape {
public:
	const char* getShapeName() const { return "Cone"; }
	ShapeCone() {}
	ShapeCone(ConeF c) :ConeF(c){}
	V2f getCenter() const { return ConeF::getCenter(); }
	float getRadius() const { return radius; }
	bool intersects(const Shape * o) const;
	bool intersectsAABB(V2f const & min, V2f const & max) const { return ConeF::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return ConeF::intersectsCircle(center, radius); }
	bool contains(V2f const & p) const { return ConeF::contains(p); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return ConeF::raycast(ray, out_rh); }
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return ConeF::getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { ConeF::glDraw(filled); }
	~ShapeCone(){}
};
class ShapePolygon : public Polygon2f, public Shape {
public:
	const char* getShapeName() const { return "Poly"; }
	ShapePolygon() {}
	ShapePolygon(Polygon2f c) :Polygon2f(c){}
	V2f getCenter() const { return Polygon2f::getCenter(); }
	float getRadius() const { return Polygon2f::getRadius(); }
	bool intersects(const Shape * o) const;
	bool intersectsAABB(V2f const & min, V2f const & max) const { return Polygon2f::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return Polygon2f::intersectsCircle(center, radius); }
	bool contains(V2f const & p) const { return Polygon2f::contains(p); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return Polygon2f::raycast(ray, out_rh); }
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
