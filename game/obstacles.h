#pragma once

#include <GL/freeglut.h>
#include "codegiraffe/v2.h"
#include "codegiraffe/circle.h"
#include "codegiraffe/box.h"
#include "codegiraffe/cone.h"
#include "codegiraffe/polygon.h"
#include "codegiraffe/line.h"
#include "codegiraffe/glutrenderingcontext.h"

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
	/** general-case collision detection */
	virtual bool raycast(Ray const & ray, RaycastHit & out_rh) const = 0;
	/** must return a point inside of the shape, presumably at it's geometric center */
	virtual V2f getCenter() const = 0;
	/** return the point about which this shape rotates */
	virtual V2f getOrigin() const = 0;
	/** @param moveDelta how to move this shape */
	virtual void translate(V2f const & moveDelta) = 0;
	/** how has this been rotated so far */
	virtual float getRotation() const = 0;
	/** return false if this shape cannot rotate */
	virtual bool rotate(const float radians) = 0;
	/** return false if this shape cannot rotate */
	virtual bool setRotation(const float radians) = 0;
	/** if the shape were generalized into a sphere, centered on getCenter, what is the minimum radius required to make sure all of the points of this shape are within the sphere */
	virtual float getRadius() const = 0;

	/** used for general-case collision, and spatially aware logic */
	virtual void getClosestRaycastHit(V2f const & point, RaycastHit & out_rh) const = 0;

	// TODO spherecast. including inverted, for things bouncing around inside

	/** draw using the given rendering context */
	virtual void draw(GLUTRenderingContext * g, bool filled) const = 0;
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
protected:
	long mask;
public:
	static const int NOTHING = 0;
	static const int EVERYTHING = -1;
	static const int VORONOI_NODE = 1 << 0;
	static const int DELAUNY_TRIANGULATION = 1 << 1;
	static const int VORONOI_EDGE = 1 << 2;
	static const int VORONOI_POLYHEDRON = 1 << 3;
	static const int STATIC = 1 << 4;
	static const int DYNAMIC = 1 << 5;
	static const int BLOCKS_LOS = 1 << 6;
	static const int MARKED_FOR_DELETE = 1 << 31;

	/** collide with everything by default (all the bits are set) */
	Obstacle() : mask(-1){}
	Obstacle(const long mask) :mask(mask){}
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & colledisionData){}
	void draw(GLUTRenderingContext * g, bool filled) const { getShape()->draw(g, filled); }

	/** allows a shape to ignore shapes that do not share bits in the mask */
	long getMask() const { return mask; }
	/** allows a shape to ignore shapes that do not share bits in the mask */
	void setMask(const long a_mask) { mask = a_mask; }
	/** binary-OR the given mask into this mask */
	void ensureMaskOverlaps(const int otherMask) { mask |= otherMask; }

	/**
	 * @param o the given shape to test collission against
	 * @param mask the layer to test. -1 for all masks
	 */
	bool intersects(const Shaped * o, const long mask) const {
		return masksCollide(mask) && getShape()->intersects(o->getShape());
	}
	/** AABB collision is common and important enough to deserve it's own access */
	bool intersectsAABB(V2f const & min, V2f const & max, const long mask) const {
		return masksCollide(mask) && getShape()->intersectsAABB(min, max);
	}
	/** Spherical collision is common and important enough to deserve it's own access */
	bool intersectsCircle(V2f const & center, const float radius, const long mask) const {
		return masksCollide(mask) && getShape()->intersectsCircle(center, radius);
	}
	/** check if the given point is inside the shape */
	bool contains(V2f const & p, const long mask) const {
		return masksCollide(mask) && getShape()->contains(p);
	}
	/** general-case collision detection */
	bool raycast(Ray const & ray, RaycastHit & out_rh, const long mask) const {
		return masksCollide(mask) && getShape()->raycast(ray, out_rh);
	}
};

class ShapePoint : public V2f, public Shape {
public:
	const char* getShapeName() const { return "Point"; }
	ShapePoint() {}
	ShapePoint(V2f & p) :V2f(p){}
	V2f getCenter() const { return *this; }
	V2f getOrigin() const { return *this; }
	virtual void translate(V2f const & moveDelta) { this->add(moveDelta); }
	bool rotate(const float radians) { return false; }
	float getRotation() const { return 0; }
	bool setRotation(const float radians) { return false; }
	float getRadius() const { return 0; }
	bool intersects(const Shape * o) const { return o->contains(*this); }
	bool intersectsAABB(V2f const & min, V2f const & max) const { return V2f::isBetweenAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return V2f::isWithin(radius, center); }
	bool contains(V2f const & p) const { return p.isEqual(*this); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { 
		if (ray.distanceFrom(*this) == 0) {
			out_rh.setFromLine(ray.start, *this);
			return V2f::dot(out_rh.normal, ray.direction) < 0;
		}
		return false;
	}
	void getClosestRaycastHit(V2f const & point, RaycastHit & out_rh) const { out_rh.setFromLine(point, *this); }
	void draw(GLUTRenderingContext * g, bool filled) const { glBegin(GL_POINTS); glVertex(); glEnd(); }
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & collisionData){}
	~ShapePoint(){}
};

class ShapeLine : public Linef, public Shape {
public:
	const char* getShapeName() const { return "Line"; }
	ShapeLine() {}
	ShapeLine(V2f * a, V2f * b) :Linef(a,b){}
	ShapeLine(Linef & line) :Linef(line){}
	V2f getCenter() const { return V2f::between(getStart(), getEnd()); }
	V2f getOrigin() const { return getStart(); }
	virtual void translate(V2f const & moveDelta) { Linef::translate(moveDelta); }
	float getRotation() const { return getDelta().normalizeIfNotZero().piRadians(); }
	bool rotate(const float radians) { return Linef::rotate(radians); }
	bool setRotation(const float radians) { return Linef::setRotation(radians); }
	float getRadius() const { return Linef::length() / 2; }
	bool intersects(const Shape * o) const;
	bool intersectsAABB(V2f const & min, V2f const & max) const { return Linef::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return Linef::intersectsCircle(center, radius); }
	bool contains(V2f const & p) const { return p.isCCW(getStart(), getEnd()) && p.distance(getCenter()) < getRadius(); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return Linef::raycast(ray, out_rh); }
	void getClosestRaycastHit(V2f const & point, RaycastHit & out_rh) const { Linef::getClosestRaycastHit(point, out_rh); }
	void draw(GLUTRenderingContext * g, bool filled) const { g->drawLine(getStart(), getEnd()); }
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & collisionData){}
	~ShapeLine(){}
}; 

class ShapeAABB : public AABBf, public Shape {
public:
	const char* getShapeName() const { return "AABB"; }
	ShapeAABB() {}
	ShapeAABB(AABBf r) :AABBf(r){}
	V2f getCenter() const { return AABBf::getCenter(); }
	V2f getOrigin() const { return min; }
	void translate(V2f const & moveDelta) { AABBf::translate(moveDelta); }
	float getRotation() const { return 0; }
	bool rotate(const float radians) { return false; }
	bool setRotation(const float radians) { return false; }
	float getRadius() const { return AABBf::getCircumscribedRadius(); }
	bool intersects(const Shape * o) const { return o->intersectsAABB(min, max); }
	bool intersectsAABB(V2f const & min, V2f const & max) const { return AABBf::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return AABBf::intersectsCircle(center, radius); }
	// this method needs BoxObject to be defined before the method can be defined
	bool contains(V2f const & p) const { return AABBf::contains(p); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return AABBf::raycast(ray, out_rh); }
	void getClosestRaycastHit(V2f const & point, RaycastHit & out_rh) const { AABBf::getClosestRaycastHit(point, out_rh); }
	void draw(GLUTRenderingContext * g, bool filled) const { g->drawRect(*this, filled); }
	void * calculateCollisionResolution(Collidable * otherObject){ return 0; }
	void resolveCollision(Collidable * o, void * & collisionData){}
	~ShapeAABB(){}
};

class ShapeCircle : public Circf, public Shape {
public:
	const char* getShapeName() const { return "Circle"; }
	ShapeCircle() {}
	ShapeCircle(Circf c) :Circf(c){}
	V2f getCenter() const {return center;}
	V2f getOrigin() const { return center; }
	void translate(V2f const & moveDelta) { center += moveDelta; }
	float getRotation() const { return 0; }
	bool rotate(const float radians) { return false; }
	bool setRotation(const float radians) { return false; }
	float getRadius() const { return radius; }
	bool intersects(const Shape * o) const { return o->intersectsCircle(center, radius); }
	bool intersectsAABB(V2f const & min, V2f const & max) const { return Circf::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return Circf::intersectsCircle(center, radius); }
	// this method needs BoxObject to be defined before the method can be defined
	bool contains(V2f const & p) const { return Circf::contains(p); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return Circf::raycast(ray, out_rh); }
	void getClosestRaycastHit(V2f const & point, RaycastHit & out_rh) const { Circf::getClosestRaycastHit(point, out_rh); }
	void draw(GLUTRenderingContext * g, bool filled) const { g->drawCircle(*this, filled); }
	~ShapeCircle(){}
};

class ShapeBox : public Boxf, public Shape {
public:
	const char* getShapeName() const { return "Box"; }
	ShapeBox() {}
	ShapeBox(Boxf b) :Boxf(b){}
	ShapeBox(AABBf r) :Boxf(r.getCenter(), r.getDimension(), 0){}
	float getRotation() const { return rotation.piRadians(); }
	V2f getCenter() const { return center; }
	V2f getOrigin() const { return center; }
	void translate(V2f const & moveDelta) { center += moveDelta; }
	bool rotate(const float radians) { rotation = V2f(rotation.piRadians() + radians); return true; }
	bool setRotation(const float radians) { rotation.rotate(radians); return true; }
	float getRadius() const { return size.magnitude() / 2; }
	bool intersects(const Shape * o) const;
	bool intersectsAABB(V2f const & min, V2f const & max) const { return Boxf::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return Boxf::intersectsCircle(center, radius); }
	bool contains(V2f const & p) const { return Boxf::contains(p); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return Boxf::raycast(ray, out_rh); }
	void getClosestRaycastHit(V2f const & point, RaycastHit & out_rh) const { Boxf::getClosestRaycastHit(point, out_rh); }
	void draw(GLUTRenderingContext * g, bool filled) const { g->drawBox(*this, filled); }
	~ShapeBox(){}
};
class ShapeCone : public Conef, public Shape {
public:
	const char* getShapeName() const { return "Cone"; }
	ShapeCone() {}
	ShapeCone(Conef c) :Conef(c){}
	V2f getCenter() const { return Conef::getCenter(); }
	float getRadius() const { return radius; }
	V2f getOrigin() const { return origin; }
	void translate(V2f const & moveDelta) { origin += moveDelta; }
	float getRotation() const { return Conef::getRotation(); }
	bool rotate(const float radians) { Conef::rotate(radians); return true; }
	bool setRotation(const float radians) { Conef::setRotation(radians); return true; }
	bool intersects(const Shape * o) const;
	bool intersectsAABB(V2f const & min, V2f const & max) const { return Conef::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return Conef::intersectsCircle(center, radius); }
	bool contains(V2f const & p) const { return Conef::contains(p); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return Conef::raycast(ray, out_rh); }
	void getClosestRaycastHit(V2f const & point, RaycastHit & out_rh) const { Conef::getClosestRaycastHit(point, out_rh); }
	void draw(GLUTRenderingContext * g, bool filled) const { Conef::glDraw(filled); }
	~ShapeCone(){}
};
class ShapePolygon : public Polygon2f, public Shape {
public:
	const char* getShapeName() const { return "Poly"; }
	ShapePolygon() {}
	ShapePolygon(Polygon2f c) :Polygon2f(c){}
	V2f getCenter() const { return Polygon2f::getCenter(); }
	float getRadius() const { return Polygon2f::getRadius(); }
	V2f getOrigin() const { return Polygon2f::getOrigin(); }
	void translate(V2f const & moveDelta) { Polygon2f::translate(moveDelta); }
	float getRotation() const { return Polygon2f::getRotation(); }
	bool rotate(const float radians) { Polygon2f::rotate(radians); return true; }
	bool setRotation(const float radians) { Polygon2f::setRotation(radians); return true; }
	bool intersects(const Shape * o) const;
	bool intersectsAABB(V2f const & min, V2f const & max) const { return Polygon2f::intersectsAABB(min, max); }
	bool intersectsCircle(V2f const & center, const float radius) const { return Polygon2f::intersectsCircle(center, radius); }
	bool contains(V2f const & p) const { return Polygon2f::contains(p); }
	bool raycast(Ray const & ray, RaycastHit & out_rh) const { return Polygon2f::raycast(ray, out_rh); }
	void getClosestRaycastHit(V2f const & point, RaycastHit & out_rh) const { Polygon2f::getClosestRaycastHit(point, out_rh); }
	void draw(GLUTRenderingContext * g, bool filled) const { Polygon2f::glDraw(filled); }
	~ShapePolygon(){}
};

#define OBSTACLE_OBJECT(NAME, SHAPETYPE, PRIMITIVESHAPE) \
class NAME : public Obstacle { \
	SHAPETYPE shape; public: Shape * getShape() { return &shape; } const Shape * getShape() const { return &shape; } \
	NAME(PRIMITIVESHAPE s, const long mask) : Obstacle(mask), shape(s) {}\
	~NAME(){} \
	SHAPETYPE * get##SHAPETYPE() { return &shape; } \
	const SHAPETYPE * get##SHAPETYPE() const { return &shape; } \
};

OBSTACLE_OBJECT(AABBObject, ShapeAABB, AABBf)
OBSTACLE_OBJECT(CircleObject, ShapeCircle, Circf)
OBSTACLE_OBJECT(BoxObject, ShapeBox, Boxf)
OBSTACLE_OBJECT(ConeObject, ShapeCone, Conef)
OBSTACLE_OBJECT(PolygonObject, ShapePolygon, Polygon2f)
OBSTACLE_OBJECT(LineObject, ShapeLine, Linef)
OBSTACLE_OBJECT(PointObject, ShapePoint, V2f)

#undef OBSTACLE_OBJECT