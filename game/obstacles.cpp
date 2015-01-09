#include "obstacles.h"
#include "agent.h"

bool ShapeBox::intersects(const Shape * o) const {
	const ShapeAABB * aabb = dynamic_cast<const ShapeAABB*>(o);		if (aabb != NULL) { return Boxf::intersectsAABB(aabb->min, aabb->max); }
	const ShapeCircle * circ = dynamic_cast<const ShapeCircle*>(o);	if (circ != NULL) { return Boxf::intersectsCircle(*circ); }
	const ShapeBox * boxs = dynamic_cast<const ShapeBox*>(o);		if (boxs != NULL) { return ((Boxf*)boxs)->intersects(*this); }
	const ShapePolygon * po = dynamic_cast<const ShapePolygon*>(o);	if (po != NULL) { return po->intersectsBox(*this); }
	const ShapeCone * cone = dynamic_cast<const ShapeCone *>(o);	if (cone != NULL) { return cone->intersectsBox(*this); }
	return false;
}

bool ShapeCone::intersects(const Shape * o) const {
	const ShapeAABB * aabb = dynamic_cast<const ShapeAABB*>(o);		if (aabb != NULL) { return Conef::intersectsAABB(aabb->min, aabb->max); }
	const ShapeCircle * circ = dynamic_cast<const ShapeCircle*>(o);	if (circ != NULL) { return Conef::intersectsCircle(circ->center, circ->radius); }
	const ShapeBox * boxs = dynamic_cast<const ShapeBox*>(o);		if (boxs != NULL) { return intersectsBox(*boxs); }
	const ShapePolygon * po = dynamic_cast<const ShapePolygon*>(o);	if (po != NULL) { return po->intersectsCone(*this); }
	const ShapeCone * cone = dynamic_cast<const ShapeCone *>(o);	if (cone != NULL) { return intersectsCone(*cone); }
	return false;
}

bool ShapePolygon::intersects(const Shape * o) const {
	const ShapeAABB * aabb = dynamic_cast<const ShapeAABB*>(o);		if (aabb != NULL) { return Polygon2f::intersectsAABB(aabb->min, aabb->max); }
	const ShapeCircle * circ = dynamic_cast<const ShapeCircle*>(o);	if (circ != NULL) { return Polygon2f::intersectsCircle(circ->center, circ->radius); }
	const ShapeBox * boxs = dynamic_cast<const ShapeBox*>(o);		if (boxs != NULL) { return Polygon2f::intersectsBox(*boxs); }
	const ShapePolygon * po = dynamic_cast<const ShapePolygon*>(o);	if (po != NULL) { return po->intersectsPolygon(*this); }
	const ShapeCone * cone = dynamic_cast<const ShapeCone *>(o);	if (cone != NULL) { return Polygon2f::intersectsCone(*cone); }
	return false;
}

bool ShapeLine::intersects(const Shape * o) const {
	RaycastHit rh;
	V2f d = getDelta();
	float len = d.magnitude();
	if (len <= 0) return false;
	bool hit = o->raycast(Ray(getStart(), d / len), rh);
	return (hit && rh.distance >= 0 && rh.distance < len);
}

bool ShapeLineP::intersects(const Shape * o) const {
	RaycastHit rh;
	V2f d = getDelta();
	float len = d.magnitude();
	if (len <= 0) return false;
	bool hit = o->raycast(Ray(getStart(), d / len), rh);
	return (hit && rh.distance >= 0 && rh.distance < len);
}