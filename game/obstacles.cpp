#include "obstacles.h"
#include "agent.h"

bool BoxObject::intersects(const Shape * o) const {
	const AABBObject * aabb = dynamic_cast<const AABBObject*>(o);
	if (aabb != NULL) {
		return BoxF::intersectsAABB(aabb->min, aabb->max);
	}
	const CircleObject * co = dynamic_cast<const CircleObject*>(o);
	if(co != NULL) {
		return BoxF::intersectsCircle(*co);
	}
	const Agent * a = dynamic_cast<const Agent*>(o);
	if (a != NULL) {
		return BoxF::intersectsCircle(a->body);
	}
	const BoxObject * bo = dynamic_cast<const BoxObject*>(o);
	if(bo != NULL) {
		BoxF * bf = (BoxF*)bo;
		return bf->intersects(*this);
	}
	const PolygonObject * po = dynamic_cast<const PolygonObject*>(o);
	if (po != NULL) {
		return po->intersectsBox(*this);
	}
	const ConeObject * con = dynamic_cast<const ConeObject *>(o);
	if (con != NULL) {
		return con->intersectsBox(*this);
	}
	return false;
}

bool ConeObject::intersects(const Shape * o) const {
	const AABBObject * aabb = dynamic_cast<const AABBObject*>(o);
	if (aabb != NULL) {
		return ConeF::intersectsAABB(aabb->min, aabb->max);
	}
	const CircleObject * co = dynamic_cast<const CircleObject*>(o);
	if (co != NULL) {
		return ConeF::intersectsCircle(co->center, co->radius);
	}
	const Agent * a = dynamic_cast<const Agent*>(o);
	if (a != NULL) {
		return ConeF::intersectsCircle(a->body.center, a->body.radius);
	}
	const BoxObject * bo = dynamic_cast<const BoxObject*>(o);
	if (bo != NULL) {
		return intersectsBox(*bo);
	}
	const PolygonObject * po = dynamic_cast<const PolygonObject*>(o);
	if (po != NULL) {
		return po->intersectsCone(*this);
	}
	const ConeObject * con = dynamic_cast<const ConeObject *>(o);
	if (con != NULL) {
		return intersectsCone(*con);
	}
	return false;
}

bool PolygonObject::intersects(const Shape * o) const {
	const AABBObject * aabb = dynamic_cast<const AABBObject*>(o);
	if (aabb != NULL) {
		return Polygon2f::intersectsAABB(aabb->min, aabb->max);
	}
	const CircleObject * co = dynamic_cast<const CircleObject*>(o);
	if (co != NULL) {
		return Polygon2f::intersectsCircle(co->center, co->radius);
	}
	const Agent * a = dynamic_cast<const Agent*>(o);
	if (a != NULL) {
		return Polygon2f::intersectsCircle(a->body.center, a->body.radius);
	}
	const BoxObject * bo = dynamic_cast<const BoxObject*>(o);
	if (bo != NULL) {
		return Polygon2f::intersectsBox(*bo);
	}
	const PolygonObject * po = dynamic_cast<const PolygonObject*>(o);
	if (po != NULL) {
		return po->intersectsPolygon(*this);
	}
	const ConeObject * con = dynamic_cast<const ConeObject *>(o);
	if (con != NULL) {
		return Polygon2f::intersectsCone(*con);
	}
	return false;
}
