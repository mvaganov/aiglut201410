#include "obstacles.h"
#include "agent.h"

// defines the intersects function after all the component classes have been defined
bool CircleObject::intersects(const Obstacle * o) const {
	const CircleObject * co = dynamic_cast<const CircleObject*>(o);
	if(co != NULL) {
		return CircF::intersects(*co);
	}
	const Agent * a = dynamic_cast<const Agent*>(o);
	if (a != NULL) {
		return CircF::intersects(a->body);
	}
	const BoxObject * bo = dynamic_cast<const BoxObject*>(o);
	if(bo != NULL) {
		BoxF * bf = (BoxF*)bo;
		return bf->intersects(*this);
	}
	const ConeObject * con = dynamic_cast<const ConeObject *>(o);
	if (con != NULL) {
		return con->intersectCircle(center, radius);
	}
	return false;
}

bool BoxObject::intersects(const Obstacle * o) const {
	const CircleObject * co = dynamic_cast<const CircleObject*>(o);
	if(co != NULL) {
		return BoxF::intersects(*co);
	}
	const Agent * a = dynamic_cast<const Agent*>(o);
	if (a != NULL) {
		return BoxF::intersects(a->body);
	}
	const BoxObject * bo = dynamic_cast<const BoxObject*>(o);
	if(bo != NULL) {
		BoxF * bf = (BoxF*)bo;
		return bf->intersects(*this);
	}
	const ConeObject * con = dynamic_cast<const ConeObject *>(o);
	if (con != NULL) {
		return con->intersectsBox(*this);
	}
	return false;
}

bool ConeObject::intersects(const Obstacle * o) const {
	const CircleObject * co = dynamic_cast<const CircleObject*>(o);
	if (co != NULL) {
		return ConeF::intersectCircle(co->center, co->radius);
	}
	const Agent * a = dynamic_cast<const Agent*>(o);
	if (a != NULL) {
		return ConeF::intersectCircle(a->body.center, a->body.radius);
	}
	const BoxObject * bo = dynamic_cast<const BoxObject*>(o);
	if (bo != NULL) {
		return intersectsBox(*bo);
	}
	const ConeObject * con = dynamic_cast<const ConeObject *>(o);
	if (con != NULL) {
		return intersectsCone(*con);
	}
	return false;
}