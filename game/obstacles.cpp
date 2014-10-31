#include "obstacles.h"
#include "agent.h"

// defines the intersects function after all the component classes have been defined
bool CircleObject::intersects(const Obstacle * o) const {
	const CircleObject * co = dynamic_cast<const CircleObject*>(o);
	if(co != NULL) {
		return CircF::intersects(*co);
	}
	const BoxObject * bo = dynamic_cast<const BoxObject*>(o);
	if(bo != NULL) {
		BoxF * bf = (BoxF*)bo;
		return bf->intersects(*this);
	}
	const Agent * a = dynamic_cast<const Agent*>(o);
	if(a != NULL) {
		return CircF::intersects(a->body);
	}
	return false;
}

bool BoxObject::intersects(const Obstacle * o) const {
	const CircleObject * co = dynamic_cast<const CircleObject*>(o);
	if(co != NULL) {
		return BoxF::intersects(*co);
	}
	const BoxObject * bo = dynamic_cast<const BoxObject*>(o);
	if(bo != NULL) {
		BoxF * bf = (BoxF*)bo;
		return bf->intersects(*this);
	}
	const Agent * a = dynamic_cast<const Agent*>(o);
	if(a != NULL) {
		return BoxF::intersects(a->body);
	}
	return false;
}
