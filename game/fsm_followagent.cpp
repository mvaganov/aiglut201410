#include "fsm_followagent.h"
#include "agent.h"
#include <stdio.h>

void FSM_FollowAgent::execute(Agent * a, int a_ms) {
	Bullet * b = a->findClosestBullet();
	if(b != NULL) {
		//a->setFSM(new FSM_FleeBullet(b));
		return;
	}
	float distance = V2f::distance(a->body.center, followed->body.center);
	if(distance > a->body.radius * 5 + followed->body.radius || !followed->playerControlled) {
		a->setFSM(new FSM_Idle());
		return;
	}
	a->acceleration = seek(followed->body.center, a);
}
