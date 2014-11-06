#include "fsm_fleebullets.h"
#include "agent.h"
#include <stdio.h>
#include "bullet.h"

void FSM_Flee::execute(Agent * a, int a_ms) {
	Bullet * b = a->findClosestBullet();
	if (b != NULL) {
		threat = b;
	}
	float distance = V2f::distance(a->body.center, threat->body.center);
	if (distance > a->body.radius * 5 + threat->body.radius || !threat->alive) {
		a->setFSM(new FSM_Idle());
		return;
	}
	a->acceleration = flee(threat->body.center, a);
}
