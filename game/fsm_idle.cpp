#include "fsm_idle.h"
#include "agent.h"
#include <stdio.h>
#include "fsm_followagent.h"
#include "fsm_fleebullets.h"
#include "bullet.h"

void FSM_Idle::execute(Agent * a, int a_ms) {
	if (a->playerControlled) {
		a->acceleration = seek(a->targetPosition, a);
		return;
	}
	Bullet * b = a->findClosestBullet();
	if(b != NULL) {
		a->setFSM(new FSM_Flee(b));
		return;
	}
	Agent * other = a->findClosestPlayerControlledAgent();
	if(other != NULL) {
		a->setFSM(new FSM_FollowAgent(other));
		return;
	}
	a->acceleration = stop(a, a_ms);
}