#include "fsm_fleebullets.h"
#include "agent.h"
#include <stdio.h>
#include "bullet.h"
#include "game.h"
#include "fsm_hide.h"

void FSM_Flee::execute(Agent * a, int a_ms) {
	Bullet * b = a->findClosestBullet();
	if (b != NULL) {
		threat = b;
	}
	float distance = V2f::distance(a->body.center, threat->body.center);

	// if there are no more bullets, calm down.
	if (a->playerControlled || !threat->alive || distance > a->body.radius * 5 + threat->body.radius) {
		a->setFSM(new FSM_Idle());
		return;
	}
	// for a close call, hide!
	if (distance < a->body.radius * 2 + b->body.radius) {
		Agent * player = NULL;
		for (int i = 0; i < a->game->agents.size(); ++i) {
			if (a->game->agents[i]->playerControlled) {
				a->setFSM(new FSM_Hide(a->game->agents[i]));
				return;
			}
		}
	}
	a->acceleration = flee(threat->body.center, a);
	a->direction = a->acceleration.normal();
}
