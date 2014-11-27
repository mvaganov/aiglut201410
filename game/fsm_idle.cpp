#include "fsm_idle.h"
#include "agent.h"
#include <stdio.h>
#include "fsm_followagent.h"
#include "fsm_fleebullets.h"
#include "bullet.h"
#include "fsm_hide.h"
#include "game.h"
#include "fsm_alignment.h"
#include "fsm_cohesion.h"
#include "fsm_separation.h"

void FSM_Idle::execute(Agent * a, int a_ms) {
	if (a->playerControlled) {
		a->setFSM(new FSM_FollowAgent(NULL));
		return;
	}
	Bullet * b = a->findClosestBullet();
	if(b != NULL) {
		a->setFSM(new FSM_Flee(b));
	}
	Agent * other = a->findClosestPlayerControlledAgent();
	if(other != NULL) {
		a->setFSM(new FSM_FollowAgent(other));
		return;
	}
	// count other agents nearby, and if there are 3 or more, go to alignment
	TemplateVector<Agent*> nearby;
	a->game->gatherListOfAgentsAt(CircF(a->body.center, alignmentRange+a->body.radius), nearby);
	if(nearby.size() >= 3) {
		a->setFSM(new FSM_Alignment(alignmentRange, alignmentRange));
//		a->setFSM(new FSM_Cohesion(alignmentRange, alignmentRange));
//		a->setFSM(new FSM_Separation(alignmentRange, alignmentRange, 1));
		return;
	}
	a->acceleration = stop(a, a_ms);
}
