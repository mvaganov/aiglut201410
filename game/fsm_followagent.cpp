#include <GL/freeglut.h>
#include "fsm_followagent.h"
#include "agent.h"
#include <stdio.h>
#include "fsm_fleebullets.h"
#include "bullet.h"
#include "game.h"
#include "astar.h"

void FSM_FollowAgent::execute(Agent * a, int a_ms) {
	if(!a->playerControlled) {
		if(followed == NULL) {
			a->setFSM(new FSM_Idle());
			return;
		}

		Bullet * b = a->findClosestBullet();
		if(b != NULL) {
			a->setFSM(new FSM_Flee(b));
			return;
		}
		if (a->playerControlled || !followed->playerControlled || !followed->alive || !a->hasLineOfSight(followed)) {
			// TODO this would be a good place to search for a path using Astar....
			a->setFSM(new FSM_Idle());
			return;
		}
		a->acceleration = seek(followed->body.center, a);
		a->direction = (followed->body.center - a->body.center).normal();
	}
	else {
		a->acceleration = seek(a->targetPosition, a);
	}
	//a->acceleration += obstacleAvoidance(&a->game->obstacles, &a->sensorArea, a, &calc) * 10;
}
void FSM_FollowAgent::draw(Agent * a, GLUTRenderingContext * g) {
	for(int i = 0; i < calc.actualHits.size(); ++i) {
		calc.actualHits[i]->getShape()->draw(g, true);
	}
	if (!a->playerControlled && followed) {
		g->drawLine(followed->body.center, a->body.center);
	}
	g->setColor(0x0088ff);
	for (int i = 0; i < calc.hitLocations.size(); ++i) {
		V2f start, end;
		start = calc.hitLocations[i];
		end = start + calc.hitNormals[i] * calc.hitForce[i];
		g->drawLine(start, end);
	}
}
