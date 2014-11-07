#include <GL/freeglut.h>
#include "fsm_followagent.h"
#include "agent.h"
#include <stdio.h>
#include "fsm_fleebullets.h"
#include "bullet.h"
#include "game.h"

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
		float distance = V2f::distance(a->body.center, followed->body.center);
		if(distance > a->body.radius * 10 + followed->body.radius
		|| a->playerControlled 
		|| !followed->alive) {
			a->setFSM(new FSM_Idle());
			return;
		}
		a->acceleration = seek(followed->body.center, a);
		a->direction = (followed->body.center - a->body.center).normal();
	}
	else {
		a->acceleration = seek(a->targetPosition, a);
	}
	a->acceleration += obstacleAvoidance(&a->game->obstacles, &a->sensorArea, a, &calc) * 10;
}
void FSM_FollowAgent::draw(Agent * a, GLUTRenderingContext * g_screen) {
	for(int i = 0; i < calc.actualHits.size(); ++i) {
		calc.actualHits[i]->glDraw(true);
	}
	g_screen->setColor(0x0088ff);
	for (int i = 0; i < calc.hitLocations.size(); ++i) {
		V2f start, end;
		start = calc.hitLocations[i];
		end = start + calc.hitNormals[i] * calc.hitForce[i];
		g_screen->drawLine(start, end);
	}
}
