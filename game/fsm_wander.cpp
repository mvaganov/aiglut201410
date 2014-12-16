#include <GL/freeglut.h>
#include "fsm_wander.h"
#include "agent.h"
#include <stdio.h>
#include "fsm_fleebullets.h"
#include "bullet.h"
#include "game.h"
#include <math.h>

void FSM_Wander::execute(Agent * a, int a_ms) {
	
	wanderpoint.center=a->body.center+a->velocity;
	wanderpoint.radius=a->body.radius*2;
	validObstacles.clear();
	if (!a->alive || a->playerControlled) {
		a->setFSM(new FSM_Idle());
		return;
	}
	Agent * other = a->findClosestPlayerControlledAgent();
	if(other != NULL) {
		a->setFSM(new FSM_Idle());
		return;
	}
	for(int i = 0; i < a->game->obstacles.size(); ++i) {
		// ignore agents as valid obstacles
		if(dynamic_cast<Agent*>(a->game->obstacles[i]) == 0) {
			float dist = V2f::distance(a->body.center, a->game->obstacles[i]->getCenter());
			if(dist < 5) {
				validObstacles.add(a->game->obstacles[i]);
			}
		}
	}
	if(a->wanderTime>=a->wanderEnd){
		float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		float theta = (2.0 * V_PI) * r;
		wandertarget.x = cos(theta);
		wandertarget.y = sin(theta);
		wandertarget=wandertarget+wanderpoint.center;
		a->wanderTime=0;
	}
	a->wanderTime++;
	a->acceleration = seek(wandertarget, a);
	//a->acceleration += obstacleAvoidance(&validObstacles, &a->sensorArea, a, 0) * 20;
	a->direction = (wandertarget - a->body.center).normal();
}

void FSM_Wander::draw(Agent * a, GLUTRenderingContext * g_screen) {
	for(int i = 0; i < validObstacles.size(); ++i) {
		validObstacles[i]->glDraw(true);
	}
	g_screen->setColor(0x0088ff);
	g_screen->drawCircle(wanderpoint.center, wanderpoint.radius, false);
	g_screen->drawCircle(wandertarget,0.5f,false);
		//g_screen->printf(hideLocations[i], "%.2f", distance[i]);
}
