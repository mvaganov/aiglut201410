#include <GL/freeglut.h>
#include "fsm_evade.h"
#include "agent.h"
#include <stdio.h>
#include "fsm_fleebullets.h"
#include "bullet.h"
#include "game.h"

void FSM_Evade::execute(Agent * a, int a_ms) {
	validObstacles.clear();
	if (!target->alive || a->playerControlled) {
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
	//get targets actual location
	V2f targLoc = target->body.center, agentLoc, delta,futurePosition;
	//get the agents actual location
	agentLoc = a->getCenter();
	//find the targets future location within one second
	futurePosition=targLoc+target->velocity;
	//find agents's distance from targets future location
	delta=futurePosition-agentLoc;
	float dist=delta.magnitude();
	bool reach=false;
	float time=1;
	while(!reach){
	//if the agent can travel the distance within one second to targets future position exit
	if(target->maximumSpeed*time>=dist){
		reach=true;
		evade=futurePosition;
	}
	//if not increase targets position by one time interval and check if preditor can travel within that time
	futurePosition=futurePosition+target->velocity;
	delta=futurePosition-agentLoc;
	dist=delta.magnitude();
	time++;
	}
	a->acceleration = flee(evade, a);
	//a->acceleration += obstacleAvoidance(&validObstacles, &a->sensorArea, a, 0) * 20;
	a->direction = (evade - a->body.center).normal();
}

void FSM_Evade::draw(Agent * a, GLUTRenderingContext * g_screen) {
	for(int i = 0; i < validObstacles.size(); ++i) {
		validObstacles[i]->glDraw(true);
	}
	g_screen->setColor(0x0088ff);
	g_screen->drawCircle(evade, a->body.radius, false);
		//g_screen->printf(hideLocations[i], "%.2f", distance[i]);
}
