#include <GL/freeglut.h>
#include "fsm_hide.h"
#include "agent.h"
#include <stdio.h>
#include "fsm_fleebullets.h"
#include "bullet.h"
#include "game.h"

void FSM_Hide::execute(Agent * a, int a_ms) {
	validObstacles.clear();
	hideLocations.clear();
	distance.clear();
	if (!predator->alive || a->playerControlled) {
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
	V2f predLoc = predator->body.center, obsLoc, delta;
	float closestDist = -1;
	V2f closestLoc;
	for(int i = 0; i < validObstacles.size(); ++i) {
		obsLoc = validObstacles[i]->getCenter();
		delta = obsLoc - predLoc;
		float dist = delta.magnitude();
		float rayDist;
		V2f outEdge, outNormal;
		validObstacles[i]->raycast(obsLoc, delta.normal(), rayDist, outEdge, outNormal);
		V2f hideLoc = outEdge + delta.normal() * a->body.radius;
		hideLocations.add(hideLoc);
		float thisDist = V2f::distance(a->body.center, hideLoc);
		distance.add(thisDist);
		if(closestDist < 0 || thisDist < closestDist) {
			closestDist = thisDist;
			closestLoc = hideLoc;
		}
	}
	if(closestDist > 0) {
		if (closestDist < a->body.radius) {
			a->setFSM(new FSM_Idle());
			return;
		}
		a->acceleration = seek(closestLoc, a);
		a->acceleration += obstacleAvoidance(&validObstacles, &a->sensorArea, a, 0) * 20;
	}
	a->direction = (closestLoc - a->body.center).normal();
}
void FSM_Hide::draw(Agent * a, GLUTRenderingContext * g_screen) {
	for(int i = 0; i < validObstacles.size(); ++i) {
		validObstacles[i]->glDraw(true);
	}
	g_screen->setColor(0x0088ff);
	for(int i = 0; i < hideLocations.size(); ++i) {
		g_screen->drawCircle(hideLocations[i], a->body.radius, false);
		g_screen->printf(hideLocations[i], "%.2f", distance[i]);
	}
}
