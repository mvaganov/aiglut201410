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
return;
	TemplateSet<Obstacle*> possibleObstacles;
//	a->game->gatherStaticObstaclesAt(Circf(a->body.center, a->body.radius + 5), possibleObstacles);
	for (int i = 0; i < possibleObstacles.size(); ++i)
		validObstacles.add(possibleObstacles[i]);
	V2f predLoc = predator->body.center, obsLoc, delta;
	float closestDist = -1;
	V2f closestLoc;
	for(int i = 0; i < validObstacles.size(); ++i) {
		obsLoc = validObstacles[i]->getShape()->getCenter();
		delta = obsLoc - predLoc;
		float dist = delta.magnitude();
//		float rayDist;
//		V2f outEdge, outNormal;
		RaycastHit rh;
		validObstacles[i]->raycast(Ray(obsLoc, delta.normal()), rh, Obstacle::STATIC);// rayDist, outEdge, outNormal);
		V2f hideLoc = rh.point + delta.normal() * a->body.radius;
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
void FSM_Hide::draw(Agent * a, GLUTRenderingContext * g) {
	for(int i = 0; i < validObstacles.size(); ++i) {
		validObstacles[i]->getShape()->draw(g, true);
	}
	g->setColor(0x0088ff);
	for(int i = 0; i < hideLocations.size(); ++i) {
		g->drawCircle(hideLocations[i], a->body.radius, false);
		g->printf(hideLocations[i], "%.2f", distance[i]);
	}
}
