#include <GL/freeglut.h>
#include "fsm_astar.h"
#include "agent.h"
#include <stdio.h>
#include "fsm_fleebullets.h"
#include "bullet.h"
#include "game.h"
#include "fsm_followagent.h"

void FSM_Astar::execute(Agent * a, int a_ms) {
	if (path->size() == 0 || a->hasLineOfSight(followed)) {
		a->setFSM(new FSM_FollowAgent(followed));
		return;
	}
	GraphNode * nextNode = path->getLast();
	a->acceleration = seek(nextNode->location, a);
	V2f delta = a->body.center - nextNode->location;
	float dist = delta.magnitude();
	V2f normal = delta / dist;
	Obstacle * obs;
	V2f hit, norm;
	Obstacle * ignoreMe = a;
	RaycastHit rh;
	if (dist < (a->body.radius * 4) && !a->game->raycast(Ray(a->body.center, normal), rh, dist, true, obs, 1)) {
		//delta.magnitude() < a->body.radius) {
		path->remove(path->size()-1);
	}
	a->acceleration += obstacleAvoidance(&a->game->obstacles, &a->sensorArea, a, 0) * 10;
}
void FSM_Astar::draw(Agent * a, GLUTRenderingContext * g_screen) {
	for (int i = 0; i < path->size(); ++i) {
		g_screen->setColor(0x008800);
		g_screen->drawCircle(path->get(i)->location, .2f, true);
		g_screen->setColor(0xffffff);
		g_screen->printf(path->get(i)->location, "%d", i);
	}
}
