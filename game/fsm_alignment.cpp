#include <GL/freeglut.h>
#include "fsm_alignment.h"
#include "fsm_idle.h"
#include "agent.h"
#include "game.h"
#include "steering.h"
void FSM_Alignment::enter(Agent * a) {
	a->game->gatherListOfAgentsAt(CircF(a->body.center, flockRadius + a->body.radius), flock);
	flock.removeData(a);
}
void FSM_Alignment::execute(Agent * a, int a_ms) {
	if(a->playerControlled) {
		a->setFSM(new FSM_Idle());
		return;
	}
	V2f force = alignment(a, flock);
	if(force.isZero()) force = stop(a, a_ms);
	a->acceleration += force;
	a->acceleration += obstacleAvoidance(&a->game->obstacles, &a->sensorArea, a, 0) * 10;
}
void FSM_Alignment::draw(Agent * a, GLUTRenderingContext * g_screen) {
	g_screen->drawCircle(a->body.center, flockRadius + a->body.radius, false);
}
