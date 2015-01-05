#include <GL/freeglut.h>
#include "fsm_separation.h"
#include "fsm_idle.h"
#include "agent.h"
#include "game.h"
#include "steering.h"
void FSM_Separation::enter(Agent * a) {
	//a->game->gatherListOfAgentsAt(Circf(a->body.center, flockRadius + a->body.radius), flock);
	flock.removeData(a);
}
void FSM_Separation::execute(Agent * a, int a_ms) {
	if(a->playerControlled) {
		a->setFSM(new FSM_Idle());
		return;
	}
	V2f force = separation(a, comfortBubble, flock);
	if(force.isZero()) force = stop(a, a_ms);
	a->acceleration += force;
	a->acceleration += obstacleAvoidance(&a->game->obstacles, &a->sensorArea, a, 0) * 10;
}
void FSM_Separation::draw(Agent * a, GLUTRenderingContext * g_screen) {
	g_screen->drawCircle(a->body.center, flockRadius + a->body.radius, false);
}
