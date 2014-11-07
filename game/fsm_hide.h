#pragma once

#include "fsm.h"
#include "codegiraffe/templatevector.h"
#include "codegiraffe/v2.h"
class Agent;
class Obstacle;

class FSM_Hide : public FiniteStateMachine {
	Agent * predator;
	TemplateVector<Obstacle *> validObstacles;
	TemplateVector<V2f> hideLocations;
	TemplateVector<float> distance;
public:
	FSM_Hide(Agent * a):predator(a){}
	void enter(Agent * a) { }
	void execute(Agent * a, int a_ms);
	void exit(Agent * a) { }
	void draw(Agent * a, GLUTRenderingContext * g_screen);
	const char * getName() { return "hide"; }
	// because FSM_Hide allocates memory (TemplateVector member variables) it must have a destructor.
	~FSM_Hide(){}
};