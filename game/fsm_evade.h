#pragma once

#include "fsm.h"
#include "codegiraffe/templatevector.h"
#include "codegiraffe/v2.h"
class Agent;
class Obstacle;

class FSM_Evade : public FiniteStateMachine {
	Agent * target;
	V2f evade;
	TemplateVector<Obstacle *> validObstacles;
public:
	FSM_Evade(Agent * a):target(a){}
	void enter(Agent * a) { }
	void execute(Agent * a, int a_ms);
	void exit(Agent * a) { }
	void draw(Agent * a, GLUTRenderingContext * g_screen);
	const char * getName() { return "evade"; }
	// because FSM_Hide allocates memory (TemplateVector member variables) it must have a destructor.
	~FSM_Evade(){}
};