#pragma once

#include "fsm.h"
#include "codegiraffe/templatevector.h"
#include "codegiraffe/v2.h"
#include "codegiraffe/circle.h"
class Agent;
class Obstacle;

class FSM_Wander : public FiniteStateMachine {
	CircF wanderpoint;
	//Agent * target;
	V2f wandertarget;
	TemplateVector<Obstacle *> validObstacles;
public:
	FSM_Wander(Agent * a){}
	void enter(Agent * a) { }
	void execute(Agent * a, int a_ms);
	void exit(Agent * a) { }
	void draw(Agent * a, GLUTRenderingContext * g_screen);
	const char * getName() { return "wander"; }
	// because FSM_Hide allocates memory (TemplateVector member variables) it must have a destructor.
	~FSM_Wander(){}
};