#pragma once

#include "fsm.h"
#include "codegiraffe/templatevector.h"
#include "codegiraffe/v2.h"
#include "steering.h"

class Agent;
class Obstacle;

class FSM_FollowAgent : public FiniteStateMachine {
	Agent * followed;
	CalculationsFor_ObstacleAvoidance calc;
public:
	FSM_FollowAgent(Agent * a):followed(a){}
	void enter(Agent * a) { }
	void execute(Agent * a, int a_ms);
	void exit(Agent * a) { }
	void draw(Agent * a, GLUTRenderingContext * g_screen);
	const char * getName() { return "follow"; }
	// because FSM_FollowAgent allocates memory (in calc) it must have a destructor.
	~FSM_FollowAgent(){}
};