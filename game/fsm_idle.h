#pragma once

#include "fsm.h"

class FSM_Idle : public FiniteStateMachine
{
	float alignmentRange;
public:
	FSM_Idle():alignmentRange(5){}
	void enter(Agent * a) { }
	void execute(Agent * a, int a_ms);
	void exit(Agent * a) { }
	const char * getName() { return "idle"; }
};
