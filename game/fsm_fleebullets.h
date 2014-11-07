#pragma once

#include "fsm.h"
class Agent;

class FSM_Flee : public FiniteStateMachine {
	Agent * threat;
public:
	FSM_Flee(Agent * a) : threat(a){}
	void enter(Agent * a) { }
	void execute(Agent * a, int a_ms);
	void exit(Agent * a) { }
	const char * getName() { return "flee"; }
};