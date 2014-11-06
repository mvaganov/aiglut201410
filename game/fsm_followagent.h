#pragma once

#include "fsm.h"
class Agent;

class FSM_FollowAgent : public FiniteStateMachine {
	Agent * followed;
public:
	FSM_FollowAgent(Agent * a):followed(a){}
	void enter(Agent * a) { }
	void execute(Agent * a, int a_ms);
	void exit(Agent * a) { }
};