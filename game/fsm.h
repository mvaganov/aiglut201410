#pragma once

class Agent;

class FiniteStateMachine {
public:
	virtual void enter(Agent * a) = 0;
	virtual void execute(Agent * a, int a_ms) = 0;
	virtual void exit(Agent * a) = 0;
};
