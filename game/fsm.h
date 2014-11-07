#pragma once

class Agent;
struct GLUTRenderingContext;

class FiniteStateMachine {
public:
	virtual void enter(Agent * a) = 0;
	virtual void execute(Agent * a, int a_ms) = 0;
	virtual void exit(Agent * a) = 0;
	virtual void draw(Agent * a, GLUTRenderingContext * g_screen){}
	virtual const char * getName() = 0;
	virtual ~FiniteStateMachine(){}
};
