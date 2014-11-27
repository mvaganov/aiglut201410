#pragma once

#include "fsm.h"
#include "codegiraffe/templatevector.h"
#include "codegiraffe/v2.h"

class Agent;

class FSM_Separation : public FiniteStateMachine {
	TemplateVector<Agent*> flock;
	float flockRadius;
	float tooFarAway;
	float comfortBubble;
public:
	FSM_Separation(float flockRadius, float tooFarAway, float comfortBubble)
		:flockRadius(flockRadius),tooFarAway(tooFarAway),comfortBubble(comfortBubble){}
	// look for fellow agents to separate with
	void enter(Agent * a);
	void execute(Agent * a, int a_ms);
	void exit(Agent * a) { }
	void draw(Agent * a, GLUTRenderingContext * g_screen);
	const char * getName() { return "separation"; }
	~FSM_Separation(){}
};
