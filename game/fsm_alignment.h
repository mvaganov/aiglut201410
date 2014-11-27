#pragma once

#include "fsm.h"
#include "codegiraffe/templatevector.h"
#include "codegiraffe/v2.h"

class Agent;

class FSM_Alignment : public FiniteStateMachine {
	TemplateVector<Agent*> flock;
	float flockRadius;
	float tooFarAway;
public:
	FSM_Alignment(float flockRadius, float tooFarAway)
		:flockRadius(flockRadius),tooFarAway(tooFarAway){}
	// look for fellow agents to align with
	void enter(Agent * a);
	void execute(Agent * a, int a_ms);
	void exit(Agent * a) { }
	void draw(Agent * a, GLUTRenderingContext * g_screen);
	const char * getName() { return "alignment"; }
	~FSM_Alignment(){}
};
