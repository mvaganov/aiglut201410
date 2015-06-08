#pragma once

#include "fsm.h"
#include "codegiraffe/templatevector.h"
#include "codegiraffe/v2.h"
#include "steering.h"

class GraphNode;
class Agent;
class Obstacle;

class FSM_Astar : public FiniteStateMachine {
	Agent * followed;
	TemplateVector<GraphNode*> * path;
public:
	FSM_Astar(Agent * a, TemplateVector<GraphNode*> * p):followed(a),path(p){}
	void enter(Agent * a) { }
	void execute(Agent * a, int a_ms);
	void exit(Agent * a) { }
	void draw(Agent * a, GLUTRenderingContext * g_screen);
	const char * getName() { return "A*"; }
	// because FSM_FollowAgent allocates memory (in calc) it must have a destructor.
	~FSM_Astar(){delete path;path=0;}
};