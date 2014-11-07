#pragma once

#include "codegiraffe/v2.h"
#include "codegiraffe/templatevector.h"
#include "obstacles.h"
#include <stdio.h>

class Agent; // class protoype, a "forward declaration"

V2f seek(V2f target, Agent * agent);

V2f stop(Agent * agent, int a_ms);

V2f flee(V2f target, Agent * agent);

class CalculationsFor_ObstacleAvoidance {
public:
	TemplateVector<Obstacle *> actualHits;
	TemplateVector<V2f> hitLocations;
	TemplateVector<V2f> hitNormals;
	TemplateVector<float> hitForce;
	void clear() {
		actualHits.clear();
		hitLocations.clear();
		hitNormals.clear();
		hitForce.clear();
	}
};

V2f obstacleAvoidance(TemplateVector<Obstacle*> * obstacles,
	Obstacle * sensorArea, Agent * agent,
	CalculationsFor_ObstacleAvoidance * calc);