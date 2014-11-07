#include "agent.h"
#include "steering.h"
#include "game.h"

//vec2 seek(vec2 target, Agent agent) {
V2f seek(V2f target, Agent * agent) {
//	vec2 desired = target - agent.position;
	V2f desired = target - agent->body.center;
//	desired *= agent.maxSpeed / desired.mag();
	//desired *= agent->maximumSpeed / desired.magnitude();
	desired.normalize();
	desired *= agent->maximumSpeed;
//	vec2 force = desired - agent.velocity;
	V2f force = desired - agent->velocity;
//	return force * (agent.maxForce / agent.maxSpeed);
	return force * (agent->maximumForce / agent->maximumSpeed);
//}
}

V2f stop(Agent * agent, int a_ms) {
	if(!agent->velocity.isZero()) {
		V2f v = agent->velocity;
		if (a_ms != 0) {
			V2f perfectAccelToStop = -v * 1000.0f / (float)a_ms;
			float mag = perfectAccelToStop.magnitude();
			if(mag <= agent->maximumForce)
				return perfectAccelToStop;
		}
		return -v * (agent->maximumForce / v.magnitude());
	}
	return V2f::ZERO();
}

V2f flee(V2f target, Agent * agent) {
	V2f fleeAccel = -seek(target, agent);
//	if(fleeAccel.magnitude() > agent->maximumSpeed)
	return fleeAccel;
}

V2f obstacleAvoidance(TemplateVector<Obstacle*> * obstacles, Obstacle * sensorArea, Agent * a, CalculationsFor_ObstacleAvoidance * calc) {
	if (calc) calc->clear();
	V2f totalForce;
	for (int i = 0; i < obstacles->size(); ++i) {
		Obstacle * actuallyHit = obstacles->get(i);
		if (actuallyHit != a && actuallyHit->intersects(sensorArea)) {
			float totalDistance = a->velocity.magnitude();
			V2f normal;
			V2f point = actuallyHit->getClosestPointOnEdge(a->body.center, normal);
			V2f closestOnVelocity;
			V2f::closestPointOnLine(
				a->body.center, a->body.center + a->velocity,
				point, closestOnVelocity);
			V2f delta = closestOnVelocity - point;
			delta.normalize();
			float distFromStart = (point - a->body.center).magnitude();
			float forceForThisHit = totalDistance - distFromStart;
			if (calc) {
				calc->actualHits.add(actuallyHit);
				calc->hitLocations.add(point);
				calc->hitNormals.add(delta);
				calc->hitForce.add(forceForThisHit);
			}
			totalForce += delta * forceForThisHit;
		}
	}
	return totalForce;
}
