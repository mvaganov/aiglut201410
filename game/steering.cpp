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

// TODO replace the obstacles list with the cellspace partition... seriosuly, do this.
V2f obstacleAvoidance(TemplateVector<Obstacle*> * obstacles, Obstacle * sensorArea, Agent * a, CalculationsFor_ObstacleAvoidance * calc) {
	if (calc) calc->clear();
	V2f totalForce;
	for (int i = 0; i < obstacles->size(); ++i) {
		Obstacle * actuallyHit = obstacles->get(i);
		if (actuallyHit != a && actuallyHit->intersects(sensorArea->getShape(), Obstacle::STATIC)) {
			float totalDistance = a->velocity.magnitude();
			V2f normal;
			RaycastHit closest;
			actuallyHit->getShape()->getClosestRaycastHit(a->body.center, closest);
			V2f closestOnVelocity;
			V2f::closestPointOnLine(
				a->body.center, a->body.center + a->velocity,
				closest.point, closestOnVelocity);
			V2f delta = closestOnVelocity - closest.point;
			delta.normalize();
			float distFromStart = (closest.point - a->body.center).magnitude();
			float forceForThisHit = totalDistance - distFromStart;
			if (calc) {
				calc->actualHits.add(actuallyHit);
				calc->hitLocations.add(closest.point);
				calc->hitNormals.add(delta);
				calc->hitForce.add(forceForThisHit);
			}
			totalForce += delta * forceForThisHit;
		}
	}
	return totalForce;
}

//vec2 Separation(Agent agent, List<Agent> neighbors) {
V2f separation(Agent * agent, float neighborRadius, TemplateVector<Agent*> & neighbors) {
	V2f totalPush;
	int countTooClose = 0;
	bool inRange = false;
	for(int i = 0; i < neighbors.size(); ++i) {
		V2f delta = agent->body.center - neighbors[i]->body.center;
		float dist = delta.magnitude() - (agent->body.radius + neighbors[i]->body.radius);
		if(dist < neighborRadius) {
			countTooClose++;
			totalPush += delta.normal();
		}
	}
	if(countTooClose > 1)
		return totalPush / (float)countTooClose;
	return V2f::ZERO();
}

//vec2 Cohesion(Agent agent, List<Agent> neighbors) {
V2f cohesion(Agent * agent, TemplateVector<Agent*> & neighbors) {
	//if (neighbors.empty()) return vec2.zero;
	if(neighbors.size() == 0) return V2f::ZERO();
	//vec2 centerOfMass = agent.position;
	V2f centerOfMass = agent->body.center;
	//for (Agent neighbor : neighbors)
	for(int i = 0; i < neighbors.size(); ++i)
		//centerOfMass += neighbor.position;
		centerOfMass += neighbors[i]->body.center;
	//centerOfMass /= neighbors.count();
	centerOfMass /= (float)neighbors.size();
	//vec2 desired = centerOfMass - agent.position;
	V2f desired = centerOfMass - agent->body.center;
	//desired *= agent.maxSpeed / desired.mag();
	desired *= agent->maximumSpeed / desired.magnitude();
	//vec2 force = desired - agent.velocity;
	V2f force = desired - agent->velocity;
	//return force * (agent.maxForce / agent.maxSpeed);
	return force * (agent->maximumForce / agent->maximumSpeed);
//}
}

//vec2 Alignment(Agent agent, List<Agent> neighbors) {
	//if (neighbors.empty()) return vec2.zero;
	//vec2 avgHeading = norm( agent.velocity );
	//for (Agent neighbor : neighbors)
		//avgHeading += norm( neighbor.velocity );
	//avgHeading /= neighbors.count();
	//vec2 desired = avgHeading * agent.maxSpeed;
	//vec2 force = desired - agent.velocity;
	//return force * (agent.maxForce / agent.maxSpeed);
//}

V2f alignment(Agent * agent, TemplateVector<Agent*> & neighbors) {
	V2f avgHeading = agent->velocity.normal();
	for(int i = 0; i < neighbors.size(); ++i) {
		avgHeading += neighbors[i]->velocity.normal();
	}
	avgHeading /= (float)neighbors.size();
	agent->direction = avgHeading.normal();
	return agent->direction;
}
