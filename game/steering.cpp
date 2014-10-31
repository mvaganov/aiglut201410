#include "agent.h"
#include "steering.h"

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
