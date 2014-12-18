#include "agent.h"
#include "game.h"
#include "bullet.h"

Bullet * Agent::findClosestBullet() {
	Bullet * closestBullet = NULL;
	float closestBulletDistance;
	float thisBulletDistance;
	TemplateVector<Agent*> nearbyAgents;
	game->gatherListOfAgentsAt(CircF(body.center, body.radius * 5), nearbyAgents);
	// if there are any bullets nearby
	for (int i = 0; i < nearbyAgents.size(); ++i) {
		Bullet * b = dynamic_cast<Bullet*>(nearbyAgents[i]);
		if(b != NULL && (thisBulletDistance = V2f::distance(b->body.center, body.center)) 
											< (body.radius*5+b->body.radius)) {
			if(closestBullet == NULL || thisBulletDistance < closestBulletDistance) {
				closestBullet = b;
				closestBulletDistance = thisBulletDistance;
			}
		}
	}
	return closestBullet;
}
Agent * Agent::findClosestPlayerControlledAgent() {
	Agent * aggroed = NULL;
	// in a loop, check if any agents 
	for(int i = 0; i < game->agents.size(); ++i) {
		Agent * a = game->agents[i];
		if(a == this)
			continue;
			// that are player controlled
		if(a->playerControlled // are within aggroRadius distance 
			// (TODO make an aggro cone, TODO make aggroRadius variable)
		&& V2f::distance(a->body.center, body.center) < (body.radius*5+a->body.radius) ) {
			aggroed = a;
			break;
		}
	}
	return aggroed;
}

void Agent::update(int a_ms) {
	if(fsm != NULL) {
		//printf("%s ", fsm->getName());
		fsm->execute(this, a_ms);
	}
	updateMovement(a_ms);
}

/**
* @param subject what this Agent is trying to see
* @return true if no Obstacles interrupt this Agent's view of the given subject.
*/
bool Agent::hasLineOfSight(Obstacle * subject) {
	Obstacle * agents[] = { this, subject };
	const int numAgents = sizeof(agents) / sizeof(agents[0]);
	V2f delta = subject->getCenter() - body.center;
	float distance = delta.magnitude();
	V2f normal = delta / distance;
	Obstacle * obs;
	float dist;
	V2f hit, norm;
	return !game->raycast(body.center, normal, distance, true, obs, dist, hit, norm, agents, numAgents);
}
