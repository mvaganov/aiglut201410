#include "agent.h"
#include "game.h"
#include "bullet.h"

Bullet * Agent::findClosestBullet() {
	Bullet * closestBullet = NULL;
	float closestBulletDistance;
	float thisBulletDistance;
	// if there are any bullets nearby
	for(int i = 0; i < game->agents.size(); ++i) {
		Bullet * b = dynamic_cast<Bullet*>(game->agents[i]);
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
		fsm->execute(this, a_ms);
	} else {
		switch (behavior) {
		case BEHAVIOR_NONE:
			break;
		case BEHAVIOR_SEEK:
			acceleration = seek(targetPosition, this);
			break;
		case BEHAVIOR_AGGRO:
			V2f seekForce, fleeForce, stopForce;
			// search the game for nearby agents
			Agent * aggroed = findClosestPlayerControlledAgent();
			Bullet * closestBullet = findClosestBullet();
			// run from the closest bullet.
			if (closestBullet != NULL) {
				fleeForce = flee(closestBullet->body.center, this);
			}
			// if so, set target location to that player
			if (aggroed != NULL) {
				this->targetPosition = aggroed->body.center;
				// seek the agent
				seekForce = seek(targetPosition, this);
			}
			else { // otherwise
				stopForce = stop(this, a_ms); // stop moving
			}
			// fuzzy logic -- applying weights to multiple steering impulses
			acceleration = seekForce + fleeForce * 10 + stopForce * 0.5f;
			break;
		}
	}
	updateMovement(a_ms);
}
