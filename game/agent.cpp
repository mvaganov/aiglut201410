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
	}
	updateMovement(a_ms);
}
