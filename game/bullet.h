#pragma once
#include "agent.h"

class Bullet : public Agent
{
public:
	Bullet(CircF circle, Game * g, V2f direction, float speed)
		:Agent(circle, g)
	{
		setFSM(NULL);
		this->maximumSpeed = speed;
		this->behavior = BEHAVIOR_NONE;
		this->direction = direction;
		this->velocity = direction * speed;
		this->color = 0x0000ff;
	}

	void update(int a_ms) {
		Agent::update(a_ms);
		float distFromCenter = V2f::distance(body.center,V2f::ZERO());
		if(distFromCenter > 20) {
			this->alive = false;
		}
	}
	void * calculateCollisionResolution(Obstacle * otherObject) {
		return (void*)(otherObject != parent);
	}
	void resolveCollision(Obstacle * o, void * collisionData) {
		this->alive = false;
		Agent * a = dynamic_cast<Agent *>(o);
		if(a != NULL) {
//			a->alive = false;
			a->color = Random::PRNG() & 0xffffff;
		}
	}
};