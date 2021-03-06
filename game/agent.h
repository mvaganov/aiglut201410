#pragma once

#include "obstacles.h"
#include "codegiraffe/glutrenderingcontext.h"

#include "steering.h"
#include "fsm.h"
#include "fsm_idle.h"
#include <stdio.h>

class Game; // class prototype, forward declaration "I promise game is a class."
class Bullet;

class Agent : public ObstacleConcrete {
public:
	const char* getShapeName() const { return "Agent"; }
	/** where is this? */
	ShapeCircle body;
	/** where is this going? */
	V2f velocity;
	/** where is this facing? (cos, sin) of the angle this is facing */
	V2f direction;
	/** a 4 byte color value */
	long color;
	/** how fast our agents are allowed to go */
	float maximumSpeed;
	/** how fast our agents are allowed to go */
	float maximumForce;
	/** the current force being applied to the velocity */
	V2f acceleration;
	/** where the agent is thinking about for steering behaviors */
	V2f targetPosition;
	/** used for physics and collision calculations */
	float mass;

	/** where did  this agent come from? */
	Obstacle * parent;
	/** whether or not to keep this in the game (or delete from lists) */
	bool alive;
	/**
	 * the flag that will actually delete an object. this allows a dead object
	 * to make it through a game loop *one last time*, so that things pointing
	 * at it can make preparations
	 */
	bool readyToDelete;
	/** true if the player is controlling this guy */
	bool playerControlled;
	/** whether or not to draw debug information */
	bool showDebugLines;
	/** reference to the game that this agent is in */
	Game * game;
	/** the finite state machine */
	FiniteStateMachine * fsm;

	BoxObject sensorArea;

	Shape * getShape() { return &body; }
	const Shape * getShape() const { return &body; }

	void setFSM(FiniteStateMachine * a_fsm) {
		if(fsm) {
			fsm->exit(this);
			delete fsm;
		}
		fsm = a_fsm;
		if(fsm) {
			fsm->enter(this);
		}
	}

	int agentLook;
	static const int LOOK_PLAIN = 0, LOOK_TRIANGLE = 1, 
		LOOK_BEAK = 2, LOOK_GOOGLY_EYES = 3, LOOK_ROCKET = 4;

	int behavior;
	static const int BEHAVIOR_NONE = 0, BEHAVIOR_SEEK = 1, BEHAVIOR_AGGRO = 2;

	Agent(Circf circle, Game * game, FiniteStateMachine * newFSMBehavior)
		:body(circle),direction(V2f::ZERO_DEGREES()),agentLook(LOOK_BEAK),parent(NULL),
		behavior(BEHAVIOR_SEEK), alive(true), readyToDelete(false), mass(1), 
		playerControlled(false), showDebugLines(true), game(game), fsm(newFSMBehavior),
		sensorArea(Boxf(circle.center, V2f::ZERO(), 0), Obstacle::STATIC)
	{ }

	~Agent() { if (fsm != NULL) delete fsm;  }

	/** agent behavior, see agent.cpp */
	virtual void update(int a_ms);

	V2f getCenter() const {return body.center; }
	float getRadius() const {return body.radius;}

	void updateMovement(int a_ms) {
		if (a_ms != 0) {
			if(!acceleration.isZero()) {
				acceleration.truncate(maximumForce);
				// check the speed, to adjust acceleration so that it doesn't break the speed limit
				float currentSpeed = velocity.magnitude();
				if (currentSpeed >= maximumSpeed) {
					// reduce the acceleration by what it would be contributing to the speed
					V2f vDir = velocity / currentSpeed;
					float accelAlignmentWithVelocity = V2f::dot(vDir, acceleration);
					if (accelAlignmentWithVelocity > 0) {
						acceleration -= vDir * accelAlignmentWithVelocity;
					}
				}
			}
			velocity += acceleration * (float)a_ms / 1000.0f;
			body.center += velocity * (float)a_ms / 1000.0f;
			sensorArea.getShapeBox()->set(velocity / 2 + body.center,
				V2f(velocity.magnitude(), body.radius*2), velocity.normal());
		}
	}

	void draw(GLUTRenderingContext & g_screen) {
		g_screen.setColor(color);
		switch(agentLook) {
		case LOOK_PLAIN:
			{
				body.glDraw(false);
				V2f center = body.center;
				V2f end = body.center + direction * body.radius;
				g_screen.drawLine(center, end);
			}
			break;
		case LOOK_TRIANGLE:
			{
				V2f center = body.center;
				V2f right = direction * body.radius;//(body.radius, 0);
				V2f top = right.perp();//(0, body.radius);
				V2f bottom = -top;//(0, -body.radius);
				//top.rotate(direction);
				//bottom.rotate(direction);
				right.rotate(direction);
				top += body.center;
				bottom += body.center;
				right += body.center;
				//glDrawCircle(right, .05f, true);
				//glDrawCircle(top, .05f, true);
				//glDrawCircle(bottom, .05f, true);
				g_screen.drawLine(top, right);
				g_screen.drawLine(top, bottom);
				g_screen.drawLine(right, bottom);
				body.glDraw(false);
			}
			break;
		case LOOK_BEAK:
			{
				V2f right = direction * body.radius;
				V2f top = right.rotated((float)V_PI / 4);
				V2f bottom = right.rotated((float)-V_PI / 4);
				float rr = body.radius * body.radius;
				float squareDiagonal = (float)sqrt( rr + rr );
				right = direction * squareDiagonal;
				right += body.center;
				top += body.center;
				bottom += body.center;
				//glDrawCircle(right, .05f, true);
				//glDrawCircle(top, .05f, true);
				//glDrawCircle(bottom, .05f, true);
				g_screen.drawLine(top, right);
				g_screen.drawLine(right, bottom);
				body.glDraw(false);
			}
			break;
		case LOOK_GOOGLY_EYES:
			{
				V2f center = body.center;
				V2f center2 = (direction * body.radius) / 2 + center;
				body.glDraw(false);
				g_screen.drawCircle(center2, body.radius/2, true);
			}
			break;
		}
	}

	void drawDebug(GLUTRenderingContext * g) {
		g->setColor(0x00ff00);
		g->drawLine(body.center, body.center + velocity);
		g->setColor(0x0000ff);
		g->drawLine(body.center + velocity,
			body.center + velocity + acceleration);
		if (fsm != NULL) {
			g->setColor(0);
			g->printf(body.center + V2f(-body.radius, 0), fsm->getName());
			g->setColor(color);
			fsm->draw(this, g);
		}
		g->setColor(0xff00ff);
		sensorArea.draw(g, false);
		if (behavior == BEHAVIOR_AGGRO) {
			g->setColor(0x0000ff);
			g->drawCircle(body.center, body.radius * 5, false);
		}
	}

	void glDraw(bool filled) const { body.glDraw(filled); }

	// TODO also calculate transfer of momentum, and reflection
	void * calculateCollisionResolution(Collidable * otherObject){
		float myResponsibilityToMove = 1;
		Agent * a = dynamic_cast<Agent*>(otherObject);
		if (a == this)
			return 0;
		if (a != NULL) {
			float massSum = a->mass + mass;
			if (massSum == 0) myResponsibilityToMove = 0.5f;
			else myResponsibilityToMove = (massSum-mass) / (massSum);
		}
		V2f normal;
		Shape * s = otherObject->getShape();
		RaycastHit closest;
		s->getClosestRaycastHit(body.center, closest);
		//printf("%f %f ", normal.x, normal.y);
		V2f positionToClipTo = closest.point + normal * body.radius;
		V2f clipDelta = positionToClipTo - body.center;
		V2f * output = (new V2f(clipDelta * myResponsibilityToMove));
//#define TRACE_MEMORY1(mem, debugmessage) \
//	printf("memID %d is %s\n", ((int*)mem)[-2], debugmessage)
//		TRACE_MEMORY1(output, "v2f otput");
//#undef TRACE_MEMORY1
		return output;
	}
	void resolveCollision(Collidable * o, void * & collisionData) {
		if (!collisionData) return;
		V2f * push = (V2f*)collisionData;
		body.center += *push;
		delete push;
		collisionData = 0;
	}

	Bullet * findClosestBullet();
	Agent * findClosestPlayerControlledAgent();

	/**
	 * @param subject what this Agent is trying to see
	 * @return true if no Obstacles interrupt this Agent's view of the given subject.
	 */
	bool hasLineOfSight(Obstacle * subject);
};
