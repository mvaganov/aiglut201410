#pragma once

#include "obstacles.h"
#include "codegiraffe/glutrenderingcontext.h"

#include "steering.h"

class Agent : public Obstacle {
public:
	/** where is this? */
	CircleObject body;
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
	void * parent;

	bool alive;

	int agentLook;
	static const int LOOK_PLAIN = 0, LOOK_TRIANGLE = 1, 
		LOOK_BEAK = 2, LOOK_GOOGLY_EYES = 3, LOOK_ROCKET = 4;

	int behavior;
	static const int BEHAVIOR_NONE = 0, BEHAVIOR_SEEK = 1;

	Agent(CircF circle)
		:body(circle),direction(V2f::ZERO_DEGREES()),agentLook(LOOK_BEAK),parent(NULL),
		behavior(BEHAVIOR_SEEK), alive(true), mass(1)
	{}

	virtual void update(int a_ms) {
		switch(behavior) {
		case BEHAVIOR_NONE:
			body.center += velocity * (float)a_ms / 1000.0f;
			break;
		case BEHAVIOR_SEEK:
			acceleration = seek(targetPosition, this);
			velocity += acceleration * (float)a_ms / 1000.0f;
			body.center += velocity * (float)a_ms / 1000.0f;
			break;
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
				float squareDiagonal = sqrt( rr + rr );
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
		g_screen.setColor(0x00ff00);
		g_screen.drawLine(body.center, body.center + velocity);
		g_screen.setColor(0x0000ff);
		g_screen.drawLine(body.center + velocity, 
			body.center + velocity + acceleration);
	}

	bool intersects(const Obstacle * o) const { return body.intersects(o); }
	bool contains(V2f const & p) const { return body.contains(p); }
	bool raycast(V2f const & rayStart, V2f const & rayDirection,
		float & out_dist, V2f & out_point, V2f & out_normal) const {
			return body.raycast(rayStart, rayDirection, out_dist, out_point, out_normal);
	}
	V2f getClosestPointOnEdge(const V2f point, V2f & out_normal) const {
		return body.getClosestPointOnEdge(point, out_normal);
	}
	void glDraw(bool filled) const { body.glDraw(filled); }
	// TODO also calculate transfer of momentum, and reflection
	void * calculateCollisionResolution(Obstacle * otherObject){ 
		float myResponsibilityToMove = 1;
		Agent * a = dynamic_cast<Agent*>(otherObject);
		if (a != NULL) {
			myResponsibilityToMove = mass / (a->mass + mass);
		}
		V2f normal;
		V2f closestPoint = otherObject->getClosestPointOnEdge(body.center, normal);
		V2f positionToClipTo = closestPoint + normal * body.radius;
		V2f clipDelta = positionToClipTo - body.center;
		return new V2f(clipDelta * myResponsibilityToMove);
	}
	void resolveCollision(Obstacle * o, void * collisionData) {
		V2f * push = (V2f*)collisionData;
		body.center += *push;
		delete push;
	}
};