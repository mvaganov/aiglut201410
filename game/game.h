#pragma once

#include "codegiraffe/v2.h"
#include "codegiraffe/templatevector.h"
#include "obstacles.h"
#include "codegiraffe/glutrenderingcontext.h"
#include "agent.h"

#include <stdlib.h>
#include <stdio.h>

class Game {
public:
	V2f mousePosition, mouseClick, mouseDragged;
	TemplateVector<Obstacle*> obstacles;
	TemplateVector<Agent*> agents;
	/** who does the user have selected */
	Agent * selected;

	Agent * getAgentAt(V2f const & click) {
		for(int i = 0; i < agents.size(); ++i) {
			if(agents[i]->body.contains(click))
				return agents[i];
		}
		return 0;
	}

// helpful for finding memory by it's ID number (thanks _CrtSetDbgFlag!)
#define TRACE_MEMORY(mem, debugmessage) printf("memID %d is %s\n", ((int*)mem)[-2], debugmessage)

	Game() {
		selected = NULL;
		int agentCount = 10;
		CircF testCircle(V2f(5,5), .75f);
		RectF aabb(V2f(0, 5), V2f(3, 2));
		BoxF box(V2f(5,1), V2f(1,3), (float)V_PI / 8);
		obstacles.add(new CircleObject(testCircle));
		TRACE_MEMORY(obstacles.getLast(), "circle object");
		obstacles.add(new BoxObject(aabb));
		TRACE_MEMORY(obstacles.getLast(), "aabb object");
		obstacles.add(new BoxObject(box));
		TRACE_MEMORY(obstacles.getLast(), "box object");

		for(int i = 0; i < agentCount; ++i) {
			float extraRadius = Random::PRNGf()*0.5f;
			CircF c(Random::PRNGf() * 5, Random::PRNGf() * 5, .1f + extraRadius);
			Agent * a = new Agent(c, this, new FSM_Idle());
			TRACE_MEMORY(a, "agent");
			a->direction = V2f::randomUnitVector();
			a->maximumSpeed = Random::PRNGf(.1f, 2);
			a->maximumForce = Random::PRNGf(.1f, 10);
			a->velocity = a->direction * a->maximumSpeed;
			a->mass = Random::PRNGf(.1f, 100);
			a->color = Random::PRNG() & 0xffffff;
			a->behavior = Agent::BEHAVIOR_AGGRO;
			agents.add(a);
			obstacles.add(a);
		}
	}
	~Game() {
		//for(int i = 0; i < agents.size(); ++i) {
		for(int i = agents.size()-1; i >= 0; --i) {
			Agent * a = agents[i];
			agents.removeData(a);
			obstacles.removeData(a);
			delete a;
		}
		for(int i = obstacles.size()-1; i >= 0; --i) {
			Obstacle * o = obstacles[i];
			obstacles.removeData(o);
			delete o;
		}
	}
	void display(GLUTRenderingContext & g_screen) {
		g_screen.printf(mousePosition, "%.2f, %.2f", mousePosition.x, mousePosition.y);
		g_screen.printf(mouseClick, "%.2f, %.2f", mouseClick.x, mouseClick.y);
		g_screen.drawCircle(mousePosition, .1f, false);
		g_screen.drawCircle(mouseClick, .05f, true);
		g_screen.setColor(0x0000ff);
		V2f delta = mouseDragged - mouseClick; // whereYouAre - whereYouWere
		delta.normalize();
		g_screen.printf(mouseDragged, "%.2f, %.2f m:%.2f", delta.x, delta.y, delta.magnitude());
		g_screen.drawLine(mouseClick, mouseClick + delta);
		V2f a(0,2), b(1,0);
		g_screen.setColor(0xffaa00);
		g_screen.drawLine(a, b);
		float dist;
		V2f collisionPoint;
		if(V2f::lineIntersection(mouseClick, mouseDragged, a, b, dist, collisionPoint)) {
			g_screen.drawCircle(collisionPoint, .1f, false);
			g_screen.printf(collisionPoint, "%.2f", dist);
		}
		g_screen.setColor(0x008800);
		V2f point, normal;
		for(int i = 0; i < obstacles.size(); ++i) {
			Obstacle * THINGY = obstacles[i];
			THINGY->glDraw(false);
			if(THINGY->contains(mousePosition)) THINGY->glDraw(true);
			//point = THINGY->getClosestPointOnEdge(mousePosition, normal);
			//g_screen.drawLine(point, point+normal);
			//if(THINGY->raycast(mouseClick, delta, dist, point, normal)) {
			//	g_screen.drawLine(point, point-delta);
			//	V2f reflected;
			//	V2f::calculateReflection(delta, normal, reflected);
			//	g_screen.drawLine(point, point+reflected);
			//}
		}
		for(int i = 0; i < agents.size(); ++i) {
			agents[i]->draw(g_screen);
		}
		if(selected != NULL) {
			g_screen.setColor(0x00ff00);
			g_screen.drawCircle(selected->body.center, selected->body.radius+.1f, false);
		}
	}
	void update(int a_ms) {
		// update agents by specified amount of time
		for(int i = 0; i < agents.size(); ++i) {
			agents[i]->update(a_ms);
		}
		// detect and calculate collisions
		struct obstaclecollision {
			// TODO check-for and handle multiple collisions to the same obstacle
			Obstacle * a;
			Obstacle * b;
			void * howAWillRespond;
		};
		TemplateVector<obstaclecollision> collisions;
		for(int a = 0; a < agents.size(); a++) {
			for(int b = 0; b < obstacles.size(); ++b) {
				if(agents[a] != obstacles[b]
				&& agents[a]->intersects(obstacles[b])) {
					obstaclecollision collision = { 
						agents[a], obstacles[b],
						agents[a]->calculateCollisionResolution(obstacles[b]),
					};
					collisions.add(collision);
				}
			}
		}
		// resolve collision seperately, so that one resolved collision will not impact others mid-loop
		for (int i = 0; i < collisions.size(); ++i) {
			obstaclecollision collision = collisions[i];
			if (collision.howAWillRespond)
				collision.a->resolveCollision(collision.b, collision.howAWillRespond);
		}
		// garbage collection
		for(int i = agents.size()-1; i >= 0; --i) {
			if (!agents[i]->alive) {
				if (agents[i]->readyToDelete) {
					Agent * a = agents[i];
					// remove it from the obstacle list too
					obstacles.removeData(a);
					agents.removeData(a);
					delete a;
				} else {
					agents[i]->readyToDelete = true;
				}
			}
		}
	}
};
