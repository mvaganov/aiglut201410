#pragma once

#include "codegiraffe/v2.h"
#include "codegiraffe/templatevector.h"
#include "obstacles.h"
#include "codegiraffe/glutrenderingcontext.h"
#include "agent.h"

#include <stdlib.h>

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

	Game() {
		selected = NULL;
		int agentCount = 10;
		CircF testCircle(V2f(5,5), .75f);
		RectF aabb(V2f(0, 5), V2f(3, 2));
		BoxF box(V2f(5,1), V2f(1,3), (float)V_PI / 8);
		obstacles.add(new CircleObject(testCircle));
		obstacles.add(new BoxObject(aabb));
		obstacles.add(new BoxObject(box));

		for(int i = 0; i < agentCount; ++i) {
			float extraRadius = Random::PRNGf()*0.5f;
			CircF c(Random::PRNGf() * 5, Random::PRNGf() * 5, .1f + extraRadius);
			Agent * a = new Agent(c);
			a->direction = V2f::randomUnitVector();
			a->velocity = a->direction;
			a->color = Random::PRNG() & 0xffffff;
			agents.add(a);
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
			point = THINGY->getClosestPointOnEdge(mousePosition, normal);
			g_screen.drawLine(point, point+normal);
			if(THINGY->raycast(mouseClick, delta, dist, point, normal)) {
				g_screen.drawLine(point, point-delta);
				V2f reflected;
				V2f::calculateReflection(delta, normal, reflected);
				g_screen.drawLine(point, point+reflected);
			}
		}
		for(int i = 0; i < agents.size(); ++i) {
			agents[i]->draw(g_screen);
		}
	}
	void update(int a_ms) {
		for(int i = 0; i < agents.size(); ++i) {
			agents[i]->update(a_ms);
		}
	}
};
