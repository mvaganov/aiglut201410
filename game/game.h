#pragma once

#include "codegiraffe/v2.h"
#include "codegiraffe/templatevector.h"
#include "obstacles.h"
#include "codegiraffe/glutrenderingcontext.h"
#include "agent.h"
#include "codegiraffe/cone.h"

#include <stdlib.h>
#include <stdio.h>

#include "graph.h"
#include "mazegen.h"
#include "astar.h"

#include "codegiraffe/polygon.h"

class Game {
public:
	V2f mousePosition, mouseClick, mouseDragged;
	TemplateVector<Obstacle*> obstacles;
	TemplateVector<Agent*> agents;

	Graph * mapGraph;
	TemplateVector<GraphNode*> * mapPath;

	/** who does the user have selected */
	Agent * selected;
	GraphNode * selectedNode;

	/** whether or not to draw debug lines for FSM steering behaviors */
	bool drawDebug;

	Agent * getAgentAt(V2f const & click) {
		for(int i = 0; i < agents.size(); ++i) {
			if(agents[i]->body.contains(click))
				return agents[i];
		}
		return 0;
	}

	/** generate a list of agents in the given area */
	void gatherListOfAgentsAt(CircF const & area, TemplateVector<Agent*> & out_agents) {
		for(int i = 0; i < agents.size(); ++i) {
			if(V2f::distance(agents[i]->body.center, area.center) 
				< area.radius + agents[i]->body.radius) {
					out_agents.add(agents[i]);
			}
		}
	}

// helpful for finding memory by it's ID number (thanks _CrtSetDbgFlag!)
#define TRACE_MEMORY(mem, debugmessage) \
	printf("memID %d is %s\n", ((int*)mem)[-2], debugmessage)

	void generateWallBoxesForGraph(Graph * g, float wallEdgeValue) {
		GraphNode * n;
		GraphEdge * e;
		const float wallWidth = 1.0f / 20.0f;
		for(int node = 0; node < g->nodes.size(); ++node) {
			n = g->nodes[node];
			for(int i = 0; i < n->edges.size(); ++i) {
				e = &n->edges[i];
				if(e->cost == wallEdgeValue) {
					V2f delta = e->to->location - e->from->location;
					V2f inBetween = e->from->location + (delta * (.5f + wallWidth/2));
					float dist = delta.magnitude();
					V2f normal = delta / dist;
					float width = dist * wallWidth;
					BoxF wall(inBetween, V2f(width, dist), normal);
					obstacles.add(new BoxObject(wall));
				}
			}
		}
	}

	// used for testing object types
	//ConeObject * testcone;
	//BoxObject * testBox;
	//PolygonObject * testPoly;

	Game() :mapGraph(0), mapPath(0), selected(0), selectedNode(0), drawDebug(true) {
		selectedNode = NULL;
		selected = NULL;
		int agentCount = 10;
		// code for testing object types
		//CircF testCircle(V2f(5,5), .75f);
		//RectF aabb(V2f(0, 5), V2f(3, 2));
		//BoxF box(V2f(5,1), V2f(1,3), (float)V_PI / 8);
		//CircleObject * cobj = new CircleObject(testCircle);
		//TRACE_MEMORY(cobj, "circle object");
		//obstacles.add(cobj);
		//obstacles.add(new BoxObject(aabb));
		//testBox = (BoxObject*)obstacles.getLast();
		//TRACE_MEMORY(obstacles.getLast(), "aabb object");
		//obstacles.add(new BoxObject(box));
		//TRACE_MEMORY(obstacles.getLast(), "box object");
		//obstacles.add(new ConeObject(ConeF(V2f(-1, 6), 1, 1.0f, 4.0f)));
		//testcone = (ConeObject*)obstacles.getLast();
		//V2f points[] = { V2f(-1.5f, .9f), V2f(0.0f, 1.0f), V2f(1.2f, .5f), V2f(.2f, -1.0f), V2f(-.9f, -.9f) };
		//const int numPoints = sizeof(points) / sizeof(points[0]);
		//testPoly = new PolygonObject(Polygon2f(V2f(4.0f, 7.0f), 1, points, numPoints));
		//obstacles.add(testPoly);

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
		mapGraph = Graph::createDenseGrid(5, 5, V2f(-10, -10), V2f(10, 10));
		randomMazeGen_prim(mapGraph, -1);
		generateWallBoxesForGraph(mapGraph, -1);

		mapPath = Astar(mapGraph->nodes[0], mapGraph->nodes.getLast());
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
		if (mapGraph) {
			delete mapGraph;
		}
		if (mapPath) {
			delete mapPath;
		}
	}

	/** assumes that the graph will always have at least one node */
	GraphNode * getGraphNodeClosestTo(V2f position) {
		GraphNode * best = mapGraph->nodes[0];
		float shortestDistance = V2f::distance(best->location, position);
		for(int i = 1; i < mapGraph->nodes.size(); ++i) {
			GraphNode * n = mapGraph->nodes[i];
			float dist = V2f::distance(n->location, position);
			if(dist < shortestDistance) {
				shortestDistance = dist;
				best = n;
			}
		}
		return best;
	}

	/** 
	 * @param start where to look from
	 * @param direction what direction to look
	 * @param maxDistance ignore hits that are further than this. Use a negative value to ignore maxDistance
	 * @param dontCareAboutObstacle if true, will return when even one obstacle is found in the way (out_obstacle is not expected to be the closest Obstacle)
	 * @param out_obstacle the obstacle hit while raycasting
	 * @param out_dist how far away the hit happened
	 * @param out_point where the hit happend
	 * @param out_normal the direction of the surface at that hit
	 * @param ignoreList a list of objects to ignore
	 * @param ignoreCount how many objects to ignore
	 * @return true if nothing was hit
	 */
	bool raycast(V2f start, V2f direction, float maxDistance, bool dontCareAboutObstacle, Obstacle * & out_obstacle, float & out_dist, V2f & out_point, V2f & out_normal, Obstacle ** ignoreList, int ignoreCount) {
		float closest = -1;
		out_obstacle = NULL;
		for(int i = 0; i < obstacles.size(); ++i) {
			Obstacle * o = obstacles[i];
			bool checkThisOne = true;
			if(ignoreList != NULL) {
				for(int ignore = 0; ignore < ignoreCount; ++ignore) {
					if(ignoreList[ignore] == o) {
						checkThisOne = false;
						break;
					}
				}
			}
			if(checkThisOne) {
				float dist;
				V2f point, normal;
				if(o->raycast(start, direction, dist, point, normal)
				&& dist > 0 && dist < maxDistance && (closest < 0 || dist < closest)) {
					out_dist = dist;
					out_obstacle = o;
					out_point = point;
					out_normal = normal;
					if (dontCareAboutObstacle) return true; // return ASAP if block/no-block is all that is required
				}
			}
		}
		return (out_obstacle != NULL);
	}

	void display(GLUTRenderingContext & g_screen) {
		g_screen.setColor(0x00aaff);
		mapGraph->glDraw(&g_screen);
		g_screen.printf(mousePosition, "%.2f, %.2f", mousePosition.x, mousePosition.y);
		g_screen.drawCircle(mousePosition, .1f, false);

		// testing cone stuff
		g_screen.drawCircle(mouseClick, .05f, true);
		g_screen.drawLine(mouseClick, mousePosition);
		V2f hit, norm;

		// used for testing object types
		//float dist;
		//Obstacle * obs = testPoly;//testcone;//testBox;
		//if (obs->raycast(mouseClick, (mousePosition - mouseClick).normal(), dist, hit, norm)) {
		//	g_screen.drawCircle(CircF(hit, .1f), false);
		//	g_screen.drawLine(hit, hit + norm);
		//	g_screen.setColor(0);
		//	g_screen.printf(hit+V2f(0,.2f), "%f", dist);
		//}
		//hit = obs->getClosestPointOnEdge(mousePosition, norm);
		//g_screen.drawCircle(CircF(hit, .05f), false);
		//g_screen.drawLine(hit, hit + norm * 0.5);

		/*
		g_screen.setColor(0x00ff00);
		V2f delta = mousePosition - mouseClick; // whereYouAre - whereYouWere
		if (delta.isZero()) delta = V2f::ZERO_DEGREES();
		float len = delta.magnitude();
		float dragLen = 1;// (mouseClick - mouseDragged).magnitude();
		norm = delta.normal();
		float piRad = norm.piRadians();
		ConeF cursorCone = ConeF(mouseClick, len, piRad, piRad + dragLen);
		bool intersect = testcone->intersectsCone(cursorCone);
		cursorCone.glDraw(intersect);
		*/

		if (mapPath != NULL) {
			g_screen.setColor(0x0000ff);
			for (int i = 0; i < mapPath->size(); ++i) {
				g_screen.printf(mapPath->get(i)->location, "%d", i);
			}
		}

		g_screen.setColor(0x008800);
		V2f point, normal;
		for(int i = 0; i < obstacles.size(); ++i) {
			Obstacle * THINGY = obstacles[i];
			THINGY->glDraw(false);
		}
		for(int i = 0; i < agents.size(); ++i) {
			agents[i]->draw(g_screen);
		}
		if (drawDebug) {
			for (int i = 0; i < agents.size(); ++i) {
				if (agents[i]->showDebugLines) {
					agents[i]->drawDebug(g_screen);
				}
			}
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
