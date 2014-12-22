#pragma once

#include "game.h"
#include "mazegen.h"

Agent * Game::getAgentAt(V2f const & click) {
	//for(int i = 0; i < agents.size(); ++i) {
	//	if(agents[i]->body.contains(click))
	//		return agents[i];
	//}
	TemplateSet<Obstacle*> result;
	printf("%.1f, %.1f\n", click.x, click.y);
	movingObstacles->gatherAt(click, result);
	if (result.size() > 0)
		return dynamic_cast<Agent*>(result[0]);
	return 0;
}

void Game::gatherStaticObstaclesAt(CircF const & area, TemplateSet<Obstacle*> & out_obstacles) {
	staticObstacles->gatherAt(area, out_obstacles);
}

/** generate a list of agents in the given area */
void Game::gatherListOfAgentsAt(CircF const & area, TemplateVector<Agent*> & out_agents) {
	//for(int i = 0; i < agents.size(); ++i) {
	//	if(V2f::distance(agents[i]->body.center, area.center) 
	//		< area.radius + agents[i]->body.radius) {
	//			out_agents.add(agents[i]);
	//	}
	//}
	TemplateSet<Obstacle*> result;
	movingObstacles->gatherAt(area, result);
	for (int i = 0; i < result.size(); ++i) {
		Agent * a = dynamic_cast<Agent*>(result[i]);
		if (a) { out_agents.add(a); }
	}
}

void Game::generateWallBoxesForGraph(AbstractGraph * g, float wallEdgeValue) {
	AbstractGraphNode * n;
	AbstractGraphEdge * e;
	const float wallWidth = 1.0f / 20.0f;
	for (int node = 0; node < g->getNodeCount(); ++node) {
		n = g->getNode(node);
		for (int i = 0; i < n->getEdgeCount(); ++i) {
			e = n->getEdge(i);
			if (e->getCost() == wallEdgeValue) {
				DelaunySet::VoronoiNode * vna = dynamic_cast<DelaunySet::VoronoiNode*>(e->getNode(0));
				DelaunySet::VoronoiNode * vnb = dynamic_cast<DelaunySet::VoronoiNode*>(e->getNode(1));
				if (vna && vnb) {
					DelaunySet::Edge * ve = dynamic_cast<DelaunySet::Edge*>(e);
					DelaunySet::VoronoiFace * vf = &ve->getFace();
					V2f a, b;
					vf->gatherOppositePoints(a, b);
					V2f delta = b - a;
					float dist = delta.magnitude();
					V2f inBetween = a + (delta / 2.0f);
					BoxF wall(inBetween, V2f(wallWidth, dist), vf->normal);
					printf("%.1f %.1f\n", wall.size.x, wall.size.y);
					static int wallID = 0;
					//						if (wallID++ < 700)
					obstacles.add(new BoxObject(wall));
				}
				else {
					GraphNode * a = dynamic_cast<GraphNode*>(e->getNode(0));
					GraphNode * b = dynamic_cast<GraphNode*>(e->getNode(1));
					if (a && b) {
						V2f delta = b->getLocation() - a->getLocation();
						V2f inBetween = a->getLocation() + (delta * (.5f + wallWidth / 2));
						float dist = delta.magnitude();
						V2f normal = delta / dist;
						float width = dist * wallWidth;
						//						BoxF wall(inBetween, V2f(width, dist), normal);
						BoxF wall(inBetween, V2f(wallWidth, 1), normal);
						obstacles.add(new BoxObject(wall));
					}
				}
			}
		}
	}
}

Game::Game() :mapGraph(NULL), mapPath(NULL), selected(NULL), selectedNode(NULL), delauny(NULL), drawDebug(true), staticObstacles(NULL), movingObstacles(NULL) {
	delaunyEdit = NULL;
	selectedNode = NULL;
	selected = NULL;
#ifdef TESTING_SHAPES
	// code for testing object types
	CircF testCircle(V2f(5, 5), .75f);
	BoxF box(V2f(5, 1), V2f(1, 3), (float)V_PI / 8);
	this->testCircle = new CircleObject(testCircle);
	//TRACE_MEMORY(cobj, "circle object");
	obstacles.add(this->testCircle);
	obstacles.add(new BoxObject(box));
	testBox = (BoxObject*)obstacles.getLast();
	//TRACE_MEMORY(obstacles.getLast(), "aabb object");
	//obstacles.add(new BoxObject(box));
	//TRACE_MEMORY(obstacles.getLast(), "box object");
	obstacles.add(new ConeObject(ConeF(V2f(-1, 6), 1, 1.0f, 4.0f)));
	testcone = (ConeObject*)obstacles.getLast();

	V2f points[] = { V2f(-1.5f, .9f), V2f(0.0f, 1.0f), V2f(1.2f, .5f), V2f(.2f, -1.0f), V2f(-.9f, -.9f) };
	const int numPoints = sizeof(points) / sizeof(points[0]);
	testPoly = new PolygonObject(Polygon2f(V2f(4.0f, 7.0f), 1, points, numPoints));
	obstacles.add(testPoly);

	//		V2f points2[] = { V2f(-2.5f, 2.9f), V2f(1.2f, 2.5f), V2f(.2f, -1.0f), V2f(-1.9f, -.9f) };
	//		const int numPoints2 = sizeof(points2) / sizeof(points2[0]);
	//		testPoly2 = new PolygonObject(Polygon2f(V2f(4.0f, 7.0f), 1, points2, numPoints2));
	BoxF box2 = BoxF(V2f(2, -3), V2f(.3f, 2), .1f);
	testPoly2 = new PolygonObject(Polygon2f(box2));
	obstacles.add(testPoly2);
#endif
	int agentCount = 10;
	for (int i = 0; i < agentCount; ++i) {
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
	V2f min(-10, -10), max(10, 10);

	//		mapGraph = Graph::createDenseGrid(5, 5, min, max);
	//
	//		randomMazeGen_prim(mapGraph, -1);
	////		generateWallBoxesForGraph(mapGraph, -1);
	//
	//		mapPath = Astar(&mapGraph->nodes[0], &mapGraph->nodes.getLast());

#ifdef DELAUNY_TESTING
	V2f boxPoints[] = { V2f(-10, 10), V2f(10, 10), V2f(10, -10), V2f(-10, -10) };
	const int boxPointCount = sizeof(boxPoints) / sizeof(boxPoints[0]);
	Polygon2f::pair boxSides[] = {
		Polygon2f::pair(0, 1),
		Polygon2f::pair(1, 2),
		Polygon2f::pair(2, 3),
		Polygon2f::pair(3, 0),
	};
	const int boxSidesCount = sizeof(boxSides) / sizeof(boxSides[0]);
	Obstacle * po =
		//new PolygonObject(Polygon2f(V2f::ZERO(), 0, boxPoints, boxPointCount, boxSides, boxSidesCount));
		new CircleObject(CircF(V2f::ZERO(), 50));
	delauny = new DelaunySet(po);
	delauny->makeRandom(10);
	//delauny->currentNodes.add(DelaunySet::VoronoiNode(V2f(4, 4)));
	//delauny->currentNodes.add(DelaunySet::VoronoiNode(V2f(5, 4)));
	//delauny->currentNodes.add(DelaunySet::VoronoiNode(V2f(4, 5)));
	//delauny->currentNodes.add(DelaunySet::VoronoiNode(V2f(5, 5)));
	//V2f p(5, -5);
	//delauny->currentNodes.add(DelaunySet::VoronoiNode(p + V2f(0.0f) * 2));
	//delauny->currentNodes.add(DelaunySet::VoronoiNode(p + V2f(1.0f) * 2));
	//delauny->currentNodes.add(DelaunySet::VoronoiNode(p + V2f(2.0f) * 2));
	//delauny->currentNodes.add(DelaunySet::VoronoiNode(p + V2f(3.0f) * 2));
	//delauny->currentNodes.add(DelaunySet::VoronoiNode(p + V2f(4.0f) * 2));
//			delaunyEdit = delauny;

	delauny->calculateAllTriangles();
	delauny->gatherVoronoi(voronoiNodes, true);
	randomMazeGen_prim(delauny, -1, 2);

	generateWallBoxesForGraph(delauny, -1);
#endif
	staticObstacles = new CellSpacePartition(RectF(V2f(0, 0), 100), V2i(10, 10));
	movingObstacles = new CellSpacePartition(RectF(V2f(0, 0), 200), V2i(5, 5));
	staticObstacles->clear();
	for (int i = 0; i < obstacles.size(); ++i) {
		if (obstacles[i] != NULL && dynamic_cast<Agent*>(obstacles[i]) == NULL)
			staticObstacles->add(obstacles[i]);
	}
}
Game::~Game() {
	//for(int i = 0; i < agents.size(); ++i) {
	for (int i = agents.size() - 1; i >= 0; --i) {
		Agent * a = agents[i];
		agents.removeData(a);
		obstacles.removeData(a);
		delete a;
	}
	for (int i = obstacles.size() - 1; i >= 0; --i) {
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
	if (delauny) {
		if (delauny->boundary){
			delete delauny->boundary;
		}
		delete delauny;
	}
	if (staticObstacles) delete staticObstacles;
	if (movingObstacles) delete movingObstacles;
}

/** assumes that the graph will always have at least one node */
GraphNode * Game::getGraphNodeClosestTo(V2f position) {
	GraphNode * best = &mapGraph->nodes[0];
	float shortestDistance = V2f::distance(best->getLocation(), position);
	for (int i = 1; i < mapGraph->nodes.size(); ++i) {
		GraphNode * n = &mapGraph->nodes[i];
		float dist = V2f::distance(n->getLocation(), position);
		if (dist < shortestDistance) {
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
bool Game::raycast(Ray ray, RaycastHit & out_rh, float maxDistance, bool dontCareAboutObstacle, Obstacle * & out_obstacle, Obstacle ** ignoreList, int ignoreCount) {
	bool hit = false;
	hit = staticObstacles->raycastContainer(ray, out_rh, maxDistance, dontCareAboutObstacle, out_obstacle, ignoreList, ignoreCount);
	if (hit) return true;
	hit = movingObstacles->raycastContainer(ray, out_rh, maxDistance, dontCareAboutObstacle, out_obstacle, ignoreList, ignoreCount);
	if (hit) return true;
	return false;
}

void Game::display(GLUTRenderingContext & g_screen) {
	if (delauny != NULL) {
		g_screen.setColor(0xdddddd);
		if (delaunyEdit != NULL) {
			delaunyEdit->glDraw(g_screen);
			return;
		}
		// draw the voronoi cells
//		delauny->glDraw(g_screen);
		Polygon2f * p2f;
		for (int i = 0; i < voronoiNodes.size(); ++i) {
			p2f = &voronoiNodes[i]->getPolygon2f(delauny->boundary);
			bool filled = p2f->contains(mousePosition);// false;
			p2f->glDraw(filled);
		}
	}

	if (staticObstacles) {
		g_screen.setColor(0x88ff88);
		staticObstacles->glDraw();
		for (int i = 0; i < staticObstacles->cells.size(); ++i) {
			g_screen.printf(staticObstacles->cells[i]->getShape()->getCenter(), "%d", i);
		}
	}
	if (movingObstacles) {
		g_screen.setColor(0xff88ff);
		movingObstacles->glDraw();
		for (int i = 0; i < movingObstacles->cells.size(); ++i) {
			g_screen.printf(movingObstacles->cells[i]->getShape()->getCenter(), "%d", i);
		}
	}

	g_screen.setColor(0x00aaff);
	if (mapGraph) mapGraph->glDraw(&g_screen);
	g_screen.printf(mousePosition, "%.2f, %.2f", mousePosition.x, mousePosition.y);
	g_screen.drawCircle(mousePosition, .1f, false);


	// testing cone stuff
	RaycastHit rh;
	V2f hit, norm;
	float dist = 0, maxDist = 0;
	Obstacle * o;
	g_screen.drawCircle(mouseClick, .05f, true);
	// TEST <-- TODO
	V2f delta = mousePosition - mouseClick;
	if (!delta.isZero()) { dist = delta.magnitude(); delta /= dist; }
	maxDist = dist;

	//TemplateSet<int> cellList;
	//staticObstacles->raycastCellList(mouseClick, delta, dist, cellList);
	//for (int i = 0; i < cellList.size(); ++i) { staticObstacles->cells[cellList[i]]->glDraw(true); }

	if (raycast(Ray(mouseClick, delta), rh, maxDist, false, o, 0, 0))
	{
		g_screen.setColor(0);
		g_screen.drawCircle(rh.point, 1, false);
		g_screen.drawLine(rh.point, rh.point + rh.normal * 3);
		g_screen.printf(rh.point, "%.2f, ", rh.distance);
		//			o->glDraw(true);
	}
	g_screen.drawLine(mouseClick, mousePosition);

#ifdef TESTING_SHAPES
	// used for testing object types
	Obstacle * obs = testcone;//testBox;//testPoly;//
	if (obs->getShape()->raycast(Ray(mouseClick, (mousePosition - mouseClick).normal()), rh)) {
		g_screen.drawCircle(CircF(rh.point, .1f), false);
		g_screen.drawLine(rh.point, rh.point + rh.normal);
		g_screen.setColor(0);
		g_screen.printf(rh.point + V2f(0, .2f), "%f", rh.distance);
	}
	hit = obs->getShape()->getClosestPointOnEdge(mousePosition, norm);
	g_screen.drawCircle(CircF(hit, .05f), false);
	g_screen.drawLine(hit, hit + norm * 0.5);


	g_screen.setColor(0x00ff00);
	delta = mousePosition - mouseClick; // whereYouAre - whereYouWere
	if (delta.isZero()) delta = V2f::ZERO_DEGREES();
	float len = delta.magnitude();
	float dragLen = 1;// (mouseClick - mouseDragged).magnitude();
	norm = delta.normal();
	testPoly->getPolygon()->setRotation(norm.piRadians());
	float piRad = norm.piRadians();
	ConeF cursorCone = ConeF(mouseClick, len, piRad, piRad + dragLen);
	bool intersect;
	intersect = cursorCone.intersectsCircle(testCircle->getCircle()->center, testCircle->getCircle()->radius);
	cursorCone.glDraw(intersect);

	//		testPoly->origin = mousePosition;

	intersect = testPoly->intersects(testBox);
	g_screen.setColor(intersect ? 0x0000ff : 0x00ff00);
	g_screen.drawCircle(mousePosition, 1, intersect);
	testBox->glDraw(intersect);

	testPoly2->glDraw(false);
#endif

	if (mapPath != NULL) {
		g_screen.setColor(0x0000ff);
		for (int i = 0; i < mapPath->size(); ++i) {
			GraphNode * gn = dynamic_cast<GraphNode*>(mapPath->get(i));
			g_screen.printf(gn->getLocation(), "%d", i);
		}
	}

	g_screen.setColor(0x008800);
	V2f point, normal;
	for (int i = 0; i < obstacles.size(); ++i) {
		Obstacle * THINGY = obstacles[i];
		THINGY->getShape()->glDraw(false);
	}
	for (int i = 0; i < agents.size(); ++i) {
		agents[i]->draw(g_screen);
	}
	if (drawDebug) {
		for (int i = 0; i < agents.size(); ++i) {
			if (agents[i]->showDebugLines) {
				agents[i]->drawDebug(g_screen);
			}
		}
	}
	if (selected != NULL) {
		g_screen.setColor(0x00ff00);
		g_screen.drawCircle(selected->body.center, selected->body.radius + .1f, false);
	}
}

void Game::gatherCollisions(Obstacle * s, TemplateSet<Obstacle*> & out_possibleObstacles) {
	staticObstacles->gatherCollisions(s, out_possibleObstacles);
	movingObstacles->gatherCollisions(s, out_possibleObstacles);
}

void Game::update(int a_ms) {
	TemplateVector<int> cellIndexList;
	movingObstacles->clear();
	for (int i = 0; i < agents.size(); ++i) {
		if (agents[i]->alive) {
			movingObstacles->add(agents[i]);
		}
	}

	// update agents by specified amount of time
	for (int i = 0; i < agents.size(); ++i) {
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
	//int uselessIters = 0;
	TemplateSet<Obstacle*> hitObstacles;
	for (int a = 0; a < agents.size(); a++) {
		Obstacle * agentObstacle = dynamic_cast<Obstacle*>(agents[a]);
		gatherCollisions(agentObstacle, hitObstacles);
		for (int i = 0; i < hitObstacles.size(); ++i) {
			Obstacle * whatWasHit = dynamic_cast<Obstacle*>(hitObstacles[i]);
			if (whatWasHit != NULL && agentObstacle->intersects(hitObstacles[i]))
			{
				obstaclecollision collision = {
					agentObstacle, whatWasHit,
					agentObstacle->calculateCollisionResolution(whatWasHit),
				};
				collisions.add(collision);
			}
			//else {
			//	printf("non-collision in list...");
			//	uselessIters++;
			//}
		}
	}
	//printf("%d ", uselessIters);

	// resolve collision seperately, so that one resolved collision will not impact others mid-loop
	for (int i = 0; i < collisions.size(); ++i) {
		obstaclecollision collision = collisions[i];
		if (collision.howAWillRespond) {
			collision.a->resolveCollision(collision.b, collision.howAWillRespond);
		}
	}
	collisions.clear();
	// garbage collection
	for (int i = agents.size() - 1; i >= 0; --i) {
		if (!agents[i]->alive) {
			if (agents[i]->readyToDelete) {
				Agent * a = agents[i];
				// remove it from the obstacle list too
				obstacles.removeData(a);
				agents.removeData(a);
				delete a;
			}
			else {
				agents[i]->readyToDelete = true;
			}
		}
	}
}
