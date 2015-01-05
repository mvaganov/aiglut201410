#pragma once

#include "game.h"
#include "mazegen.h"
#include "bullet.h"
#include "obstacles.h"

Obstacle* Game::getSelectedAtPosition(V2f position, long mask) {
	for (int i = 0; i < selectedObstacles.size(); ++i) {
		Obstacle * obs = selectedObstacles[i];
		if (obs->contains(position, mask))
			return obs;
	}
	return 0;
}

void translate(TemplateVector<Obstacle*> & obstacles, V2f moved) {
	for (int i = 0; i < obstacles.size(); ++i) {
		obstacles[i]->getShape()->translate(moved);
	}
}

void rotate(TemplateVector<Obstacle*> & obstacles, float rotationPiRadians) {
	for (int i = 0; i < obstacles.size(); ++i) {
		obstacles[i]->getShape()->rotate(rotationPiRadians);
	}
}

void addToMap(CellSpacePartition * map, TemplateVector<Obstacle*> & obstacles) {
	for (int i = 0; i < obstacles.size(); ++i) {
		map->add(obstacles[i]);
	}
}

void Game::mouse(V2f position, int button, int state) {
	if (button == -1 && state == -1) {
		if (doingOperation) {
			rotationLast = rotationCurrent;
			if (mousePosition != position) {
				V2f opDelta = position - operationPoint;
				printf("%.2f %.2f : ", opDelta.x, opDelta.y);
				if (!opDelta.isZero()) {
					opDelta.normalize();
					//rotationCurrent = opDelta.piRadians(V2f::ZERO());
					rotationCurrent = opDelta.piRadians();
				}
				printf("%f\n", rotationCurrent);
				rotationDelta = rotationCurrent - rotationLast;
			}
			else {
				rotationDelta = 0;
			}
			// TODO make this an event processed by Update...
			if (doingRotation && rotationDelta != 0) {
				objectMap.remove(selectedObstacles.getRawList(), selectedObstacles.size());
				rotate(*selectedObstacles.asVector(), rotationDelta);
				objectMap.add(selectedObstacles.getRawList(), selectedObstacles.size());
			}
		}
		mousePosition = position;
	}

	switch (button) {
	case 3: g_screen->zoom(1 / 1.03125f, position); break;
	case 4: g_screen->zoom(1.03125f, position); break;
	case GLUT_MIDDLE_BUTTON:
		switch (state) {
		case GLUT_DOWN: g_screen->scrollStart(position); break;
		case GLUT_UP:   g_screen->scrollStop(position);  break;
		case 2: g_screen->scrollDrag(position); break;
		}
		break;
	case GLUT_LEFT_BUTTON:
		switch (state) {
		case GLUT_DOWN: {
			if (doingOperation) {
				doingOperation = false;
				break;
			}
			selectedArea.set(position, position);
			Obstacle * o = NULL;
			o = getSelectedAtPosition(position, Obstacle::EVERYTHING);
			if (o == NULL) {
				TemplateSet<Obstacle*> atPosition;
				objectMap.gatherAtPosition(position, atPosition, Obstacle::EVERYTHING);
				if (atPosition.size() != 0) {
					o = atPosition[0];
					// TODO if ctrl is not down
					selectedObstacles.clear();
					selectedObstacles.add(o);
				}
			}
			draggingSelected = (o != NULL);
			if (draggingSelected) {
				mouseDragged = position;
			}
			if (!draggingSelected) {
				selectedObstacles.clear();
				drawingSelectionBox = true;
			}
			refreshSelection();
		} break;
		case GLUT_UP: {
			if (!selectedArea.isValid()) {
				selectedArea.validate();
			}
			if (selectedArea.getWidth() > 0 && selectedArea.getHeight() > 0) {
				objectMap.gatherAtAABB(selectedArea, selectedObstacles, Obstacle::EVERYTHING);
			}
			else {
				objectMap.gatherAtPosition(position, selectedObstacles, Obstacle::EVERYTHING);
			}
			refreshSelection();
			selectedArea.invalidate();
			drawingSelectionBox = false;
		}break;
		case 2: {
			if (drawingSelectionBox)
				selectedArea.max = position;
			else if (draggingSelected) {
				V2f moved = position - mouseDragged;
				if (!moved.isZero() && selectedObstacles.size() > 0) {
					objectMap.remove(selectedObstacles.getRawList(), selectedObstacles.size());
					translate(*selectedObstacles.asVector(), moved);
					objectMap.add(selectedObstacles.getRawList(), selectedObstacles.size());
					refreshSelection();
				}
				mouseDragged = position;
			}
		} break;
		} break;
	case GLUT_RIGHT_BUTTON:
		switch (state) {
		case GLUT_DOWN:
			userRay.start = position;
			break;
		case GLUT_UP: break;
		case 2:
			userRay.direction = position - userRay.start;
			if (userRay.direction.isZero()) {
				userRayDistance = 0;
				userRay.direction = V2f::ZERO_DEGREES();
				//printf("ZERO\n");
			}
			else {
				userRayDistance = userRay.direction.magnitude();
				userRay.direction /= userRayDistance;
				//printf("%.2f\n", userRayDistance);
			}
			break;
		} break;
	}
}

/**
* @param button keycode
* @param state 0: pressed, 1: released, 2: dragged, -1:not actually active
*/
void Game::key(V2f position, int key, int state) {
	V2f delta;
	switch (key) {
	case 27:	gameLoopRunning = false;  break;
	case ' ':
	{
	}
		break;
	case '\t':
		drawDebug = !drawDebug;
		break;
	case 'p':
		paused = !paused;
		break;
	case 'r':	// TODO rotate shapes like in Blender
		operationPoint = selectedOrigin;
		doingOperation = true;
		doingRotation = true;
		delta = position - operationPoint;
		if (!delta.isZero()) {
			delta.normalize();
			rotationCurrent = delta.piRadians();
		} else
			rotationCurrent = 0;
		break;
	case 's':	// TODO scale shapes like in Blender
		operationPoint = position;
		doingOperation = true;
		doingScale = true;
		break;
	case 'g':	// TODO grab/translate shapes like in Blender
		operationPoint = position;
		doingOperation = true;
		doingTranslation = true;
		break;
	case 'd':	// TODO if control is down, duplicate the shapes like in blender
		break;
	}
}

Game::Game() ://mapGraph(NULL), mapPath(NULL), selected(NULL), selectedNode(NULL), delauny(NULL), 
drawDebug(true), objectMap(AABBf(V2f(-50,-50),V2f(50,50)), V2i(7,7)) {
	gameLoopRunning = true;
	paused = false;

	V2f min(-10, -10), max(10, 10);
}

void Game::addObstacle(Obstacle * obstacle) {
	obstacles.add(obstacle);
	objectMap.add(obstacle);
}


void Game::init() {
	objectMap.addLayer(Obstacle::STATIC, AABBf(V2f(-60, -60), V2f(60, 60)), V2i(16, 16));
	addObstacle(new CircleObject(Circf(V2f(-4, -4), 2), Obstacle::STATIC));
	addObstacle(new AABBObject(AABBf(V2f(-4, 1), V2f(-2, 5)), Obstacle::STATIC));
	addObstacle(new BoxObject(Boxf(V2f(4, 0), V2f(2, 5), V2f(1)), Obstacle::STATIC));
	addObstacle(new ConeObject(Conef(V2f(0, 1), 1.5f, (float)V_2PI * 1 / 8, (float)V_2PI * 4 / 8), Obstacle::STATIC)); // make the radius of the cone correct based on the calculated center
	V2f trapezoid[] = { V2f(-.5f, .5f), V2f(.5f, 1.2f), V2f(.5f, -1.5f), V2f(-.5f, -.8f) };
	int trapezoidCount = sizeof(trapezoid) / sizeof(trapezoid[0]);
	addObstacle(new PolygonObject(Polygon2f(V2f(5, -7), 0, trapezoid, trapezoidCount), Obstacle::STATIC)); // TODO raycast does not always hit the shape correctly from all sides
	V2f a(-2, -2), b(5, -4);
	addObstacle(new LineObject(Linef(&a, &b), Obstacle::STATIC));
	addObstacle(new PointObject(V2f(7, 7), Obstacle::STATIC));
}

Game::~Game() {
	for (int i = obstacles.size() - 1; i >= 0; --i) {
		Obstacle * o = obstacles[i];
		obstacles.removeData(o);
		delete o;
	}
}

/**
* @param start where to look from
* @param direction what direction to look
* @param maxDistance ignore hits that are further than this. Use a negative value to ignore maxDistance
* @param dontCareAboutObstacle if true, will return when even one obstacle is found in the way (out_obstacle is not expected to be the closest Obstacle)
* @param out_obstacle the obstacle hit while raycasting
* @param out_rh where the hit happened, how far away the hit happened, the direction of the surface at that hit
* @param mask which layers to check (ignoring all others) for this operation
* @return true if nothing was hit
*/
bool Game::raycast(Ray ray, RaycastHit & out_rh, float maxDistance, bool dontCareAboutObstacle, Obstacle * & out_obstacle, long mask) {
	return objectMap.raycastContainer(ray, out_rh, maxDistance, dontCareAboutObstacle, out_obstacle, mask);
}

void drawRaycastHit(RaycastHit & rh, float scale) {
	glDrawCircle(rh.point, scale / 4, false);
	rh.point.glDrawTo(rh.point + rh.normal*scale);
}

void Game::refreshSelection() {
	selectedOrigin = V2f::ZERO();
	for (int i = 0; i < selectedObstacles.size(); ++i) {
		selectedOrigin += selectedObstacles[i]->getShape()->getOrigin();
	}
	if (selectedObstacles.size() > 1)
		selectedOrigin /= (float)selectedObstacles.size();
}


void Game::draw(GLUTRenderingContext * g) {

	AABBf selectionUIVisual = selectedArea;
	selectionUIVisual.validate();
	g->setColor(0x88ff88);
	selectionUIVisual.glDraw(false);

	if (doingOperation) {
		g->setColor(0xaaaaff);
		g->drawLine(operationPoint, mousePosition);
	}

	// draw selected obstacle indicator
	if (selectedObstacles.size() > 0) {
		g->setColor(0xaaffaa);
		for (int i = 0; i < selectedObstacles.size(); ++i) {
			Shape * s = selectedObstacles[i]->getShape();
			g->drawCircle(s->getCenter(), s->getRadius(), false);
			s->draw(g, true);
		}
		g->drawCircle(selectedOrigin, .1f, true);
	}

	g->setColor(0x00aaff);
	g->printf(mousePosition, "%.2f, %.2f", mousePosition.x, mousePosition.y);
	g->drawCircle(mousePosition, .1f, false);
	int layerColors[] = { 0x0088aa, 0x66aa88, 0x88aa22, 0x123456 };
	int layerColorsCount = sizeof(layerColors) / sizeof(layerColors[0]);
	objectMap.draw(g, layerColors, layerColorsCount);

	// testing cone stuff
	RaycastHit rh;
	V2f hit, norm;
	float dist = 0, maxDist = 0;
	Obstacle * o;
	g->drawCircle(userRay.start, .05f, true);
	// TEST <-- TODO
	V2f delta = mousePosition - userRay.start;
	if (!delta.isZero()) { dist = delta.magnitude(); delta /= dist; }
	maxDist = dist;

	//TemplateSet<int> cellList;
	//staticObstaclesMap->raycastCellList(mouseClick, delta, dist, cellList);
	//for (int i = 0; i < cellList.size(); ++i) { staticObstaclesMap->cells[cellList[i]]->glDraw(true); }

	if (raycast(userRay, rh, userRayDistance, false, o, -1))
	{
		g->setColor(0);
		drawRaycastHit(rh, 1);
		g->printf(rh.point, "%.2f, ", rh.distance);
	}
	g->drawLine(userRay.start, userRay.start + userRay.direction * userRayDistance);//userRay.glDraw();

	for (int i = 0; i < obstacles.size(); ++i) {
		Obstacle * obs = obstacles[i];
		g->setColor(0x333333);
		obs->getShape()->draw(g, false);
		obs->getShape()->getClosestRaycastHit(mousePosition, rh);
		g->setColor(0x888888);
		drawRaycastHit(rh, .2f);
	}
}

void Game::gatherCollisions(Obstacle * s, TemplateSet<Obstacle*> & out_possibleObstacles) {
	objectMap.gatherCollisions(s, out_possibleObstacles);
}

void Game::update(int a_ms) {
	if (paused) return;
	TemplateVector<int> cellIndexList;
	//movingObstaclesMap->clear();
	//for (int i = 0; i < agents.size(); ++i) {
	//	if (agents[i]->alive) {
	//		movingObstaclesMap->add(agents[i]);
	//	}
	//}

	//// update agents by specified amount of time
	//for (int i = 0; i < agents.size(); ++i) {
	//	agents[i]->update(a_ms);
	//}
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
	//for (int a = 0; a < agents.size(); a++) {
	//	Obstacle * agentObstacle = agents[a];//dynamic_cast<Obstacle*>(agents[a]);
	//	gatherCollisions(agentObstacle, hitObstacles);
	//	for (int i = 0; i < hitObstacles.size(); ++i) {
	//		Obstacle * whatWasHit = dynamic_cast<Obstacle*>(hitObstacles[i]);
	//		if (whatWasHit != NULL && agentObstacle->intersects(hitObstacles[i], hitObstacles[i]->getMask()))
	//		{
	//			obstaclecollision collision = {
	//				agentObstacle, whatWasHit,
	//				agentObstacle->calculateCollisionResolution(whatWasHit),
	//			};
	//			collisions.add(collision);
	//		}
	//		//else {
	//		//	printf("non-collision in list...");
	//		//	uselessIters++;
	//		//}
	//	}
	//}
	//printf("%d ", uselessIters);

	// resolve collision seperately, so that one resolved collision will not impact others mid-loop
	for (int i = 0; i < collisions.size(); ++i) {
		obstaclecollision collision = collisions[i];
		if (collision.howAWillRespond) {
			collision.a->resolveCollision(collision.b, collision.howAWillRespond);
		}
	}
	collisions.clear();
	//// garbage collection
	//for (int i = agents.size() - 1; i >= 0; --i) {
	//	if (!agents[i]->alive) {
	//		if (agents[i]->readyToDelete) {
	//			Agent * a = agents[i];
	//			// remove it from the obstacle list too
	//			obstacles.removeData(a);
	//			agents.removeData(a);
	//			delete a;
	//		}
	//		else {
	//			agents[i]->readyToDelete = true;
	//		}
	//	}
	//}
}

// maze generation

//void Game::generateWallBoxesForGraph(AbstractGraph * g, float wallEdgeValue) {
//	AbstractGraphNode * n;
//	AbstractGraphEdge * e;
//	const float wallWidth = 1.0f / 20.0f;
//	for (int node = 0; node < g->getNodeCount(); ++node) {
//		n = g->getNode(node);
//		for (int i = 0; i < n->getEdgeCount(); ++i) {
//			e = n->getEdge(i);
//			if (e->getCost() == wallEdgeValue) {
//				DelaunySet::VoronoiNode * vna = dynamic_cast<DelaunySet::VoronoiNode*>(e->getNode(0));
//				DelaunySet::VoronoiNode * vnb = dynamic_cast<DelaunySet::VoronoiNode*>(e->getNode(1));
//				if (vna && vnb) {
//					DelaunySet::Edge * ve = dynamic_cast<DelaunySet::Edge*>(e);
//					DelaunySet::VoronoiFace * vf = &ve->getFace();
//					V2f a, b;
//					vf->gatherOppositePoints(a, b);
//					V2f delta = b - a;
//					float dist = delta.magnitude();
//					V2f inBetween = a + (delta / 2.0f);
//					Boxf wall(inBetween, V2f(wallWidth, dist), vf->normal);
//					static int wallID = 0;
//					//						if (wallID++ < 700)
//					obstacles.add(new BoxObject(wall));
//				}
//				else {
//					GraphNode * a = dynamic_cast<GraphNode*>(e->getNode(0));
//					GraphNode * b = dynamic_cast<GraphNode*>(e->getNode(1));
//					if (a && b) {
//						V2f delta = b->getLocation() - a->getLocation();
//						V2f inBetween = a->getLocation() + (delta * (.5f + wallWidth / 2));
//						float dist = delta.magnitude();
//						V2f normal = delta / dist;
//						float width = dist * wallWidth;
//						//						Boxf wall(inBetween, V2f(width, dist), normal);
//						Boxf wall(inBetween, V2f(wallWidth, 1), normal);
//						obstacles.add(new BoxObject(wall));
//					}
//				}
//			}
//		}
//	}
//}
