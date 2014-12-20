#pragma once

#include "codegiraffe/v2.h"
#include "codegiraffe/templatevector.h"
#include "obstacles.h"
#include "codegiraffe/glutrenderingcontext.h"
#include "agent.h"
#include "codegiraffe/cone.h"
#include "cellspacepartition.h"

#include <stdlib.h>
#include <stdio.h>

#include "graph.h"
#include "astar.h"
#include "codegiraffe/polygon.h"
#include "delauny.h"
#define DELAUNY_TESTING
//#define TESTING_SHAPES

class Game {
public:
	V2f mousePosition, mouseClick, mouseDragged;
	TemplateVector<Obstacle*> obstacles;
	TemplateVector<Agent*> agents;

	CellSpacePartition<Obstacle> * staticObstacles;
	CellSpacePartition<Obstacle> * movingObstacles;

	Graph * mapGraph;
	TemplateVector<AbstractGraphNode*> * mapPath;

	/** who does the user have selected */
	Agent * selected;
	AbstractGraphNode * selectedNode;

	DelaunySet * delauny;
	DelaunySet * delaunyEdit;
	TemplateVector<DelaunySet::VoronoiNode*> voronoiNodes;

	/** whether or not to draw debug lines for FSM steering behaviors */
	bool drawDebug;

	Agent * getAgentAt(V2f const & click);

	void gatherStaticObstaclesAt(CircF const & area, TemplateSet<Obstacle*> & out_obstacles);

	/** generate a list of agents in the given area */
	void gatherListOfAgentsAt(CircF const & area, TemplateVector<Agent*> & out_agents);
	
// helpful for finding memory by it's ID number (thanks _CrtSetDbgFlag!)
#define TRACE_MEMORY(mem, debugmessage) \
	printf("memID %d is %s\n", ((int*)mem)[-2], debugmessage)

	void generateWallBoxesForGraph(AbstractGraph * g, float wallEdgeValue);

#ifdef TESTING_SHAPES
	// used for testing object types
	BoxObject * testBox;
	CircleObject * testCircle;
	ConeObject * testcone;
	PolygonObject * testPoly;
	PolygonObject * testPoly2;
#endif

	Game();
	~Game();

	/** assumes that the graph will always have at least one node */
	GraphNode * getGraphNodeClosestTo(V2f position);

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
	bool raycast(V2f start, V2f direction, float maxDistance, bool dontCareAboutObstacle, 
		Obstacle * & out_obstacle, float & out_dist, V2f & out_point, V2f & out_normal, Obstacle ** ignoreList, int ignoreCount);

	void display(GLUTRenderingContext & g_screen);

	void gatherCollisions(Obstacle * s, TemplateSet<Obstacle*> & out_possibleObstacles);

	void update(int a_ms);
};
