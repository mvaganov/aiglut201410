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
#define TESTING_SHAPES

class Game {
public:
	GLUTRenderingContext * g_screen;
	V2f mousePosition, mouseDragged;
	Ray userRay;
	float userRayDistance;

	float rotationLast, rotationCurrent, rotationDelta;
	V2f operationPoint;

	AABBf selectedArea;
	TemplateSet<Obstacle*> selectedObstacles;
	TemplateSet<Obstacle*> selectedObstacleSingle;
	V2f selectedOrigin;

	void refreshSelection();

	Obstacle* getSelectedAtPosition(V2f position, long mask);

	void addObstacle(Obstacle * obstacle);

	TemplateVector<Obstacle*> obstacles;
//	TemplateVector<Agent*> agents;

//	CellSpacePartition * staticObstaclesMap;
//	CellSpacePartition * movingObstaclesMap;
	LayeredPartitions objectMap;

	//Graph * mapGraph;
	//TemplateVector<AbstractGraphNode*> * mapPath;

	///** who does the user have selected */
	//Agent * selected;
	//AbstractGraphNode * selectedNode;
	//Agent * getAgentAt(V2f const & click);
	///** generate a list of agents in the given area */
	//void gatherListOfAgentsAt(Circf const & area, TemplateVector<Agent*> & out_agents);


	//DelaunySet * delauny;
	//DelaunySet * delaunyEdit;
	//TemplateVector<DelaunySet::VoronoiNode*> voronoiNodes;

	/** whether or not to draw debug lines for FSM steering behaviors */
	bool drawDebug;

	bool gameLoopRunning;
	bool paused;
	bool drawingSelectionBox;
	bool draggingSelected;

	bool doingOperation;
	bool doingRotation, doingTranslation, doingScale;

//	void gatherStaticObstaclesAt(Circf const & area, TemplateSet<Obstacle*> & out_obstacles);

//	void gatherObstaclesAt(Circf const & area, TemplateSet<Obstacle*> & out_obstacles);

	
// helpful for finding memory by it's ID number (thanks _CrtSetDbgFlag!)
#define TRACE_MEMORY(mem, debugmessage) \
	printf("memID %d is %s\n", ((int*)mem)[-2], debugmessage)

	void generateWallBoxesForGraph(AbstractGraph * g, float wallEdgeValue);

	Game();
	void init();
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
	 * @return true if nothing was hit
	 */
	bool raycast(Ray ray, RaycastHit & out_rh, float maxDistance, bool dontCareAboutObstacle, 
		Obstacle * & out_obstacle, long mask);

	void draw(GLUTRenderingContext * g);

	void gatherCollisions(Obstacle * s, TemplateSet<Obstacle*> & out_possibleObstacles);

	void update(int a_ms);

	/**
	 * @param button GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, GLUT_RIGHT_BUTTON
	 * @param state 0: pressed, 1: released, 2: dragged, -1:not actually active
	 */
	void mouse(V2f position, int button, int state);

	/**
	* @param button keycode
	* @param state 0: pressed, 1: released, 2: dragged, -1:not actually active
	*/
	void key(V2f position, int button, int state);
};
