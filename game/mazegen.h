#pragma once

#include "graph.h"
#include "codegiraffe/random.h"
#include "codegiraffe/templateset.h"
#include <stdio.h>
#include "delauny.h"

// from http://en.wikipedia.org/wiki/Maze_generation_algorithm#Randomized_Prim.27s_algorithm
inline void randomMazeGen_prim(AbstractGraph * g, float edgeWallValue, float minimumOpeningSize) {
    // Start with a grid full of walls.
	AbstractGraphNode * n;
	for(int node = 0; node < g->getNodeCount(); ++node) {
		n = g->getNode(Random::PRNG() % g->getNodeCount()); //n = &g->nodes[node];
		printf("%d(%d) ", node, n->getEdgeCount());
		for (int i = 0; i < n->getEdgeCount(); ++i) {
			n->getEdge(i)->setCost(edgeWallValue);
		}
	}
    // Pick a cell, mark it as part of the maze. Add the walls of the cell to the wall list.
	TemplateSet<AbstractGraphNode*> theMaze;
	n = g->getNode(Random::PRNG() % g->getNodeCount());// &g->nodes[Random::PRNG() % g->getNodeCount()];
	theMaze.add(n);
	TemplateSet<AbstractGraphEdge*> wallList;
	for (int i = 0; i < n->getEdgeCount(); ++i) { wallList.add(n->getEdge(i)); }
    // While there are walls in the list:
	while(wallList.size() > 0) {
        // Pick a random wall from the list. If the cell on the opposite side isn't in the maze yet:
		AbstractGraphEdge * randomWall = wallList[Random::PRNG() % wallList.size()];
		// if we can determine that the wall is too small to act like a corridor
		DelaunySet::Edge * dedge = dynamic_cast<DelaunySet::Edge*>(randomWall);
		if (dedge != NULL && dedge->getFace().radius < minimumOpeningSize) {
			// it cannot be an open corridor. remove it from te wallList, so it won't be turned into an opening
			wallList.removeData(randomWall);
			continue;
		}
		AbstractGraphNode * from = randomWall->getNode(0), *to = randomWall->getNode(1);
		if (theMaze.indexOf(from) == -1) { from = randomWall->getNode(1); to = randomWall->getNode(0); }
		if (theMaze.indexOf(to) == -1) {
            // Make the wall a passage and mark the cell on the opposite side as part of the maze.
			GraphEdge * e = dynamic_cast<GraphEdge*>(randomWall);
			if (e) e->calculateCost();
			AbstractGraphEdge * otherEdge = to->getEdgeTo(from);
			e = dynamic_cast<GraphEdge*>(otherEdge);
			if(e && e != randomWall) e->calculateCost();
			theMaze.add(to);
			// Add the neighboring walls of the cell to the wall list.
			n = to;
			for (int i = 0; i < n->getEdgeCount(); ++i) { wallList.add(n->getEdge(i)); }
		}
        // If the cell on the opposite side already was in the maze, remove the wall from the list.
		wallList.removeData(randomWall);
	}
}