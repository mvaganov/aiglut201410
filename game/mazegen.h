#pragma once

#include "graph.h"
#include "codegiraffe/random.h"
// from http://en.wikipedia.org/wiki/Maze_generation_algorithm#Randomized_Prim.27s_algorithm
inline void randomMazeGen_prim(Graph * g, float edgeWallValue) {
    // Start with a grid full of walls.
	GraphNode * n;
	for(int node = 0; node < g->nodes.size(); ++node) {
		n = g->nodes[node];
		for(int i = 0; i < n->edges.size(); ++i) { n->edges[i].cost = edgeWallValue; }
	}
    // Pick a cell, mark it as part of the maze. Add the walls of the cell to the wall list.
	TemplateVector<GraphNode*> theMaze;
	n = g->nodes[Random::PRNG() % g->nodes.size()];
	theMaze.add(n);
	TemplateVector<GraphEdge*> wallList;
	for(int i = 0; i < n->edges.size(); ++i) { wallList.add(&n->edges[i]); }
    // While there are walls in the list:
	while(wallList.size() > 0) {
        // Pick a random wall from the list. If the cell on the opposite side isn't in the maze yet:
		GraphEdge * randomWall = wallList[Random::PRNG() % wallList.size()];
		if(theMaze.indexOf(randomWall->to) == -1) {
            // Make the wall a passage and mark the cell on the opposite side as part of the maze.
			randomWall->calculateCost();
			randomWall->to->getEdgeTo(randomWall->from)->calculateCost();
			theMaze.add(randomWall->to);
            // Add the neighboring walls of the cell to the wall list.
			n = randomWall->to;
			for(int i = 0; i < n->edges.size(); ++i) { wallList.add(&n->edges[i]); }
		}
        // If the cell on the opposite side already was in the maze, remove the wall from the list.
		wallList.removeData(randomWall);
	}
}