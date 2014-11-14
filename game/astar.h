#pragma once
#include "codegiraffe/templatevector.h"
#include "graph.h"
#include <map>

inline GraphNode * nodeInOpenSetWithLowestScore(
	TemplateVector<GraphNode*> & openset, 
	std::map<GraphNode*, float> & f_score)
{
	int lowestIndex = 0;
	GraphNode * best = openset[lowestIndex];
	float lowestF = f_score[best];
	for(int i = 1; i < openset.size(); ++i) {
		GraphNode * n = openset[i];
		if(f_score[n] < lowestF) {
			best = n;
			lowestIndex = i;
			lowestF = f_score[n];
		}
	}
	return best;
}

inline TemplateVector<GraphNode*> * reconstruct_path(
	std::map<GraphNode *, GraphNode *> & came_from, 
	GraphNode* current);

// the larger this is, the more direct A* tries to be.
static const float A_STAR_WEIGHT = 10; // 0 will cause a breadth-first search

inline float heuristic_cost_estimate(GraphNode * start, GraphNode * goal) {
	return V2f::distance(start->location, goal->location) * A_STAR_WEIGHT;
}

//function A*(start,goal)
inline TemplateVector<GraphNode*> * Astar(GraphNode * start, GraphNode * goal) {
    //closedset := the empty set    // The set of nodes already evaluated.
	TemplateVector<GraphNode*> closedset;
    //openset := {start}    // The set of tentative nodes to be evaluated, initially containing the start node
	TemplateVector<GraphNode*> openset;
	openset.add(start);
    //came_from := the empty map    // The map of navigated nodes.
	std::map<GraphNode *, GraphNode *> came_from;
	came_from[start] = NULL;
 
    //g_score[start] := 0    // Cost from start along best known path.
	std::map<GraphNode*, float> g_score;
	g_score[start] = 0;
    // Estimated total cost from start to goal through y.
    //f_score[start] := g_score[start] + heuristic_cost_estimate(start, goal)
	std::map<GraphNode*, float> f_score;
	f_score[start] = g_score[start] + heuristic_cost_estimate(start, goal);
 
    //while openset is not empty
	while(openset.size() > 0) {
        //current := the node in openset having the lowest f_score[] value
		GraphNode * current = nodeInOpenSetWithLowestScore(openset, f_score);
        //if current = goal
		if(current == goal) {
            //return reconstruct_path(came_from, goal)
			return reconstruct_path(came_from, goal);
		}
        //remove current from openset
		openset.removeData(current);
        //add current to closedset
		closedset.add(current);
        //for each neighbor in neighbor_nodes(current)
		for(int i = 0; i < current->edges.size(); ++i) {
			GraphNode * neighbor = current->edges[i].to;
            //if neighbor in closedset
			if(closedset.indexOf(neighbor) >= 0) {
                //continue
				continue;
			}
			// ignore negative cost edges
			float edgeCost = current->edges[i].cost;
			if(edgeCost < 0) continue;

			//tentative_g_score := g_score[current] + dist_between(current,neighbor)
			float tentative_g_score = g_score[current] + edgeCost;
			
			//if neighbor not in openset or tentative_g_score < g_score[neighbor] 
			if(openset.indexOf(neighbor) < 0 || tentative_g_score < g_score[neighbor]) {
				//came_from[neighbor] := current
				came_from[neighbor] = current;
				//g_score[neighbor] := tentative_g_score
				g_score[neighbor] = tentative_g_score;
				//f_score[neighbor] := g_score[neighbor] + heuristic_cost_estimate(neighbor, goal)
				f_score[neighbor] = g_score[neighbor] + heuristic_cost_estimate(neighbor, goal);
				//if neighbor not in openset
				if(openset.indexOf(neighbor) < 0) {
					//add neighbor to openset
					openset.add(neighbor);
				}
			}
		}
	}
    //return failure
	return 0;
}
 
//function reconstruct_path(came_from,current)
inline TemplateVector<GraphNode*> * reconstruct_path(
	std::map<GraphNode *, GraphNode *> & came_from, 
	GraphNode* current)
{
    //total_path := [current]
	TemplateVector<GraphNode*> * total_path = new TemplateVector<GraphNode*>();
	total_path->add(current);
    //while current in came_from:
	while(came_from[current] != 0) {
        //current := came_from[current]
		current = came_from[current];
		//total_path.append(current)
		total_path->add(current);
	}
    //return total_path
	return total_path;
}
