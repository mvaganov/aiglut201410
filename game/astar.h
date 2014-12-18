#pragma once
#include "codegiraffe/templatevector.h"
#include "graph.h"
#include "codegiraffe/templatehashmap.h" // TemplateHashMap<K,V> is a replacement for std::map<K,V>

inline AbstractGraphNode * nodeInOpenSetWithLowestScore(
	TemplateVector<AbstractGraphNode*> & openset,
	TemplateHashMap<AbstractGraphNode*, float> & f_score)
{
	int lowestIndex = 0;
	AbstractGraphNode * best = openset[lowestIndex];
	float lowestF = f_score[best];
	for(int i = 1; i < openset.size(); ++i) {
		AbstractGraphNode * n = openset[i];
		if(f_score[n] < lowestF) {
			best = n;
			lowestIndex = i;
			lowestF = f_score[n];
		}
	}
	return best;
}

inline TemplateVector<AbstractGraphNode*> * reconstruct_path(
	TemplateHashMap<AbstractGraphNode*, AbstractGraphNode*> & came_from,
	AbstractGraphNode* current);

// the larger this is, the more direct A* tries to be.
static const float A_STAR_WEIGHT = 10; // 0 will cause a breadth-first search

inline float heuristic_cost_estimate(AbstractGraphNode * start, AbstractGraphNode * goal) {
	GraphNode * s = dynamic_cast<GraphNode*>(start);
	GraphNode * g = dynamic_cast<GraphNode*>(goal);
	if (s && g) return V2f::distance(s->getLocation(), g->getLocation()) * A_STAR_WEIGHT;
	return 0;
}

//function A*(start,goal)
inline TemplateVector<AbstractGraphNode*> * Astar(AbstractGraphNode * start, AbstractGraphNode * goal) {
    //closedset := the empty set    // The set of nodes already evaluated.
	TemplateVector<AbstractGraphNode*> closedset;
    //openset := {start}    // The set of tentative nodes to be evaluated, initially containing the start node
	TemplateVector<AbstractGraphNode*> openset;
	openset.add(start);
    //came_from := the empty map    // The map of navigated nodes.
	TemplateHashMap<AbstractGraphNode*, AbstractGraphNode*>	came_from;
	came_from[start] = NULL;
 
    //g_score[start] := 0    // Cost from start along best known path.
	TemplateHashMap<AbstractGraphNode*, float> g_score;
	g_score[start] = 0;
    // Estimated total cost from start to goal through y.
    //f_score[start] := g_score[start] + heuristic_cost_estimate(start, goal)
	TemplateHashMap<AbstractGraphNode*, float> f_score;
	f_score[start] = g_score[start] + heuristic_cost_estimate(start, goal);
 
    //while openset is not empty
	while(openset.size() > 0) {
		//for (int i = 0; i < closedset.size(); ++i) { GraphNode * n = closedset[i]; printf("%03x<%03x ", (((int)n >> 4) & 0xfff), (((int)came_from[n] >> 4) & 0xfff)); } printf("|");
		//for (int i = 0; i < openset.size(); ++i) { GraphNode * n = openset[i]; printf(" %03x<%03x", (((int)n >> 4) & 0xfff), (((int)came_from[n] >> 4) & 0xfff)); } printf("\n");

        //current := the node in openset having the lowest f_score[] value
		AbstractGraphNode * current = nodeInOpenSetWithLowestScore(openset, f_score);
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
		for(int i = 0; i < current->getEdgeCount(); ++i) {
			AbstractGraphNode * neighbor = current->getEdge(i)->otherNode(current);
            //if neighbor in closedset
			if(closedset.indexOf(neighbor) >= 0) {
                //continue
				continue;
			}
			// ignore negative cost edges
			float edgeCost = current->getEdge(i)->getCost();
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
inline TemplateVector<AbstractGraphNode*> * reconstruct_path(
	TemplateHashMap<AbstractGraphNode*, AbstractGraphNode*> & came_from,
	AbstractGraphNode* current)
{
    //total_path := [current]
	TemplateVector<AbstractGraphNode*> * total_path = new TemplateVector<AbstractGraphNode*>();
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
