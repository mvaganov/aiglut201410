#pragma once
#include <GL/freeglut.h>
#include "codegiraffe/v2.h"
#include "codegiraffe/templatevector.h"
#include "codegiraffe/box.h"
#include "codegiraffe/glutrenderingcontext.h"
class GraphNode;

struct GraphEdge {
	GraphNode * from;
	GraphNode * to;
	float cost;
	GraphEdge():from(0),to(0),cost(-1){}
	GraphEdge(GraphNode * f, GraphNode * t):from(f),to(t),cost(-1){}
	float calculateCost();
};

class GraphNode {
public:
	V2f location;
	TemplateVector<GraphEdge> edges;
	int getEdgeIndex(GraphNode * toThisNode) {
		for (int i = 0; i < edges.size(); ++i) {
			if (edges[i].to == toThisNode)
				return i;
		}
		return -1;
	}
	GraphEdge * getEdgeTo(GraphNode * n) {
		int index = getEdgeIndex(n);
		if (index != -1) return &edges[index];
		return 0;
	}
	void addNeighbor(GraphNode * n) {
		if (getEdgeIndex(n) == -1) {
			edges.add(GraphEdge(this, n));
		}
	}
	bool removeNeighbor(GraphNode * n) {
		int index = getEdgeIndex(n);
		if (index >= 0) {
			edges.remove(index);
			return true;
		}
		return false;
	}
	void glDraw(GLUTRenderingContext * g_screen) {
		BoxF(location, V2f(.5f, .5f), 0).glDraw();
		for(int i = 0; i < edges.size(); ++i) {
			if(edges[i].cost >= 0) {
				V2f a = edges[i].from->location, 
					b = edges[i].to->location;
				a.glDrawTo(b);
			}
		}
	}
};

inline float GraphEdge::calculateCost() { return cost = (to->location - from->location).magnitude(); }

class Graph { // Graph * g = Graph::createDenseGrid(5, 5, V2f(3,3), V2f(9,9));
public:
	TemplateVector<GraphNode*> nodes;
	static Graph * createDenseGrid(int rows, int cols, V2f min, V2f max) {
		float dx = (max.x-min.x) / (cols-1), dy = (max.y-min.y) / (rows-1);
		int totalNodes = rows * cols;
		Graph * g = new Graph;
		g->nodes.setSize(totalNodes);
		int index = 0;
		GraphNode * n;
		V2f cursor = min;
		for(int i = 0; i < g->nodes.size(); ++i) {
			g->nodes[i] = new GraphNode;
		}
		for(int row = 0; row < rows; ++row) {
			cursor.x = min.x;
			for(int col = 0; col < cols; ++col) {
				n = g->nodes[index];
				if(row != 0)		{ n->addNeighbor(g->nodes[index-cols]); } // not on bottom
				if(row != rows-1)	{ n->addNeighbor(g->nodes[index+cols]); } // not on top
				if(col != 0)		{ n->addNeighbor(g->nodes[index-1]); } // not on left
				if(col != cols-1)	{ n->addNeighbor(g->nodes[index+1]); } // not on right
				n->location = cursor;
				g->nodes[index++] = n;
				cursor.x += dx;
			}
			cursor.y += dy;
		}
		index = 0;
		return g;
	}
	void glDraw(GLUTRenderingContext * g_screen) {
		for(int i = 0; i < nodes.size(); ++i) {
			nodes[i]->glDraw(g_screen);
		}
	}
	~Graph() { for(int i = 0; i < nodes.size(); ++i) delete nodes[i]; }
};