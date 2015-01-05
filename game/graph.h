#pragma once
#include <GL/freeglut.h>
#include "codegiraffe/v2.h"
#include "codegiraffe/templatevector.h"
#include "codegiraffe/templatevectorlist.h"
#include "codegiraffe/box.h"
#include "codegiraffe/glutrenderingcontext.h"
#include <stdio.h>

class AbstractGraphNode;

class AbstractGraphEdge {
public:
	virtual AbstractGraphNode * getNode(const int n) const = 0;
	virtual AbstractGraphNode * otherNode(const AbstractGraphNode * n) const = 0;
	virtual float getCost() const = 0;
	virtual void setCost(const float cost) = 0;
	virtual ~AbstractGraphEdge(){}
	virtual bool has(const AbstractGraphNode * n) const = 0;
	virtual bool hasBoth(const AbstractGraphNode * n, const AbstractGraphNode * nn) const = 0;
	virtual bool operator==(AbstractGraphEdge const & e) const = 0;
};

class AbstractGraphNode {
public:
	virtual AbstractGraphEdge * getEdgeTo(const AbstractGraphNode * n) const = 0;
	virtual void addEdge(AbstractGraphEdge * e) = 0;
	virtual int getEdgeCount() const = 0;
	virtual AbstractGraphEdge * getEdge(const int index) = 0;
	virtual ~AbstractGraphNode(){}
};

class AbstractGraph {
public:
	virtual void connectNodes(AbstractGraphNode * from, AbstractGraphNode * to) = 0;
	virtual AbstractGraphEdge * getEdge(AbstractGraphNode * a, AbstractGraphNode * b, bool createIfNotThere) = 0;
	virtual AbstractGraphNode * createNode() = 0;
	virtual void destroyNode(AbstractGraphNode * n) = 0;
	virtual int getNodeCount() const = 0;
	virtual AbstractGraphNode * getNode(const int index) = 0;
	virtual ~AbstractGraph(){}
};

class GraphNode;
class Graph;

class GraphEdge : public AbstractGraphEdge {
protected:
	AbstractGraphNode * nodes[2];
	float cost;
public:
	void init() { nodes[0] = nodes[1] = NULL; cost = -1; }
	void set(AbstractGraphNode * a, AbstractGraphNode * b) { nodes[0] = a; nodes[1] = b; cost = -1; }
	GraphEdge() { init(); }
	GraphEdge(AbstractGraphNode * f, AbstractGraphNode * t) { set(f, t); }
	float calculateCost();
	/** if n is neither a nor b, returns NULL */
	AbstractGraphNode * otherNode(const AbstractGraphNode * n) const { return (n == nodes[0]) ? nodes[1] : ((n == nodes[1]) ? nodes[0] : NULL); }
	/** @param n 0 or 1 */
	AbstractGraphNode * getNode(const int n) const { return nodes[n]; }
	float getCost() const { return cost; }
	void setCost(const float cost) { this->cost = cost; }
	bool has(const AbstractGraphNode * n) const { return n == getNode(0) || n == getNode(1); }
	bool hasBoth(const AbstractGraphNode * n, const AbstractGraphNode * nn) const { return (n == getNode(0) && nn == getNode(1)) || (n == getNode(1) && nn == getNode(0)); }
	bool operator==(AbstractGraphEdge const & e) const { return hasBoth(e.getNode(0), e.getNode(1)); }
};

class GraphNode : public AbstractGraphNode {
protected:
	TemplateVector<AbstractGraphEdge*> edges;
	V2f location;
public:
	~GraphNode(){}
	const V2f & getLocation() const { return location; }
	void setLocation(const V2f loc) { location = loc; }
	int getEdgeCount() const { return edges.size(); }
	AbstractGraphEdge * getEdge(const int index) { return edges[index]; }
	int getEdgeIndex(const AbstractGraphNode * toThisNode) const {
		for (int i = 0; i < edges.size(); ++i) {
			if (edges[i]->otherNode(this) == toThisNode)
				return i;
		}
		return -1;
	}
	void addEdge(AbstractGraphEdge * e) { edges.add(e); }
	AbstractGraphEdge * getEdgeTo(const AbstractGraphNode * n) const {
		int index = getEdgeIndex(n);
		if (index != -1) return edges[index];
		return 0;
	}
	bool removeNeighbor(AbstractGraphNode * n) {
		int index = getEdgeIndex(n);
		if (index >= 0) {
			edges.remove(index);
			return true;
		}
		return false;
	}
	void glDraw(GLUTRenderingContext * g_screen) {
		Boxf(location, V2f(.5f, .5f), 0).glDraw();
		for(int i = 0; i < edges.size(); ++i) {
			if(edges[i]->getCost() >= 0) {
				V2f a = ((GraphNode*)edges[i]->getNode(0))->location,
					b = ((GraphNode*)edges[i]->getNode(1))->location;
				a.glDrawTo(b);
			}
		}
	}
};

inline float GraphEdge::calculateCost() {
	GraphNode * a = (GraphNode*)nodes[0];
	GraphNode * b = (GraphNode*)nodes[1];
	return cost = (b->getLocation() - a->getLocation()).magnitude();
}

class Graph : public AbstractGraph { // Graph * g = Graph::createDenseGrid(5, 5, V2f(3,3), V2f(9,9));
public:
	/** memory pool for nodes */
	TemplateVectorList<GraphNode> nodes;
	/** memory pool for edges */
	TemplateVectorList<GraphEdge> edges;
	bool edgesSameBothWays;

	AbstractGraphNode * createNode() {
		nodes.add();
		return &nodes.getLast();
	}
	int getNodeCount() const { return nodes.size(); }
	AbstractGraphNode * getNode(const int index) { return &nodes[index]; }

	void destroyNode(AbstractGraphNode * n) { }

	Graph(bool edgesSameBothWays) :edgesSameBothWays(edgesSameBothWays){}
	~Graph() {}

	AbstractGraphEdge * getEdge(AbstractGraphNode * a, AbstractGraphNode * b, bool createIfNotThere) {
		for (int i = 0; i < edges.size(); ++i) {
			if (                     (edges[i].getNode(0) == a && edges[i].getNode(1) == b)
			|| (edgesSameBothWays && (edges[i].getNode(0) == b && edges[i].getNode(1) == a)))
				return &edges[i];
		}
		if (createIfNotThere) {
			edges.add();
			GraphEdge * e = &edges.getLast();
			e->set(a, b);
			return e;
		}
		return NULL;
	}

	void connectNodes(AbstractGraphNode * a_from, AbstractGraphNode * a_to) {
		GraphNode * from = (GraphNode*)a_from;
		GraphNode * to = (GraphNode*)a_to;
		if (from->getEdgeIndex(to) == -1) {
			AbstractGraphEdge * e = NULL;
			if (edgesSameBothWays && (e = to->getEdgeTo(from)) != NULL) {
				from->addEdge(e);
			} else {
				e = getEdge(from, to, true);
				from->addEdge(e);
				if (edgesSameBothWays) { to->addEdge(e); }
			}
		}
	}

	static Graph * createDenseGrid(int rows, int cols, V2f min, V2f max) {
		float dx = (max.x-min.x) / (cols-1), dy = (max.y-min.y) / (rows-1);
		int totalNodes = rows * cols;
		Graph * g = new Graph(true);
		g->nodes.setSize(totalNodes);
		int index = 0;
		GraphNode * n;
		V2f cursor = min;
		for(int row = 0; row < rows; ++row) {
			cursor.x = min.x;
			for(int col = 0; col < cols; ++col) {
				n = &g->nodes[index];
				if (row != 0)		{ g->connectNodes(n, &g->nodes[index - cols]); printf("%d->%d ", index, index-cols); } // not on bottom
				if (row != rows - 1){ g->connectNodes(n, &g->nodes[index + cols]); printf("%d->%d ", index, index + cols); } // not on top
				if (col != 0)		{ g->connectNodes(n, &g->nodes[index - 1]); printf("%d->%d ", index, index - 1); } // not on left
				if (col != cols - 1){ g->connectNodes(n, &g->nodes[index + 1]); printf("%d->%d ", index, index + 1); } // not on right
				n->setLocation(cursor);
				index++;
				cursor.x += dx;
			}
			cursor.y += dy;
		}
		index = 0;
		return g;
	}
	void glDraw(GLUTRenderingContext * g_screen) {
		for(int i = 0; i < nodes.size(); ++i) {
			nodes[i].glDraw(g_screen);
		}
	}
};