#pragma once

#include "delauny.h"
#include <stdio.h>

DelaunySet::VoronoiFace::VoronoiFace(){}

void DelaunySet::VoronoiFace::refresh(Edge * e) {
	GraphNode * a = dynamic_cast<GraphNode*>(e->getNode(0));
	GraphNode * b = dynamic_cast<GraphNode*>(e->getNode(1));
	normal = a->getLocation() - b->getLocation();
	normal.normalize();
	center.setZero();
	points.setSize(e->getNeighborTriangulationCount());
	for (int i = 0; i < e->getNeighborTriangulationCount(); ++i) {
		points[i] = e->getNeighborTriangulation(i)->circum.center;
		center += points[i];
	}
	center /= (float)points.size();
	V2f va, vb;
	gatherOppositePoints(va, vb);
	radius = (vb - va).magnitude();
}
void DelaunySet::VoronoiFace::draw(GLUTRenderingContext * g) const {
	for (int i = 0; i < points.size(); ++i) {
		g->drawLine(center, points[i]);
		//center.glDrawTo(points[i]);
	}
	g->drawLine(center, center + normal);
	//center.glDrawTo(center + normal);
}
/** positive value on one side, negative on the other. not sure which is which. */
float DelaunySet::VoronoiFace::sideValue(V2f const & p) { return V2f::dot(normal, p - center); }

void DelaunySet::VoronoiFace::gatherOppositePoints(V2f & out_a, V2f & out_b) {
	out_a = center;
	out_b = center;
	if (points.size() == 1) {
		out_b = points[0];
		return;
	}
	float dist, maxDist = -1;
	for (int a = 0; a < points.size(); ++a) {
		for (int b = a + 1; b < points.size(); ++b) {
			dist = V2f::distance(points[a], points[b]);
			if (maxDist < 0 || dist > maxDist) {
				out_a = points[a];
				out_b = points[b];
				maxDist = dist;
			}
		}
	}
}

DelaunySet::VoronoiFace & DelaunySet::Edge::getFace() { return voronoiFace; }

int DelaunySet::Edge::getNeighborTriangulationCount() const { return neighborTri.size(); }

DelaunySet::Triangulation* DelaunySet::Edge::getNeighborTriangulation(const int i) const { return neighborTri[i]; }

void DelaunySet::Edge::refreshVoronoiFace() { voronoiFace.refresh(this); }

void DelaunySet::Edge::addNeighborTriangulation(Triangulation * nt) {
	bool added = neighborTri.add(nt);
	if (added) refreshVoronoiFace();
}
bool DelaunySet::Edge::removeNeighborTriangulation(Triangulation * nt) {
	bool removed = neighborTri.removeData(nt);
	if (removed) refreshVoronoiFace();
	return removed;
}

void DelaunySet::Edge::invalidateEdge() {
	if (getNode(0)) {
		VoronoiNode * vn = dynamic_cast<VoronoiNode*>(getNode(0));
		if (!vn->removeEdge(this)) {
			printf("could not remove edge from node?\n"); ABORT
		}
	}
	if (getNode(1)) {
		VoronoiNode * vn = dynamic_cast<VoronoiNode*>(getNode(1));
		if (!vn->removeEdge(this)) {
			printf("could not remove edge from node?\n"); ABORT
		}
	}
	nodes[0] = nodes[1] = NULL;
}
bool DelaunySet::Edge::isValid() const { return getNode(0) != NULL && getNode(1) != NULL; }
DelaunySet::Edge::Edge() { nodes[0] = nodes[1] = NULL; invalidateEdge(); }
DelaunySet::Edge::Edge(VoronoiNode *a, VoronoiNode *b) : GraphEdge(a, b), shape(&a->getLocation(), &b->getLocation()){}// { nodes[0] = a; nodes[1] = b; }
void DelaunySet::Edge::draw(GLUTRenderingContext * g) const {
	GraphNode * a = (GraphNode*)getNode(0), *b = (GraphNode*)getNode(1);
	a->getLocation().glDrawTo(b->getLocation());
}
bool DelaunySet::Edge::isNeighbor(Triangulation * const t) const { return neighborTri.indexOf(t) >= 0; }

void DelaunySet::Triangulation::drawEdges(GLUTRenderingContext * g) const { for (int i = 0; i < edges.size(); ++i) { edges[i]->draw(g); } }

DelaunySet::Triangulation & DelaunySet::Triangulation::operator=(Triangulation const & t) {
	circum = t.circum;
	edges = t.edges;
	startingNode = t.startingNode;
	centerMass = t.centerMass;
	nodeSet = t.nodeSet;
	return *this;
}

void DelaunySet::Triangulation::set(Edge* const * list, const int listCount, Circf const & circumscription) {
	startingNode = NULL;
	circum = circumscription;
	edges.allocateToSize(listCount);
	for (int i = 0; i < listCount; ++i) { edges[i] = list[i]; }
	calculate();
}

DelaunySet::Triangulation::Triangulation(Edge* const * list, const int listCount, Circf const & circumscription) : startingNode(NULL), mask(Obstacle::DELAUNY_TRIANGULATION) {
	set(list, listCount, circumscription);
}

/** create an invalid triangulation by default */
DelaunySet::Triangulation::Triangulation() : startingNode(NULL), circum(Circf(0, 0, -1)), mask(Obstacle::DELAUNY_TRIANGULATION) { }

int DelaunySet::Triangulation::getEdgeCount() const { return edges.size(); }
DelaunySet::Edge * DelaunySet::Triangulation::getEdge(const int i) const { return edges[i]; }

/**
* calculate center of mass, and order edges in clockwise fashion.
* @return true if this is a valid Triangulation
*/
bool DelaunySet::Triangulation::calculate() {
	gatherNodesInto(nodeSet);
	TemplateArray<VoronoiNode*> inOrder;
	inOrder.allocateToSize(nodeSet.size());
	orderNodesClockwise(nodeSet.getRawList(), nodeSet.size(), centerMass, inOrder.getRawList());
	startingNode = inOrder[0];
	// now order edges to match inOrder
	for (int i = 0; i < edges.size(); ++i) {
		VoronoiNode * n0 = inOrder[i];
		VoronoiNode * n1 = inOrder[(i + inOrder.size() - 1) % inOrder.size()];
		if (n0 == n1) { continue; }
		// find the edge that should be here
		int goesHere = this->findEdgeIndex(n0, n1);
		// swap it with what is in here
		if (i != goesHere)
			edges.swap(i, goesHere);
	}
	return true;
}

/** @return if this is a valid triangulation */
bool DelaunySet::Triangulation::isValidTriangle() const { return circum.radius >= 0; }

/** make this an invalid triangulation. put any edges freed (nolonger referenced by any triangles) into the given set */
void DelaunySet::Triangulation::invalidate(MemPool<Edge> & edgePool) {
	circum.radius = -1;
	startingNode = NULL;
	for (int i = edges.size() - 1; i >= 0; --i) {
		Edge * e = edges[i];
		if (!e->removeNeighborTriangulation(this)) {
#ifndef NO_TEST
			printf("uhhh... could not remove this triangle, it's not a neighbor...\n"); ABORT
#endif
		}
		if (e->getNeighborTriangulationCount() == 0) {
			e->invalidateEdge();
			edgePool.markFree(e);
		}
		edges[i] = NULL; // don't de-allocate the edges, this triangulation will probably be used again, and might need the same number of edges
	}
}

/** put all of the nodes from referenced edges into the given set */
void DelaunySet::Triangulation::gatherNodesInto(TemplateSet<VoronoiNode*> & nodes) const {
	for (int i = 0; i < edges.size(); ++i) {
		nodes.add((VoronoiNode*)edges[i]->getNode(0));
		nodes.add((VoronoiNode*)edges[i]->getNode(1));
	}
}

const TemplateSet<DelaunySet::VoronoiNode*> & DelaunySet::Triangulation::getNodeSet() const { return nodeSet; }

/** only works for 2D triangulation */
void DelaunySet::Triangulation::gatherPointsClockwise(TemplateArray<V2f> & points) const {
	points.allocateToSize(edges.size());
	int index = 0;
	points[index++] = startingNode->getLocation();
	VoronoiNode * cursor = startingNode;
	while (index < edges.size()) {
		cursor = (VoronoiNode*)edges[index]->otherNode(cursor);
		points[index] = cursor->getLocation();
		index++;
	}
}

/** only works for 2D triangulation */
bool DelaunySet::Triangulation::polygonContains(V2f const & p) const {
	if (!this->isValidTriangle()) return false;
	if (p.distance(circum.center) > circum.radius) return false;
	TemplateArray<V2f> points;
	gatherPointsClockwise(points);
	return p.isInsidePolyCW(points.getRawList(), points.size());
}

/** @return true if this triangulation is using the given edge */
bool DelaunySet::Triangulation::hasEdge(const Edge * const edge) const {
	for (int i = 0; i < edges.size(); ++i) {
		if (edges[i] == edge) return true;
	} return false;
}

/** @return true if this triangle referenced the given triangle as a neighbor */
bool DelaunySet::Triangulation::isNeighborsWith(Triangulation * t) const {
	for (int i = 0; i < edges.size(); ++i) {
		if (edges[i]->isNeighbor(t)) return true;
	}
	return false;
}

/** get the edge that connects the two given nodes. useful when the list of nodes is known but not the edges */
int DelaunySet::Triangulation::findEdgeIndex(VoronoiNode * n0, VoronoiNode * n1) const {
	for (int i = 0; i < edges.size(); ++i) {
		if (edges[i]->hasBoth(n0, n1)) return i;
	}
	return -1;
}

/** called after a Triangulation is created, to ensure that it's edges are aware of it */
void DelaunySet::Triangulation::recognizeNeighbors() {
	for (int i = 0; i < edges.size(); ++i) {
		edges[i]->addNeighborTriangulation(this);
	}
}

/** check this triangulation's circumscription */
bool DelaunySet::Triangulation::circumscriptionContains(V2f const & p) const { return circum.contains(p); }

/** @return if this Triangulation is related to the given node */
bool DelaunySet::Triangulation::hasNode(const VoronoiNode * node) const {
	for (int i = 0; i < edges.size(); ++i) { if (edges[i]->has(node)) return true; }
	return false;
}
/** @return if this Triangulation has every node in the given series */
bool DelaunySet::Triangulation::hasAllNodes(VoronoiNode ** list, const int listSize) const {
	for (int i = 0; i < listSize; ++i) {
		if (!hasNode(list[i])) return false;
	}
	return true;
}
/** @return if this node has exactly the nodes in this given series */
bool DelaunySet::Triangulation::hasExactlyTheseNodes(VoronoiNode ** list, const int listSize) const {
	const TemplateSet<VoronoiNode*> & nodes = getNodeSet();
	//gatherNodesInto(nodes);
	if (nodes.size() != listSize) return false;
	for (int i = 0; i < listSize; ++i) {
		if (nodes.indexOf(list[i]) == -1) return false;
	}
	return true;
}

/** a triangulation is equal if it is brokering between the exact same nodes */
bool DelaunySet::Triangulation::equals(const Triangulation & t) const {
	if (&t == this) return true;
	return getNodeSet() == t.getNodeSet();
}
bool DelaunySet::Triangulation::operator==(Triangulation const & t) const { return equals(t); }
bool DelaunySet::Triangulation::operator!=(Triangulation const & t) const { return !equals(t); }


DelaunySet::VoronoiNode::VoronoiNode() : atBoundaryOf(NULL), valid(false), needsModelRecalculated(true), mask(Obstacle::VORONOI_NODE) {}
DelaunySet::VoronoiNode::VoronoiNode(V2f const & p) : mask(Obstacle::VORONOI_NODE) { set(p); }

V2f & DelaunySet::VoronoiNode::getLocation() { return location; }
const V2f & DelaunySet::VoronoiNode::getLocation() const { return location; }

/** @return whether the polygon needs to be recalculated */
bool DelaunySet::VoronoiNode::isDirty() { return needsModelRecalculated; }
/** identify that the polygon needs to be recalculated */
void DelaunySet::VoronoiNode::setDirty() { needsModelRecalculated = true; }

bool DelaunySet::VoronoiNode::isBorderPolygon() const { return atBoundaryOf != NULL; }

int DelaunySet::VoronoiNode::getEdgeCount() const { return edges.size(); }
//		Edge * getEdge(const int i) const { return edges[i]; }
bool DelaunySet::VoronoiNode::removeEdge(Edge * e) { setDirty(); return edges.removeData(e); }
//		bool addEdge(Edge * e) { setDirty(); return edges.add(e); }

void DelaunySet::VoronoiNode::set(V2f const & p) { point.set(p); location = p; valid = true; setDirty(); atBoundaryOf = NULL; }

bool DelaunySet::VoronoiNode::isValidNode() const { return valid; }

void DelaunySet::VoronoiNode::invalidate(DelaunySet & D, TemplateSet<VoronoiNode*> & changedNodes) {
	while (edges.size() > 0) {
		Edge * e = dynamic_cast<Edge*>(edges.getLast());
		int sizeAtStart = edges.size();
		for (int i = e->getNeighborTriangulationCount() - 1; i >= 0; --i) {
			Triangulation * t = e->getNeighborTriangulation(i);
			t->gatherNodesInto(changedNodes);
			t->invalidate(D.edgePool);
			D.triangulationPool.markFree(t);//D.freeTriangles.add(t);
		}
		int sizeAtEnd = edges.size();
		if (sizeAtStart == sizeAtEnd)
			edges.pop();
	}
	valid = false;
	setDirty();
}

/** calculate the Voronoi polygon*/
void DelaunySet::VoronoiNode::calculatePolygon(TemplateVector<V2f> & out_points, V2f & out_center) const {
	TemplateSet<Triangulation*> triangles;
	TemplateSet<VoronoiNode*> nodes;
	gatherTriangles(triangles);
	for (int i = 0; i < triangles.size(); ++i) {
		out_points.add(triangles[i]->circum.center);
	}
	Polygon2f::calculatePolygonCW(out_points.getRawList(), out_points.size(), out_center);
}

/** calculate the Voronoi polygon. TODO make boundary a Shape rather than an obstacle. */
Polygon2f & DelaunySet::VoronoiNode::getPolygon2f(Obstacle * boundary) {
	if (needsModelRecalculated) {
		atBoundaryOf = NULL;
		//printf("-------------------\n");
		//TemplateVector<LineF> lines;
		calculatedPolygon.points.clear();
		TemplateSet<VoronoiNode*> nodes;
		// get all the triangulations that this node is bordering
		TemplateSet<Triangulation*> triangles;
		gatherTriangles(triangles);
		// to make the polygon, we need the centers of each triangle bordering this node
		TemplateVector<V2f> points;
		TemplateSet<Polygon2f::pair> pairs;
		points.setSize(triangles.size());
		for (int i = 0; i < triangles.size(); ++i) {
			points[i] = triangles[i]->circum.center;
		}
		//for (int i = 0; i < points.size(); ++i) { printf("(%.1f %.1f)", points[i].x, points[i].y); } printf("\n");
		// for each triangle
		TemplateSet<Triangulation*> neighborTriangles;
		for (int triIndex = 0; triIndex < triangles.size(); ++triIndex) {
			Triangulation * t = triangles[triIndex];
			neighborTriangles.clear();
			// go through the edges
			for (int n = 0; n < t->getEdgeCount(); ++n) {
				Edge * e = t->getEdge(n);
				// and see what other triangles border both this triangle (obiovusly), and this node.
				if (e->has(this)) {
					for (int i = 0; i < e->getNeighborTriangulationCount(); ++i) {
						neighborTriangles.add(e->getNeighborTriangulation(i));
					}
				}
				// if this triangle borders the entire system (it doesn't have enough triangles to fully surround the edge)
				if (e->getNeighborTriangulationCount() == 1 && e->has(this)) {
					// find the dividing plane
					V2f rayDir = e->getFace().normal.perp(), end, otherPoint;//-(delta.perp().normal()), end;
					GraphNode * a = (GraphNode*)e->getNode(0);
					GraphNode * b = (GraphNode*)e->getNode(1);
					V2f::closestPointOnLine(a->getLocation(), b->getLocation(), t->centerMass, otherPoint);
					V2f correctDirection = otherPoint - t->centerMass;
					float alignment = V2f::dot(rayDir, correctDirection);
					if (alignment < 0) { rayDir *= -1.0f; }
					// raycast the dividing plane to the border of the boundary
					if (!boundary) {
						end = t->circum.center + rayDir * 100;
					}
					else if (boundary->getShape()->contains(t->circum.center)) {
						RaycastHit rh;
//						V2f hitNorm;
//						float dist;
						if (boundary->getShape()->raycast(Ray(t->circum.center, rayDir), rh)) {//dist, end, hitNorm)) {
							end = rh.point;
							//printf("should have hit from %.1f %.1f\n");
							// add that location to the list of points (noting it's index)
							points.add(end);
							//printf("A %.1f %.1f\n", end.x, end.y);
							int index = points.size() - 1;
							// add a line to that boundary point
							pairs.add(Polygon2f::pair(triIndex, index));
						}
						else {
							end = t->circum.center;
							//printf("should have hit from %.1f %.1f\n");
						}
					}
					atBoundaryOf = boundary;
				}
			}
			// draw lines from this triangle to all other node-adjacent triangles
			for (int n = 0; n < neighborTriangles.size(); ++n) {
				int otherTriIndex = triangles.indexOf(neighborTriangles[n]);
				// a single point cannot make a line...
				if (triIndex == otherTriIndex) continue;
				if (boundary != NULL) {
					V2f thisPoint = triangles[triIndex]->circum.center;
					V2f thatPoint = triangles[otherTriIndex]->circum.center;
					bool thisInBoundary = boundary->getShape()->contains(thisPoint);
					bool thatInBoundary = boundary->getShape()->contains(thatPoint);
					if (!thisInBoundary && !thatInBoundary) continue;
					if (!thisInBoundary || !thatInBoundary) {
						atBoundaryOf = boundary;
					}
				}
				pairs.add(Polygon2f::pair(triIndex, otherTriIndex));

			}
		}
		// to make the list of points relative to the node center (since polygon models are relative to a center)
		//for (int i = 0; i < points.size(); ++i) { printf("(%.1f %.1f)", points[i].x, points[i].y); } printf("\n");
		for (int i = 0; i < points.size(); ++i) {
			points[i] -= this->location;
			//printf("(%.1f %.1f)", points[i].x, points[i].y);
		}
		//printf("\n");
		// that should be all the data needed for a polygon
		calculatedPolygon.set(this->location, V2f::ZERO_DEGREES(), points.getRawListConst(), points.size(), pairs.getRawListConst(), pairs.size(), true);
		//for (int i = 0; i < polygon.points.size(); ++i) { printf("(%.1f %.1f)", polygon.points[i].x, polygon.points[i].y); } printf("\n");
		needsModelRecalculated = false;
	}
	return calculatedPolygon;
}

bool DelaunySet::VoronoiNode::polyhedronContains(V2f const & a_position) const {
	for (int i = 0; i < edges.size(); ++i) {
		Edge * e = (Edge*)edges[i];
		float sideOfCenter = e->getFace().sideValue(location);
		float sideOfPosition = e->getFace().sideValue(a_position);
		if ((sideOfCenter < 0) != (sideOfPosition < 0)) {
			return false;
		}
	}
	return true;
}

/** @return the edge from this node to the given node */
DelaunySet::Edge * DelaunySet::VoronoiNode::getEdgeTo(DelaunySet::VoronoiNode * node) const {
	Edge * e;
	for (int i = 0; i < edges.size(); ++i) {
		e = (Edge*)edges[i];
		if (e->has(node))
		{
#ifndef NO_TEST
			if (!e->hasBoth(node, this)) {
				printf("0x%08x 0x%08x 0x%08x 0x%08x", node, this, e->getNode(0), e->getNode(1));
				printf("weird... this edge doesn't include *this"); ABORT
			}
#endif
			return e;
		}
	}
	return NULL;
}

/** gather a list of every triangulation attached to this node */
void DelaunySet::VoronoiNode::gatherTriangles(TemplateSet<Triangulation*> & triangles) const {
	for (int i = 0; i < edges.size(); ++i) {
		Edge * e = (Edge*)edges[i];
		for (int n = 0; n < e->getNeighborTriangulationCount(); ++n) {
			triangles.add(e->getNeighborTriangulation(n));
		}
	}
}

void DelaunySet::VoronoiNode::gatherEdges(TemplateSet<Edge*> & edges) const {
	for (int i = 0; i < edges.size(); ++i) {
		edges.add(edges[i]);
	}
}

/**
* @param nodes what nodes to calculate the clockwise order of
* @param nodeCount how many nodes there are
* @param out_centerMass the average location of the nodes
* @param out_inOrder where to put the list re-arranged in order
*/
void DelaunySet::orderNodesClockwise(VoronoiNode*const* nodes, const int nodeCount, V2f  & out_centerMass, VoronoiNode** out_inOrder) {
	TemplateArray<V2f > points;
	points.allocateToSize(nodeCount);
	for (int i = 0; i < points.size(); ++i){
		points[i] = nodes[i]->getLocation();
	}
	Polygon2f::calculatePolygonCW(points.getRawList(), points.size(), out_centerMass);
	for (int index = 0; index < nodeCount; ++index) {
		out_inOrder[index] = 0;
		for (int i = 0; i < points.size(); ++i) {
			if (points[index] == nodes[i]->getLocation()) {
				out_inOrder[index] = nodes[i];
			}
		}
#ifndef NO_TEST
		if (out_inOrder[index] == 0) { ABORT }
#endif
	}
}

int DelaunySet::getNodeCount() const { return m_nodePool.countUsed(); }//currentNodes.size(); }
AbstractGraphNode * DelaunySet::getNode(const int index) {
	return &m_nodePool.get(index);
}//currentNodes[index]; }

void DelaunySet::connectNodes(AbstractGraphNode * from, AbstractGraphNode * to) {
	getEdge(from, to, true);
}
AbstractGraphEdge * DelaunySet::getEdge(AbstractGraphNode * a, AbstractGraphNode * b, bool createIfNotThere) {
	VoronoiNode* va = (VoronoiNode*)a, *vb = (VoronoiNode*)b;
	Edge * e = findEdge(va, vb);
	if (e == NULL && createIfNotThere) {
		e = addEdge(Edge(va, vb));
	}
	a->addEdge(e);
	b->addEdge(e);
	return e;
}

DelaunySet::DelaunySet(Obstacle * boundary) :selectedNode(0), boundary(boundary), csp(0){
	V2f c = boundary->getShape()->getCenter();
	float rad = boundary->getShape()->getRadius();
	V2f r(rad, rad);
	AABBf area(c - r, c + r);
	csp = new CellSpacePartition(area, V2i(4, 4));
}

DelaunySet::~DelaunySet() { 
	if (csp) {
		delete csp;
		csp = 0;
	}
}

DelaunySet::Triangulation * DelaunySet::addTriangle(Triangulation const & t) {
	//Triangulation * tri;
	//if (freeTriangles.size() > 0) {
	//	tri = freeTriangles.pop();
	//}
	//else {
	//	currentTriangles.add();
	//	tri = &currentTriangles.getLast();
	//}
	//*tri = t;
	//return tri;
	return &triangulationPool.add(t);
}

void DelaunySet::removeTriangle(Triangulation * toRemove) {
	toRemove->invalidate(edgePool);
	triangulationPool.markFree(toRemove);//freeTriangles.add(toRemove);
	csp->remove(toRemove);
#ifndef NO_TEST
	// check to see if anything references the bad triangle still.
	for (int i = 0; i < triangulationPool.size(); ++i) {
		Triangulation * t = &triangulationPool[i];
		if (t->isValidTriangle() && t != toRemove) {
			if (t->isNeighborsWith(toRemove)) {
				printf("bad triangle neighbor.\n"); ABORT
			}
		}
	}
#endif
}

void DelaunySet::destroyNode(AbstractGraphNode * n) {
	removeNode((VoronoiNode*)n);
}

void DelaunySet::removeNode(VoronoiNode * toRemove) {
	TemplateSet<VoronoiNode*> changedNodes;
	removeNode(toRemove, changedNodes);
	calculateTriangles(changedNodes.getRawList(), changedNodes.size());
}

void DelaunySet::removeNode(VoronoiNode * toRemove, TemplateSet<VoronoiNode*> & changedNodes) {
	toRemove->invalidate(*this, changedNodes);
	csp->remove(toRemove);
	m_nodePool.markFree(toRemove);//freeNodes.add(toRemove);
	changedNodes.removeData(toRemove);
}

void DelaunySet::moveNode(VoronoiNode * toMove, V2f newPosition) {
	TemplateSet<VoronoiNode*> changedNodes;
	removeNode(toMove, changedNodes);
	addNode(newPosition, changedNodes);
	calculateTriangles(changedNodes.getRawList(), changedNodes.size());
}

DelaunySet::Edge * DelaunySet::findEdge(VoronoiNode * a, VoronoiNode * b) {
	//DelaunySet::Edge * eA = a->getEdgeTo(b);
	//DelaunySet::Edge * eB = a->getEdgeTo(b);
	//if (eA != eB) {
	//	ABORT
	//}
	//if (eA != NULL)
	//	return eA;
	for (int i = 0; i < edgePool.size(); ++i) {
		if (edgePool[i].hasBoth(a, b)) return &edgePool[i];
	}
	return NULL;
//	return a->getEdgeTo(b);
}

DelaunySet::Edge* DelaunySet::addEdge(Edge const & e) {
	//if (freeEdges.size() > 0) {
	//	Edge * ed = freeEdges.pop();
	//	*ed = e;
	//	return ed;
	//}
	return &edgePool.add(e); //return &edgePool.getLast();
}

DelaunySet::Edge* DelaunySet::marshalEdge(VoronoiNode * a, VoronoiNode * b) {
	Edge * e = findEdge(a, b);
	if (e == NULL) {
		e = addEdge(Edge(a, b));
	}
	a->addEdge(e);
	b->addEdge(e);
	return e;
}

DelaunySet::Triangulation * DelaunySet::getTriangle(TemplateSet<VoronoiNode*> & nodeCluster) {
	VoronoiNode * nod;
	Edge * ed;
	Triangulation *foundTriangulation = NULL, *tri;
	bool foundInHere = false;
	// for each node
	for (int n = 0; n < nodeCluster.size(); ++n) {
		nod = nodeCluster[n];
		// check all of it's edges' neighbor triangles
		for (int e = 0; e < nod->getEdgeCount(); ++e) {
			ed = (Edge*)nod->getEdge(e);
			foundInHere = false;
			for (int t = ed->getNeighborTriangulationCount() - 1; t >= 0; --t) {
				tri = ed->getNeighborTriangulation(t);
				// find out what nodes this triangle has
				const TemplateSet<VoronoiNode*> & nodes = tri->getNodeSet();
				// if we're looking for a triangle with the same nodes as this triangle
				if (nodes.size() >= nodeCluster.size() && nodes.containsAll(nodeCluster.getRawListConst(), nodeCluster.size())){
					// that is the triangle. unless we've found another before it, and it's different...
#ifndef NO_TEST
					if (!foundTriangulation) { foundTriangulation = tri; foundInHere = true; }
					else if (foundTriangulation == tri) { foundInHere = true; }
#else
					return tri;
#endif
				}
				// or if this triangle we're looking for is a bigger triangle, that has all of these nodes
				else if (nodeCluster.size() > nodes.size() && nodeCluster.containsAll(nodes.getRawListConst(), nodes.size())) {
					// this one will replace the smaller one
					tri->invalidate(edgePool);
					triangulationPool.markFree(tri);//freeTriangles.add(tri);
				}
			}
		}
	}
	return foundTriangulation;
}

/** call this when you know that this is a valid triangulation... */
DelaunySet::Triangulation * DelaunySet::marshalTriangulation(TemplateSet<VoronoiNode *> & nodeCluster, Circf circumscription) {
	Triangulation * t = getTriangle(nodeCluster);
	if (t == NULL) {
		TemplateVector<Edge*> edges;
		edges.setSize(nodeCluster.size());
		// TODO order the nodeCluster in clockwise order
		TemplateArray<VoronoiNode*> ordered(nodeCluster.size());
		V2f centerMass;
		orderNodesClockwise(nodeCluster.getRawList(), nodeCluster.size(), centerMass, ordered.getRawList());
#ifndef NO_TEST
		if (ordered.containsDuplicates()) { ABORT }
#endif
		// add the edges of the node cluster
		for (int i = 0; i < ordered.size(); ++i) {
			edges[i] = marshalEdge(ordered[i], ordered[(i + 1) % ordered.size()]);
		}
		Triangulation tri = Triangulation(edges.getRawList(), edges.size(), circumscription);
		t = addTriangle(tri);
	}
	if (t != NULL) {
		t->recognizeNeighbors();
	}
#ifndef NO_TEST
	if (t != NULL) {
		// make sure all connections are there.
		for (int i = 0; i < nodeCluster.size(); ++i) {
			if (!t->hasNode(nodeCluster[i])) { ABORT }
		}
	}
#endif
	return t;
}

/** @return true if the given node cluster can be triangulated into a circumscription */
bool DelaunySet::triangulateFor(TemplateSet<VoronoiNode*> & nodeCluster, float floatingPointRounding, Circf & out_circ) {
	switch (nodeCluster.size()) {
	case 0: case 1: { return false; }
	case 2: {
				V2f radius = (nodeCluster[1]->getLocation() - nodeCluster[0]->getLocation()) / 2.0f;
				if (radius.isZero()) return false;
				out_circ.set(nodeCluster[0]->getLocation() + radius, radius.magnitude());
				return true;
	}
	case 3: {
				return V2f::circumcenter(nodeCluster[0]->getLocation(), nodeCluster[1]->getLocation(), nodeCluster[2]->getLocation(), out_circ.center, out_circ.radius);
	}
	default:
		bool circumscriptionWorks;
		// calculate circumscription for each trio of points
		TemplateArray<Circf> circs(nodeCluster.size());
		float minRad = -1, maxRad = -1;
		for (int i = 0; i < circs.size(); ++i) {
			circumscriptionWorks = V2f::circumcenter(nodeCluster[i]->getLocation(),
				nodeCluster[(i + 1) % nodeCluster.size()]->getLocation(),
				nodeCluster[(i + 2) % nodeCluster.size()]->getLocation(),
				circs[i].center, circs[i].radius);
			if (!circumscriptionWorks) return false;
			// keep track of the min/max circles
			if (minRad < 0 || circs[i].radius < minRad) { minRad = circs[i].radius; }
			if (maxRad < 0 || circs[i].radius > maxRad) { maxRad = circs[i].radius; }
		}
		// check if the circles are similar enough (the difference between the min and max circles is small)
		if ((maxRad - minRad) < minRad*floatingPointRounding) {
			// use the average circle location
			out_circ.center = V2f::ZERO();
			for (int i = 0; i < circs.size(); ++i) { out_circ.center += circs[i].center; }
			out_circ.center /= (float)circs.size();
			// with the maximum radius
			out_circ.radius = 0;
			for (int i = 0; i < nodeCluster.size(); ++i) {
				float rad = V2f::distance(out_circ.center, nodeCluster[0]->getLocation());
				if (rad < minRad || rad > maxRad) return false;
				if (rad > out_circ.radius) { out_circ.radius = rad; }
			}
			return true;
		}
		return false;
	}
}

/**
* @param node what to create triangulations around
* @param createdTriangles if not null, adds all new triangles to this set
*/
void DelaunySet::createTriangulationsFor(VoronoiNode* node, TemplateSet<Triangulation*> * createdTriangles) {
	for (int i = 0; i < m_nodePool.size(); ++i) {
		m_nodePool[i].setMask(Obstacle::VORONOI_NODE);
	}
	TemplateSet<VoronoiNode*> nodeCluster;
//	node->setMask(Obstacle::NOTHING); // ensure nodes don't collide with their own triangulations
	nodeCluster.add(node);
	createTriangulationInternal(nodeCluster, 0, createdTriangles, Circf(node->getLocation(), -1));
	for (int i = 0; i < m_nodePool.size(); ++i) {
		m_nodePool[i].setMask(Obstacle::VORONOI_NODE);
	}
}

/** 
 * @param nodeCluster what nodes to try to triangulate with. If there arent enough nodes for triangulation, more will be added. 
 * @param startIndex <-- remove. use CSP instead
 * @param createdTriangles if not null, adds all new triangles to this set
 * @param whereToLookForNodes where to look for more nodes for this triangulation. if the radius is -1, an expansion algorithm will try to find more nodes
 */
bool DelaunySet::createTriangulationInternal(TemplateSet<VoronoiNode*> & nodeCluster, int startIndex, TemplateSet<Triangulation*> * createdTriangles, Circf whereToLookForNodes) {
	const float floatingPointRounding = 1 / 32768.0f;
	Circf circumscription;
	bool triangulationCreated = false;
	bool triangulationClear = false;
	bool atLeastOneMoreComplexCreated = false;
	float dist;
	TemplateSet<Obstacle*> searchResult;
	TemplateVector<Obstacle*> clusterObstacles;
	clusterObstacles.ensureCapacity(nodeCluster.size());
	for (int i = 0; i < nodeCluster.size(); ++i) {
		clusterObstacles.add(nodeCluster[i]);
	}
	for (int i = startIndex; i < m_nodePool.size()/*currentNodes.size()*/; ++i) {
		// grab the next node in the list
		VoronoiNode * n = &m_nodePool[i];// currentNodes[i];
		// don't consider invalid nodes or duplicates in the node cluster
		if (!n->isValidNode() || nodeCluster.indexOf(n) >= 0) continue;
		// don't consider nodes that are too far away
		if (whereToLookForNodes.radius > 0 && (dist = V2f::distance(whereToLookForNodes.center, n->getLocation())) > whereToLookForNodes.radius) continue;
		// add the next non-duplicate node
		n->setMask(Obstacle::NOTHING); // ensure nodes don't collide with their own triangulations
		nodeCluster.add(n);
		bool is3orMorePonts = nodeCluster.size() >= 3;
		// if a triangulation can be made with the current nodes
		if (triangulateFor(nodeCluster, floatingPointRounding, circumscription)) {
			if (nodeCluster.size() >= 4) {
				whereToLookForNodes = circumscription;
				whereToLookForNodes.radius *= 2;
			}
			triangulationClear = false;
			// see if the triangulation is otherwise clear so we can make a triangle here
			searchResult.clear();
			CircleObject circObj(circumscription, Obstacle::VORONOI_NODE);

			Obstacle * blockingNode =
//				getNodeAt(circumscription, 0, nodeCluster.getRawList(), nodeCluster.size());
				getNodeAt(circumscription, searchResult);
//				csp->collidesWithSomething(&circObj);
			// this is only needed if uing getNodeAt(circumscription, searchResult);
			if (blockingNode){
				// ignore the nodes being checked
				searchResult.removeListFast(clusterObstacles);
				if (searchResult.size() == 0)
					blockingNode = NULL;
			}
			triangulationClear = blockingNode == NULL;
			// try to make it a more complex triangulation!
			bool moreComplexPossible = (!is3orMorePonts || triangulationClear) &&
				createTriangulationInternal(nodeCluster, i + 1, createdTriangles, whereToLookForNodes);
			atLeastOneMoreComplexCreated |= moreComplexPossible;
			// if a more complex triangulation was not possible, but this one is clear, use this as a triangulation
			if (triangulationClear && !moreComplexPossible && is3orMorePonts) {
				Triangulation * t = marshalTriangulation(nodeCluster, circumscription);
				if (createdTriangles) createdTriangles->add(t);
				triangulationCreated = true;
			}
		}
		// if the triangulation wasn't clear, remove this recent node addition and try the next one.
		nodeCluster.removeData(n);
		n->setMask(Obstacle::VORONOI_NODE);
	}
	return triangulationCreated;
}

//bool DelaunySet::createTriangulationInternal(TemplateSet<VoronoiNode*> & nodeCluster, int startIndex, TemplateSet<Triangulation*> * createdTriangles, Circf whereToLookForNodes) {
//	const float floatingPointRounding = 1 / 32768.0f;
//	Circf circumscription;
//	bool triangulationCreated = false;
//	bool triangulationClear = false;
//	bool atLeastOneMoreComplexCreated = false;
//	float dist;
//	for (int i = startIndex; i < m_nodePool.size()/*currentNodes.size()*/; ++i) {
//		// grab the next node in the list
//		VoronoiNode * n = &m_nodePool[i];// currentNodes[i];
//		// don't consider invalid nodes or duplicates in the node cluster
//		if (!n->isValidNode() || nodeCluster.indexOf(n) >= 0) continue;
//		// don't consider nodes that are too far away
//		if (whereToLookForNodes.radius > 0 && (dist = V2f::distance(whereToLookForNodes.center, n->getLocation())) > whereToLookForNodes.radius) continue;
//		// add the next non-duplicate node
//		nodeCluster.add(n);
//		bool is3orMorePonts = nodeCluster.size() >= 3;
//		// if a triangulation can be made with the current nodes
//		if (triangulateFor(nodeCluster, floatingPointRounding, circumscription)) {
//			if (nodeCluster.size() >= 4) {
//				whereToLookForNodes = circumscription;
//				whereToLookForNodes.radius *= 2;
//			}
//			triangulationClear = false;
//			// see if the triangulation is otherwise clear so we can make a triangle here
//			VoronoiNode * blockingNode = getNodeAt(circumscription, NULL, nodeCluster.getRawList(), nodeCluster.size());
//			triangulationClear = blockingNode == NULL;
//			// try to make it a more complex triangulation!
//			bool moreComplexPossible = (!is3orMorePonts || triangulationClear) &&
//				createTriangulationInternal(nodeCluster, i + 1, createdTriangles, whereToLookForNodes);
//			atLeastOneMoreComplexCreated |= moreComplexPossible;
//			// if a more complex triangulation was not possible, but this one is clear, use this as a triangulation
//			if (triangulationClear && !moreComplexPossible && is3orMorePonts) {
//				Triangulation * t = marshalTriangulation(nodeCluster, circumscription);
//				if (createdTriangles) createdTriangles->add(t);
//				//printf("t%d(%d) ", nodeCluster.size(), currentTriangles.indexOf(t));
//				//for (int x = 0; x < nodeCluster.size(); ++x) { printf("%d ", nodes.indexOf(nodeCluster[x])); }printf("\n");
//				triangulationCreated = true;
//				//					break;
//			}
//		}
//		// if the triangulation wasn't clear, remove this recent node addition and try the next one.
//		nodeCluster.removeData(n);
//	}
//	return triangulationCreated;
//}

void DelaunySet::calculateAllTriangles() {

	if (/*currentNodes*/m_nodePool.size() > 2) {
		csp->clear();
		for (int i = 0; i < /*currentNodes*/m_nodePool.size(); ++i) {
			csp->add(&m_nodePool[i]);
		}

		for (int i = 0; i < /*currentNodes*/m_nodePool.size(); ++i) {
			createTriangulationsFor(&/*currentNodes*/m_nodePool[i], NULL);
		}
#ifndef NO_TEST
		dupTest();
#endif
	}
#ifndef NO_TEST
	printf("%d triangles\n", triangulationPool.size());
	bigValidationTest();
#endif
}

void DelaunySet::calculateTriangles(VoronoiNode** list, const int listCount) {
	for (int i = 0; i < listCount; ++i) {
		createTriangulationsFor(list[i], NULL);
	}
}

/** @param count how many random points to make in this diagram's boundary */
void DelaunySet::makeRandom(int count) {
	V2f center = boundary->getShape()->getCenter();
	float rad = boundary->getShape()->getRadius();
	V2f r(rad, rad);
	V2f min = center - r;
	V2f max = center + r;
	V2f d = max - min;
	for (int i = 0; i < count; ++i) {
		float x = Random::PRNGf(), y = Random::PRNGf();
		V2f randomPoint = V2f(x*d.x + min.x, y*d.y + min.y);
		if (boundary->getShape()->contains(randomPoint)) {
			/*currentNodes*/m_nodePool.add(VoronoiNode(randomPoint));
		}
		else {
			--i;
		}
	}
}

void DelaunySet::gatherTrianglesContaining(V2f const & point, TemplateSet<Triangulation*> & out_triangleSet) {
	// TODO use triangl
	for (int i = 0; i < triangulationPool.size(); ++i) {
		Triangulation * t = &triangulationPool[i];
		if (t->circumscriptionContains(point)) {
			out_triangleSet.add(t);
		}
	}
}

int DelaunySet::indexOf(const Triangulation * t, TemplateVector<const Triangulation*> & triangleList) {
	for (int i = 0; i < triangleList.size(); ++i) {
		if (triangleList[i]->equals(*t)) return i;
	}
	return -1;
}

bool DelaunySet::findDuplicateTriangle(int startIndex, int & a, int & b) {
	Triangulation * ta, *tb;
	for (a = startIndex; a < triangulationPool.size(); ++a) {
		ta = &triangulationPool[a];
		if (ta->isValidTriangle()) {
			for (b = a + 1; b < triangulationPool.size(); ++b) {
				tb = &triangulationPool[b];
				if (tb->isValidTriangle() && ta->equals(*tb)) {
					return true;
				}
			}
		}
	}
	return false;
}

void DelaunySet::dupTest() {
#ifndef NO_TEST
	int a, b, index = -1;
	bool found;
	do {
		found = findDuplicateTriangle(index + 1, a, b);
		if (found) {
			printf("duplicate triangles: %d and %d  ", a, b);
			printf("\n");
			index = a;
		}
	} while (found);
#endif
}

void DelaunySet::bigValidationTest() {
#ifndef NO_TEST
	// make sure all deleted triangles and edges are marked as such
	for (int i = 0; i < triangulationPool.size(); ++i) {
		Triangulation * t = &triangulationPool[i];
		bool valid = t->isValidTriangle();
		bool inFreeSet = triangulationPool.isMarkedFree(t);
		if (!valid && !inFreeSet) {
			printf("fail. invalid triangle not in the free triangle set."); { ABORT }
		}
		else if (valid && inFreeSet) {
			printf("fail. valid triangle in the free triangle set."); { ABORT }
		}
		if (valid) {
			// make sure crrent valid triangles and edges don't point at invalid ones
			for (int j = 0; j < t->getEdgeCount(); ++j) {
				Edge * e = t->getEdge(j);
				if (!e->isValid()) {
					printf("fail. valid triangle pointing at invalid edge."); { ABORT }
				}
			}
		}
	}
	for (int i = 0; i < edgePool.size(); ++i) {
		Edge * e = &edgePool[i];
		bool valid = e->isValid();
		bool inFreeSet = edgePool.isMarkedFree(e);
		if (!valid && !inFreeSet) {
			printf("fail. invalid edge not in the free edge set."); { ABORT }
		}
		else if (valid && inFreeSet) {
			printf("fail. valid edge in the free edge set."); { ABORT }
		}
		if (valid) {
			// make sure crrent valid triangles and edges don't point at invalid ones
			for (int j = 0; j < e->getNeighborTriangulationCount(); ++j) {
				Triangulation * t = e->getNeighborTriangulation(j);
				if (!t->isValidTriangle()) {
					printf("fail. valid edge pointing at invalid triangle."); { ABORT }
				}
			}
		}
	}
	// check if nodes are pointing at any invalid nodes
	for (int i = 0; i < m_nodePool.size(); ++i) {
		VoronoiNode * vn = &m_nodePool[i];
		for (int n = 0; n < vn->getEdgeCount(); ++n){
			Edge * e = dynamic_cast<Edge*>(vn->getEdge(n));
			if (!e->isValid()) {
				printf("fail. node pointing at invalid edge."); { ABORT }
			}
		}
	}
	// check if every VALID triangle's nodes agree that they are in the same triangle
	for (int i = 0; i < triangulationPool.size(); ++i) {
		Triangulation * t = &triangulationPool[i];
		if (!t->isValidTriangle()) continue;
		//			if (t->edges.size() != 3) { printf("fail. triangles must have 3 edges."); ABORT }
		for (int j = 0; j < t->getEdgeCount(); ++j) {
			Edge * e = dynamic_cast<Edge*>(t->getEdge(j));
			if (!e->isNeighbor(t)) {
				printf("fail. edge is not back-connecting to the triangle (%08x).\n", t);
				printf("edge %08x connects to node %08x and %08x\n", e, e->getNode(0), e->getNode(1));
				for (int i = 0; i < e->getNeighborTriangulationCount(); ++i) {
					printf("neighbor tri %08x\n", e->getNeighborTriangulation(i));
				}
				ABORT
			}
			VoronoiNode * n0 = dynamic_cast<VoronoiNode*>(e->getNode(0)), *n1 = dynamic_cast<VoronoiNode*>(e->getNode(1));
			if (n0->getEdgeTo(n1) != e) {
				printf("fail. node in edge is weirdly not connected correctly."); ABORT
			}
			if (n1->getEdgeTo(n0) != e) {
				printf("fail. node in edge is weirdly not connected correctly."); ABORT
			}
		}
	}
#endif
}

AbstractGraphNode * DelaunySet::createNode() {
	//VoronoiNode * newNode;
	//if (freeNodes.size() > 0) {
	//	newNode = freeNodes.pop();
	//}
	//else {
	//	currentNodes.add();
	//	newNode = &currentNodes.getLast();
	//}
	//return newNode;
	return m_nodePool.add();
}

bool DelaunySet::addNode(V2f point) {
	TemplateSet<VoronoiNode*> changedNodes;
	bool added = addNode(point, changedNodes);
	// recalculate triangles that were modified
	calculateTriangles(changedNodes.getRawList(), changedNodes.size());
#ifndef NO_TEST
	bigValidationTest();
	dupTest();
#endif
	return added;
}

bool DelaunySet::addNode(V2f point, TemplateSet<VoronoiNode*> & changedNodes) {
#ifndef NO_TEST
	bigValidationTest();
#endif
	// if this node is not in the bondaries of this diagram, ignore it
	if (boundary && !boundary->getShape()->contains(point)) return false;

	// if this node already exists, ignore it
	//for (int i = 0; i < /*currentNodes*/nodePool.size(); ++i) {
	//	if (/*currentNodes*/nodePool[i].getLocation() == point){
	//		return false;
	//	}
	//}
	TemplateSet<Obstacle*> result;
	csp->gatherAtPosition(point, result, Obstacle::VORONOI_NODE);
	if (result.size() > 0) return false;

	VoronoiNode * newNode = (VoronoiNode*)createNode();
	newNode->set(point);

	csp->add(newNode);
#ifndef NO_TEST
	bigValidationTest();
#endif

	// calculate the triangles that are going to be destroyed
	TemplateSet<Triangulation*> brokenTrianglesSet;
	gatherTrianglesContaining(point, brokenTrianglesSet);
	// calculate the nodes that will need to re-calculate triangles
	changedNodes.add(newNode);
	for (int i = 0; i < brokenTrianglesSet.size(); ++i) {
		brokenTrianglesSet[i]->gatherNodesInto(changedNodes);
	}
#ifndef NO_TEST
	bigValidationTest();
#endif
	// remove broken triangles
	for (int i = 0; i < brokenTrianglesSet.size(); ++i) {
		removeTriangle(brokenTrianglesSet[i]);
	}
#ifndef NO_TEST
	bigValidationTest();
#endif
	return true;
}

/**
* @param location where to look for nodes
* @param out_allNodesHere a list to put all of the discovered nodes in. If NULL, the function will stop when the first node is found
* @param ignoreList what nodes to ignore
* @param ignoreListCount how many nodes to ignore
* @return a node found in the given location
*/
DelaunySet::VoronoiNode * DelaunySet::getNodeAt(Circf const & location, TemplateSet<VoronoiNode*> * out_allNodesHere, VoronoiNode** ignoreList, const int ignoreListCount) {
	VoronoiNode * n, *foundOne = NULL;
	for (int i = 0; i < /*currentNodes*/m_nodePool.size(); ++i) {
		n = &/*currentNodes*/m_nodePool[i];
		if (ignoreList && TemplateArray<VoronoiNode*>::indexOf(n, ignoreList, 0, ignoreListCount) >= 0) continue;
		if (n->isValidNode() && location.contains(n->getLocation())) {
			foundOne = n;
			if (!out_allNodesHere) break;
			out_allNodesHere->add(n);
		}
	}
	return foundOne;
}

/**
* @param location where to look for nodes
* @param out_allNodesHere a list to put all of the discovered nodes in. If NULL, the function will stop when the first node is found
* @return a node found in the given location
*/
Obstacle * DelaunySet::getNodeAt(Circf const & location, TemplateSet<Obstacle*> & out_allNodesHere) {
	csp->gatherAtCircle(location, out_allNodesHere, Obstacle::VORONOI_NODE);
	return (out_allNodesHere.size() == 0) ? NULL : out_allNodesHere.getLast();
}

//DelaunySet::VoronoiNode * DelaunySet::getNodePolyhedronContains(V2f const & point) {
//	for (int i = 0; i < /*currentNodes*/nodePool.size(); ++i) {
//		VoronoiNode * n = &/*currentNodes*/nodePool[i];
//		if (n->isValidNode() && n->polyhedronContains(point)) {
//			return n;
//		}
//	}
//	return NULL;
//}
//
//void DelaunySet::gatherVoronoi(TemplateVector<VoronoiNode*> & nodes, bool includeBorderPolygons) {
//	for (int i = 0; i < /*currentNodes*/nodePool.size(); ++i) {
//		VoronoiNode * vn = &/*currentNodes*/nodePool[i];
//		if (vn->isValidNode()) {
//			vn->getPolygon2f(boundary);
//			if (includeBorderPolygons || !vn->isBorderPolygon()) {
//				nodes.add(vn);
//			}
//		}
//	}
//}

void DelaunySet::draw(GLUTRenderingContext * g) const {
	if (boundary) boundary->draw(g, false);
	csp->draw(g);
#ifndef NO_TEST
	for (int i = 0; i < triangulationPool.size(); ++i) {
		if (!triangulationPool[i].isValidTriangle()) continue;
		g->setColor(0xffddff);
		triangulationPool[i].circum.glDraw(false);
		g->setColor(0xff00ff);
		g->printf(triangulationPool[i].centerMass, "%d", i);
		triangulationPool[i].drawEdges(g);
	}
#endif

	for (int i = 0; i < edgePool.size(); ++i) {
		const Edge * e = &edgePool[i];
		if (e->isValid()) {
			if (e->getCost() > 0) {
				g->setColor(0x88ff88);
				e->draw(g);
#ifndef NO_TEST
				// draw which triangles are neighboring the edges. useful for debugging, if you think triangles might be overlapping
				g->setColor(0x00aa00);
				VoronoiNode * a = (VoronoiNode*)e->getNode(0);
				VoronoiNode * b = (VoronoiNode*)e->getNode(1);
				V2f c = V2f::between(a->getLocation(), b->getLocation());
				g->printf(c, "%d %.1f", e->getNeighborTriangulationCount(), e->getCost());
				for (int n = 0; n < e->getNeighborTriangulationCount(); ++n) {
					Triangulation * t = e->getNeighborTriangulation(n);
					V2f delta = t->centerMass - c;
					g->drawLine(c, c + delta*0.9f);
				}
#endif
			}
		}
	}

	if (selectedNode) {
		const Triangulation * t;
		TemplateSet<Triangulation*> trianglefaces;
		selectedNode->gatherTriangles(trianglefaces);
		TemplateSet<VoronoiNode*> triangleNodes;
		for (int i = 0; i < trianglefaces.size(); ++i) {
			t = trianglefaces[i];
			g->setColor(0xdddddd);
			if (t->getEdgeCount() >= 3) {
				triangleNodes.clear();
				t->gatherNodesInto(triangleNodes);
				// TODO draw a triangle fan from the centerMass, along the edges in clock-wise order
				glBegin(GL_TRIANGLE_FAN);
				for (int p = 0; p < triangleNodes.size(); ++p) {
					triangleNodes[p]->getLocation().glVertex();
				}
				glEnd();
			}
		}
		g->setColor(0xff00ff);
		for (int i = 0; i < selectedNode->getEdgeCount(); ++i) {
			const Edge * e = (Edge*)selectedNode->getEdge(i);
#ifndef NO_TEST
			if (!e->isValid()) {
				ABORT
			}
#endif
			e->draw(g);
		}
	}
	const Triangulation * selectedTriangle = NULL;
	for (int i = 0; i < triangulationPool.size(); ++i) {
		const Triangulation* t = &triangulationPool[i];
		if (t->isValidTriangle() && t->polygonContains(specialCursor)) {
			selectedTriangle = t;
		}
	}

	g->setColor(0xaa0000);
	for (int i = 0; i < /*currentNodes*/m_nodePool.size(); ++i) {
		const VoronoiNode * n = &/*currentNodes*/m_nodePool[i];
		if (n->isValidNode()) {
			glDrawCircle(n->getLocation(), 0.2f, false);
			g->printf(n->getLocation(), "%d", i);
		}
	}

	TemplateVector<V2f> edgeBorder;
	for (int i = 0; i < triangulationPool.size(); ++i) {
		if (!triangulationPool[i].isValidTriangle()) continue;

		glDrawCircle(triangulationPool[i].circum.center, .05f, false);

		const Triangulation * t = &triangulationPool[i];

		for (int n = 0; n < t->getEdgeCount(); ++n) {
			Edge * e = t->getEdge(n);

			edgeBorder.clear();
			V2f edgeCenter;
			for (int edgeNeighbor = 0; edgeNeighbor < e->getNeighborTriangulationCount(); ++edgeNeighbor) {
				const Triangulation * ent = e->getNeighborTriangulation(edgeNeighbor);
				edgeBorder.add(ent->circum.center);
				edgeCenter += ent->circum.center;
			}
			edgeCenter /= (float)edgeBorder.size();
			VoronoiNode * a = (VoronoiNode*)e->getNode(0);
			VoronoiNode * b = (VoronoiNode*)e->getNode(1);
			V2f delta = a->getLocation() - b->getLocation();
			// TODO if this is 3D, order the points to be clockwise around the edge vector
			// draw the separation
			if (edgeBorder.size() > 1) {
				for (int i = 0; i < edgeBorder.size(); ++i) {
					edgeCenter.glDrawTo(edgeBorder[(i + 1) % edgeBorder.size()]);
				}
			}
			
			// if this edge has area that is not triangulated
			if (edgeBorder.size() == 1) {
				g->setColor(0x0000aa);
				g->drawCircle(edgeCenter, .1f, false);
				g->setColor(0xaa0000);
				V2f rayDir = e->getFace().normal.perp(), end;//-(delta.perp().normal()), end;
				V2f otherPoint;
				a = (VoronoiNode*)e->getNode(0);
				b = (VoronoiNode*)e->getNode(1);
				V2f::closestPointOnLine(a->getLocation(), b->getLocation(), t->centerMass, otherPoint);
				V2f correctDirection = otherPoint - t->centerMass;
				float alignment = V2f::dot(rayDir, correctDirection);
				if (alignment < 0) { rayDir *= -1.0f; }
				if (!boundary) {
					end = edgeCenter + rayDir * 100;
					edgeCenter.glDrawTo(end);
				}
				else if (boundary->getShape()->contains(edgeCenter)) {
					V2f end;
					RaycastHit rh;
					if (boundary->getShape()->raycast(Ray(edgeCenter, rayDir), rh)){// dist, end, hitNorm)) {
						end = rh.point;//edgeCenter;
					}
					edgeCenter.glDrawTo(end);
				}
			}
			
		}
	}

	if (selectedNode) {

		TemplateVector<V2f> points;
		V2f center;
		selectedNode->calculatePolygon(points, center);
		if (points.size() > 0) {
			g->setColor(0xff00ff);
			glBegin(GL_TRIANGLE_FAN);
			center.glVertex();
			V2f::glVertexList(points.getRawList(), points.size());
			points[0].glVertex();
			glEnd();

			for (int i = 0; i < points.size(); ++i) {
				V2f delta = center - points[i];
				delta.normalize();
			}
		}
		Polygon2f * calcedPoly = &selectedNode->getPolygon2f(boundary);
		g->setColor(selectedNode->isBorderPolygon() ? 0x0000ff : 0xffff00);
		calcedPoly->glDraw(true);
	}
}