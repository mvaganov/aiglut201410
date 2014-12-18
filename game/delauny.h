#pragma once

#include <GL/freeglut.h>
#include "codegiraffe/v2.h"
#include "codegiraffe/templatevector.h"
#include "codegiraffe/circle.h"
#include "codegiraffe/templatevectorlist.h"
#include "codegiraffe/glutrenderingcontext.h"
#include "codegiraffe/templateset.h"
#include "obstacles.h"
#include "graph.h"

// TODO voronoi around obstacles to make a nav mesh
// TODO profile and optimize

#define ABORT { int NUMBER = 0; NUMBER = 1 / NUMBER; }

#define NO_TEST

/** a structure to manage delauny triangulation */
class DelaunySet : public AbstractGraph {
public:
	// forward declarations
	class VoronoiNode;
	class Edge;
	class Triangulation;

	/** the faces of a VoronoiNode, which lie on the edges between nodes TODO make this a polygon... */
	class VoronoiFace {
	public:
		/** the points that determine the ends of the separating edge. in a 3D shape these will be the coordinates of the polygon face between voronoi nodes */
		TemplateVector<V2f> points;
		/** where the separating face is, and it's normal */
		V2f center, normal;
		float radius;

		VoronoiFace(){}

		void refresh(Edge * e) {
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
		void glDraw() const {
			for (int i = 0; i < points.size(); ++i) {
				center.glDrawTo(points[i]);
			}
			center.glDrawTo(center + normal);
		}
		/** positive value on one side, negative on the other. not sure which is which. */
		float sideValue(V2f const & p) { return V2f::dot(normal, p - center); }

		void gatherOppositePoints(V2f & out_a, V2f & out_b) {
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
	};

	/** the graph edge of the voronoi node, used to calculate the delauny triangulation */
	class Edge : public GraphEdge {
		/** the triangulations that neighbor this edge */
		TemplateSet<Triangulation*> neighborTri;
		/** the face that was calculated for this edge */
		VoronoiFace voronoiFace;
	public:
		VoronoiFace & getFace() { return voronoiFace; }

		int getNeighborTriangulationCount() const { return neighborTri.size(); }

		Triangulation* getNeighborTriangulation(const int i) const { return neighborTri[i]; }

		void refreshVoronoiFace() { voronoiFace.refresh(this); }

		void addNeighborTriangulation(Triangulation * nt) {
			bool added = neighborTri.add(nt);
			if (added) refreshVoronoiFace();
		}
		bool removeNeighborTriangulation(Triangulation * nt) {
			bool removed = neighborTri.removeData(nt);
			if (removed) refreshVoronoiFace();
			return removed;
		}

		void invalidateEdge() {
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
		bool isValid() const { return getNode(0) != NULL && getNode(1) != NULL; }
		Edge() { nodes[0] = nodes[1] = NULL; invalidateEdge(); }
		Edge(VoronoiNode *a, VoronoiNode *b) : GraphEdge(a, b){}// { nodes[0] = a; nodes[1] = b; }
		void glDraw() const {
			GraphNode * a = (GraphNode*)getNode(0), *b = (GraphNode*)getNode(1);
			a->getLocation().glDrawTo(b->getLocation());
		}
		bool isNeighbor(Triangulation * const t) const { return neighborTri.indexOf(t) >= 0; }
	};

	class Triangulation {
	private:
		/** not a set because the order matters. edges are stored in clockwise fashion. */
		TemplateArray<Edge*> edges;
		/** cached calculation of nodes in this Triangulation */
		TemplateSet<VoronoiNode*> nodeSet;
	public:
		/** circumscription of all edges */
		CircF circum;
		/** the average point of all of the node locations */
		V2f centerMass;
		/** which node to start drawing from */
		VoronoiNode * startingNode;
		void glDrawEdges() const { for (int i = 0; i < edges.size(); ++i) { edges[i]->glDraw(); } }

		Triangulation & operator=(Triangulation const & t) {
			circum = t.circum;
			edges = t.edges;
			startingNode = t.startingNode;
			centerMass = t.centerMass;
			nodeSet = t.nodeSet;
			return *this;
		}

		void set(Edge* const * list, const int listCount, CircF const & circumscription) {
			startingNode = NULL;
			circum = circumscription;
			edges.allocateToSize(listCount);
			for (int i = 0; i < listCount; ++i) { edges[i] = list[i]; }
			calculate();
		}

		Triangulation(Edge* const * list, const int listCount, CircF const & circumscription) : startingNode(NULL) {
			set(list, listCount, circumscription);
		}

		/** create an invalid triangulation by default */
		Triangulation() : startingNode(NULL), circum(0, 0, -1) { }

		int getEdgeCount() const { return edges.size(); }
		Edge * getEdge(const int i) const { return edges[i]; }

		/**
		 * calculate center of mass, and order edges in clockwise fashion.
		 * @return true if this is a valid Triangulation
		 */
		bool calculate() {
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
		bool isValidTriangle() const { return circum.radius >= 0; }

		/** make this an invalid triangulation. put any edges freed (nolonger referenced by any triangles) into the given set */
		void invalidate(TemplateSet<Edge*> & freeEdges) {
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
					freeEdges.add(e);
				}
				edges[i] = NULL; // don't de-allocate the edges, this triangulation will probably be used again, and might need the same number of edges
			}
		}

		/** put all of the nodes from referenced edges into the given set */
		void gatherNodesInto(TemplateSet<VoronoiNode*> & nodes) const {
			for (int i = 0; i < edges.size(); ++i) {
				nodes.add((VoronoiNode*)edges[i]->getNode(0));
				nodes.add((VoronoiNode*)edges[i]->getNode(1));
			}
		}

		const TemplateSet<VoronoiNode*> & getNodeSet() const { return nodeSet; }

		/** only works for 2D triangulation */
		void gatherPointsClockwise(TemplateArray<V2f> & points) const {
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
		bool polygonContains(V2f const & p) const {
			if (!this->isValidTriangle()) return false;
			if (p.distance(circum.center) > circum.radius) return false;
			TemplateArray<V2f> points;
			gatherPointsClockwise(points);
			return p.isInsidePolyCW(points.getRawList(), points.size());
		}

		/** @return true if this triangulation is using the given edge */
		bool hasEdge(const Edge * const edge) const {
			for (int i = 0; i < edges.size(); ++i) {
				if (edges[i] == edge) return true;
			} return false;
		}

		/** @return true if this triangle referenced the given triangle as a neighbor */
		bool isNeighborsWith(Triangulation * t) const {
			for (int i = 0; i < edges.size(); ++i) {
				if (edges[i]->isNeighbor(t)) return true;
			}
			return false;
		}

		/** get the edge that connects the two given nodes. useful when the list of nodes is known but not the edges */
		int findEdgeIndex(VoronoiNode * n0, VoronoiNode * n1) const {
			for (int i = 0; i < edges.size(); ++i) {
				if (edges[i]->hasBoth(n0, n1)) return i;
			}
			return -1;
		}

		/** called after a Triangulation is created, to ensure that it's edges are aware of it */
		void recognizeNeighbors() {
			for (int i = 0; i < edges.size(); ++i) {
				edges[i]->addNeighborTriangulation(this);
			}
		}

		/** check this triangulation's circumscription */
		bool circumscriptionContains(V2f const & p) const { return circum.contains(p); }

		/** @return if this Triangulation is related to the given node */
		bool hasNode(const VoronoiNode * node) const {
			for (int i = 0; i < edges.size(); ++i) { if (edges[i]->has(node)) return true; }
			return false;
		}
		/** @return if this Triangulation has every node in the given series */
		bool hasAllNodes(VoronoiNode ** list, const int listSize) const {
			for (int i = 0; i < listSize; ++i) {
				if (!hasNode(list[i])) return false;
			}
			return true;
		}
		/** @return if this node has exactly the nodes in this given series */
		bool hasExactlyTheseNodes(VoronoiNode ** list, const int listSize) const {
			const TemplateSet<VoronoiNode*> & nodes = getNodeSet();
			//gatherNodesInto(nodes);
			if (nodes.size() != listSize) return false;
			for (int i = 0; i < listSize; ++i) {
				if (nodes.indexOf(list[i]) == -1) return false;
			}
			return true;
		}

		/** a triangulation is equal if it is brokering between the exact same nodes */
		bool equals(const Triangulation & t) const {
			if (&t == this) return true;
			return getNodeSet() == t.getNodeSet();
		}
		bool operator==(Triangulation const & t) const { return equals(t); }
		bool operator!=(Triangulation const & t) const { return !equals(t); }
	};

	class VoronoiNode : public GraphNode {
		Polygon2f polygon;
		/** identifies if this is on the edge of the voronoi diagram (near the boundaries) */
		Obstacle * atBoundaryOf;
		/** */
		bool valid;
		/** whenever data changes, the node should flag itself as needed recalculation */
		bool needsModelRecalculated;

	public:
		VoronoiNode() : atBoundaryOf(NULL), valid(false), needsModelRecalculated(true){}
		VoronoiNode(V2f const & p) { set(p); }

		const V2f & getLocation() const { return location; }

		/** @return whether the polygon needs to be recalculated */
		bool isDirty() { return needsModelRecalculated; }
		/** identify that the polygon needs to be recalculated */
		void setDirty() { needsModelRecalculated = true; }

		bool isBorderPolygon() const { return atBoundaryOf != NULL; }

		int getEdgeCount() const { return edges.size(); }
		//		Edge * getEdge(const int i) const { return edges[i]; }
		bool removeEdge(Edge * e) { setDirty(); return edges.removeData(e); }
		//		bool addEdge(Edge * e) { setDirty(); return edges.add(e); }

		void set(V2f const & p) { location = p; valid = true; setDirty(); atBoundaryOf = NULL; }

		bool isValidNode() const { return valid; }

		void invalidate(DelaunySet & D, TemplateSet<VoronoiNode*> & changedNodes) {
			while (edges.size() > 0) {
				Edge * e = dynamic_cast<Edge*>(edges.getLast());
				for (int i = e->getNeighborTriangulationCount() - 1; i >= 0; --i) {
					Triangulation * t = e->getNeighborTriangulation(i);
					t->gatherNodesInto(changedNodes);
					t->invalidate(D.freeEdges);
					D.freeTriangles.add(t);
				}
			}
			valid = false;
			setDirty();
		}

		/** calculate the Voronoi polygon*/
		void calculatePolygon(TemplateVector<V2f> & out_points, V2f & out_center) const {
			TemplateSet<Triangulation*> triangles;
			TemplateSet<VoronoiNode*> nodes;
			gatherTriangles(triangles);
			for (int i = 0; i < triangles.size(); ++i) {
				out_points.add(triangles[i]->circum.center);
			}
			Polygon2f::calculatePolygonCW(out_points.getRawList(), out_points.size(), out_center);
		}

		/** calculate the Voronoi polygon*/
		Polygon2f & getPolygon2f(Obstacle * boundary) {
			if (needsModelRecalculated) {
				atBoundaryOf = NULL;
				//printf("-------------------\n");
				//TemplateVector<LineF> lines;
				polygon.points.clear();
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
							V2f ray = e->getFace().normal.perp(), end, otherPoint;//-(delta.perp().normal()), end;
							GraphNode * a = (GraphNode*)e->getNode(0);
							GraphNode * b = (GraphNode*)e->getNode(1);
							V2f::closestPointOnLine(a->getLocation(), b->getLocation(), t->centerMass, otherPoint);
							V2f correctDirection = otherPoint - t->centerMass;
							float alignment = V2f::dot(ray, correctDirection);
							if (alignment < 0) { ray *= -1.0f; }
							// raycast the dividing plane to the border of the boundary
							if (!boundary) {
								end = t->circum.center + ray * 100;
							}
							else if (boundary->contains(t->circum.center)) {
								V2f hitNorm;
								float dist;
								if (boundary->raycast(t->circum.center, ray, dist, end, hitNorm)) {
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
							bool thisInBoundary = boundary->contains(thisPoint);
							bool thatInBoundary = boundary->contains(thatPoint);
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
				polygon.set(this->location, 0, points.getRawListConst(), points.size(), pairs.getRawListConst(), pairs.size(), true);
				//for (int i = 0; i < polygon.points.size(); ++i) { printf("(%.1f %.1f)", polygon.points[i].x, polygon.points[i].y); } printf("\n");
				needsModelRecalculated = false;
			}
			return polygon;
		}

		bool polyhedronContains(V2f const & a_position) const {
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
		Edge * getEdgeTo(VoronoiNode * node) const {
			Edge * e;
			for (int i = 0; i < edges.size(); ++i) {
				e = (Edge*)edges[i];
				if (e->has(node))
				{
#ifndef NO_TEST
					if (!e->equals(node, this)) {
						printf("weird... this edge doesn't include *this"); ABORT
					}
#endif
					return e;
				}
			}
			return NULL;
		}

		/** gather a list of every triangulation attached to this node */
		void gatherTriangles(TemplateSet<Triangulation*> & triangles) const {
			for (int i = 0; i < edges.size(); ++i) {
				Edge * e = (Edge*)edges[i];
				for (int n = 0; n < e->getNeighborTriangulationCount(); ++n) {
					triangles.add(e->getNeighborTriangulation(n));
				}
			}
		}

		void gatherEdges(TemplateSet<Edge*> & edges) const {
			for (int i = 0; i < edges.size(); ++i) {
				edges.add(edges[i]);
			}
		}
	};

	/**
	* @param nodes what nodes to calculate the clockwise order of
	* @param nodeCount how many nodes there are
	* @param out_centerMass the average location of the nodes
	* @param out_inOrder where to put the list re-arranged in order
	*/
	static void orderNodesClockwise(VoronoiNode*const* nodes, const int nodeCount, V2f  & out_centerMass, VoronoiNode** out_inOrder) {
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

	int getNodeCount() const { return currentNodes.size(); }
	AbstractGraphNode * getNode(const int index) { return &currentNodes[index]; }

	void connectNodes(AbstractGraphNode * from, AbstractGraphNode * to) {
		getEdge(from, to, true);
	}
	AbstractGraphEdge * getEdge(AbstractGraphNode * a, AbstractGraphNode * b, bool createIfNotThere) {
		VoronoiNode* va = (VoronoiNode*)a, *vb = (VoronoiNode*)b;
		Edge * e = findEdge(va, vb);
		if (e == NULL && createIfNotThere) {
			e = addEdge(Edge(va, vb));
		}
		a->addEdge(e);
		b->addEdge(e);
		return e;
	}

	TemplateVectorList<VoronoiNode> currentNodes;
	TemplateVectorList<Edge> currentEdges;
	TemplateVectorList<Triangulation> currentTriangles;

	TemplateSet<VoronoiNode*> freeNodes;
	TemplateSet<Triangulation*> freeTriangles;
	TemplateSet<Edge*> freeEdges;
	VoronoiNode * selectedNode;
	Obstacle * boundary;

	DelaunySet(Obstacle * boundary) :selectedNode(0), boundary(boundary){}

	~DelaunySet() { }

	Triangulation * addTriangle(Triangulation const & t) {
		Triangulation * tri;
		if (freeTriangles.size() > 0) {
			tri = freeTriangles.pop();
		}
		else {
			currentTriangles.add();
			tri = &currentTriangles.getLast();
		}
		*tri = t;
		return tri;
	}

	void removeTriangle(Triangulation * toRemove) {
		toRemove->invalidate(freeEdges);
		freeTriangles.add(toRemove);
#ifndef NO_TEST
		// check to see if anything references the bad triangle still.
		for (int i = 0; i < currentTriangles.size(); ++i) {
			Triangulation * t = &currentTriangles[i];
			if (t->isValidTriangle() && t != toRemove) {
				if (t->isNeighborsWith(toRemove)) {
					printf("bad triangle neighbor.\n"); ABORT
				}
			}
		}
#endif
	}

	void destroyNode(AbstractGraphNode * n) {
		removeNode((VoronoiNode*)n);
	}

	void removeNode(VoronoiNode * toRemove) {
		TemplateSet<VoronoiNode*> changedNodes;
		removeNode(toRemove, changedNodes);
		calculateTriangles(changedNodes.getRawList(), changedNodes.size());
	}

	void removeNode(VoronoiNode * toRemove, TemplateSet<VoronoiNode*> & changedNodes) {
		toRemove->invalidate(*this, changedNodes);
		freeNodes.add(toRemove);
		changedNodes.removeData(toRemove);
	}

	void moveNode(VoronoiNode * toMove, V2f newPosition) {
		TemplateSet<VoronoiNode*> changedNodes;
		removeNode(toMove, changedNodes);
		addNode(newPosition, changedNodes);
		calculateTriangles(changedNodes.getRawList(), changedNodes.size());
	}

	Edge * findEdge(VoronoiNode * a, VoronoiNode * b) {
		for (int i = 0; i < currentEdges.size(); ++i) {
			if (currentEdges[i].hasBoth(a, b)) return &currentEdges[i];
		}
		return NULL;
	}

	Edge* addEdge(Edge const & e) {
		if (freeEdges.size() > 0) {
			Edge * ed = freeEdges.pop();
			*ed = e;
			return ed;
		}
		currentEdges.add(e); return &currentEdges.getLast();
	}

	Edge* marshalEdge(VoronoiNode * a, VoronoiNode * b) {
		Edge * e = findEdge(a, b);
		if (e == NULL) {
			e = addEdge(Edge(a, b));
		}
		a->addEdge(e);
		b->addEdge(e);
		return e;
	}

	Triangulation * getTriangle(TemplateSet<VoronoiNode*> & nodeCluster) {
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
				for (int t = ed->getNeighborTriangulationCount()-1; t >= 0; --t) {
					tri = ed->getNeighborTriangulation(t);
					// find out what nodes this triangle has
					const TemplateSet<VoronoiNode*> & nodes = tri->getNodeSet();
					// if we're looking for a triangle with the same nodes as this triangle
					if (nodes.size() >= nodeCluster.size() && nodes.containsAll(nodeCluster.getRawListConst(), nodeCluster.size()) ){
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
						tri->invalidate(freeEdges);
						freeTriangles.add(tri);
					}
				}
			}
		}
		return foundTriangulation;
	}

	/** call this when you know that this is a valid triangulation... */
	Triangulation * marshalTriangulation(TemplateSet<VoronoiNode *> & nodeCluster, CircF circumscription) {
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
	bool triangulateFor(TemplateSet<VoronoiNode*> & nodeCluster, float floatingPointRounding, CircF & out_circ) {
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
			TemplateArray<CircF> circs(nodeCluster.size());
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
	void createTriangulationsFor(VoronoiNode* node, TemplateSet<Triangulation*> * createdTriangles) {
		TemplateSet<VoronoiNode*> nodeCluster;
		nodeCluster.add(node);
		createTriangulationInternal(nodeCluster, 0, createdTriangles, CircF(node->getLocation(), -1));
	}

	bool createTriangulationInternal(TemplateSet<VoronoiNode*> & nodeCluster, int startIndex, TemplateSet<Triangulation*> * createdTriangles, CircF whereToLookForNodes) {
		const float floatingPointRounding = 1 / 32768.0f;
		CircF circumscription;
		bool triangulationCreated = false;
		bool triangulationClear = false;
		bool atLeastOneMoreComplexCreated = false;
		float dist;
		for (int i = startIndex; i < currentNodes.size(); ++i) {
			// grab the next node in the list
			VoronoiNode * n = &currentNodes[i];
			// don't consider invalid nodes or duplicates in the node cluster
			if (!n->isValidNode() || nodeCluster.indexOf(n) >= 0) continue;
			// don't consider nodes that are too far away
			if (whereToLookForNodes.radius > 0 && (dist = V2f::distance(whereToLookForNodes.center, n->getLocation())) > whereToLookForNodes.radius) continue;
			// add the next non-duplicate node
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
				VoronoiNode * blockingNode = getNodeAt(circumscription, NULL, nodeCluster.getRawList(), nodeCluster.size());
				triangulationClear = blockingNode == NULL;
				// try to make it a more complex triangulation!
				bool moreComplexPossible = (!is3orMorePonts || triangulationClear) &&
					createTriangulationInternal(nodeCluster, i + 1, createdTriangles, whereToLookForNodes);
				atLeastOneMoreComplexCreated |= moreComplexPossible;
				// if a more complex triangulation was not possible, but this one is clear, use this as a triangulation
				if (triangulationClear && !moreComplexPossible && is3orMorePonts) {
					Triangulation * t = marshalTriangulation(nodeCluster, circumscription);
					if (createdTriangles) createdTriangles->add(t);
//printf("t%d(%d) ", nodeCluster.size(), currentTriangles.indexOf(t));
//for (int x = 0; x < nodeCluster.size(); ++x) { printf("%d ", nodes.indexOf(nodeCluster[x])); }printf("\n");
					triangulationCreated = true;
//					break;
				}
			}
			// if the triangulation wasn't clear, remove this recent node addition and try the next one.
			nodeCluster.removeData(n);
		}
		return triangulationCreated;
	}

	void calculateAllTriangles() {
		if (currentNodes.size() > 2) {
			for (int i = 0; i < currentNodes.size(); ++i) {
				createTriangulationsFor(&currentNodes[i], NULL);
			}
#ifndef NO_TEST
			dupTest();
#endif
		}
#ifndef NO_TEST
		printf("%d triangles\n", currentTriangles.size());
		bigValidationTest();
#endif
	}

	void calculateTriangles(VoronoiNode** list, const int listCount) {
		for (int i = 0; i < listCount; ++i) {
			createTriangulationsFor(list[i], NULL);
		}
	}

	/** @param count how many random points to make in this diagram's boundary */
	void makeRandom(int count) {
		V2f center = boundary->getCenter();
		float rad = boundary->getRadius();
		V2f r(rad, rad);
		V2f min = center - r;
		V2f max = center + r;
		V2f d = max - min;
		for (int i = 0; i < count; ++i) {
			float x = Random::PRNGf(), y = Random::PRNGf();
			V2f randomPoint = V2f(x*d.x + min.x, y*d.y + min.y);
			if (boundary->contains(randomPoint)) {
				currentNodes.add(VoronoiNode(randomPoint));
			} else {
				--i;
			}
		}
	}

	void gatherTrianglesContaining(V2f const & point, TemplateSet<Triangulation*> & out_triangleSet) {
		for (int i = 0; i < currentTriangles.size(); ++i) {
			Triangulation * t = &currentTriangles[i];
			if (t->circumscriptionContains(point)) {
				out_triangleSet.add(t);
			}
		}
	}

	TemplateVector<VoronoiNode*> _p;
	TemplateVector<Triangulation*> _t;

	int indexOf(const Triangulation * t, TemplateVector<const Triangulation*> & triangleList) {
		for (int i = 0; i < triangleList.size(); ++i) {
			if (triangleList[i]->equals(*t)) return i;
		}
		return -1;
	}

	bool findDuplicateTriangle(int startIndex, int & a, int & b) {
		Triangulation * ta, *tb;
		for (a = startIndex; a < currentTriangles.size(); ++a) {
			ta = &currentTriangles[a];
			if (ta->isValidTriangle()) {
				for (b = a + 1; b < currentTriangles.size(); ++b) {
					tb = &currentTriangles[b];
					if (tb->isValidTriangle() && ta->equals(*tb)) {
						return true;
					}
				}
			}
		}
		return false;
	}

	void dupTest() {
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

	void bigValidationTest() {
#ifndef NO_TEST
		// make sure all deleted triangles and edges are marked as such
		for (int i = 0; i < currentTriangles.size(); ++i) {
			Triangulation * t = &currentTriangles[i];
			bool valid = t->isValidTriangle();
			bool inFreeSet = freeTriangles.has(t);
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
		for (int i = 0; i < currentEdges.size(); ++i) {
			Edge * e = &currentEdges[i];
			bool valid = e->isValid();
			bool inFreeSet = freeEdges.has(e);
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
		for (int i = 0; i < currentNodes.size(); ++i) {
			VoronoiNode * vn = &currentNodes[i];
			for (int n = 0; n < vn->getEdgeCount(); ++n){
				Edge * e = vn->getEdge(n);
				if (!e->isValid()) {
					printf("fail. node pointing at invalid edge."); { ABORT }
				}
			}
		}
		// check if every VALID triangle's nodes agree that they are in the same triangle
		for (int i = 0; i < currentTriangles.size(); ++i) {
			Triangulation * t = &currentTriangles[i];
			if (!t->isValidTriangle()) continue;
//			if (t->edges.size() != 3) { printf("fail. triangles must have 3 edges."); ABORT }
			for (int j = 0; j < t->getEdgeCount(); ++j) {
				Edge * e = t->getEdge(j);
				if (!e->isNeighbor(t)) {
					printf("fail. edge is not back-connecting to the triangle (%08x).\n", t);
					printf("edge %08x connects to node %08x and %08x\n", e, e->n0, e->n1);
					for (int i = 0; i < e->getNeighborTriangulationCount(); ++i) {
						printf("neighbor tri %08x\n", e->getNeighborTriangulation(i));
					}
					ABORT
				}
				VoronoiNode * n0 = e->n0, *n1 = e->n1;
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

	AbstractGraphNode * createNode() {
		VoronoiNode * newNode;
		if (freeNodes.size() > 0) {
			newNode = freeNodes.pop();
		} else {
			currentNodes.add();
			newNode = &currentNodes.getLast();
		}
		return newNode;
	}

	bool addNode(V2f point) {
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

	bool addNode(V2f point, TemplateSet<VoronoiNode*> & changedNodes) {
#ifndef NO_TEST
		bigValidationTest();
#endif
		// if this node is not in the bondaries of this diagram, ignore it
		if (boundary && !boundary->contains(point)) return false;
		// if this node already exists, ignore it
		for (int i = 0; i < currentNodes.size(); ++i) {
			if (currentNodes[i].getLocation() == point){
				return false;
			}
		}
		VoronoiNode * newNode = (VoronoiNode*)createNode();
		newNode->set(point);

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

	V2f specialCursor;

	/**
	 * @param location where to look for nodes
	 * @param out_allNodesHere a list to put all of the discovered nodes in. If NULL, the function will stop when the first node is found
	 * @param ignoreList what nodes to ignore
	 * @param ignoreListCount how many nodes to ignore
	 * @return a node found in the given location
	 */
	VoronoiNode * getNodeAt(CircF const & location, TemplateSet<VoronoiNode*> * out_allNodesHere, VoronoiNode** ignoreList, const int ignoreListCount) {
		VoronoiNode * n, * foundOne = NULL;
		for (int i = 0; i < currentNodes.size(); ++i) {
			n = &currentNodes[i];
			if (ignoreList && TemplateArray<VoronoiNode*>::indexOf(n, ignoreList, 0, ignoreListCount) >= 0) continue;
			if (n->isValidNode() && location.contains(n->getLocation())) {
				foundOne = n;
				if (!out_allNodesHere) break;
				out_allNodesHere->add(n);
			}
		}
		return foundOne;
	}

	VoronoiNode * getNodePolyhedronContains(V2f const & point) {
		for (int i = 0; i < currentNodes.size(); ++i) {
			VoronoiNode * n = &currentNodes[i];
			if (n->isValidNode() && n->polyhedronContains(point)) {
				return n;
			}
		}
		return NULL;
	}

	void gatherVoronoi(TemplateVector<VoronoiNode*> & nodes, bool includeBorderPolygons) {
		for (int i = 0; i < currentNodes.size(); ++i) {
			VoronoiNode * vn = &currentNodes[i];
			if (vn->isValidNode()) {
				vn->getPolygon2f(boundary);
				if (includeBorderPolygons || !vn->isBorderPolygon()) {
					nodes.add(vn);
				}
			}
		}
	}

	void glDraw(GLUTRenderingContext & g_screen) const {
		if (boundary) boundary->glDraw(false);

#ifndef NO_TEST
		for (int i = 0; i < currentTriangles.size(); ++i) {
			if (!currentTriangles[i].isValidTriangle()) continue;
			glColor3f(1.0f, .9f, 1.0f);
			currentTriangles[i].circum.glDraw(false);
			glColor3f(1.0f, .0f, 1.0f);
			g_screen.printf(currentTriangles[i].centerMass, "%d", i);
		}
#endif

		for (int i = 0; i < currentEdges.size(); ++i) {
			const Edge * e = &currentEdges[i];
			if (e->isValid()) {
				if (e->getCost() > 0) {
					glColor3f(0.5f, 1.0f, 0.5f);
					e->glDraw();
#ifndef NO_TEST
					// draw which triangles are neighboring the edges. useful for debugging, if you think triangles might be overlapping
					glColor3f(0, 0.7f, 0);
					VoronoiNode * a = (VoronoiNode*)e->getNode(0);
					VoronoiNode * b = (VoronoiNode*)e->getNode(1);
					V2f c = V2f::between(a->getLocation(), b->getLocation());
					g_screen.printf(c, "%d %.1f", e->getNeighborTriangulationCount(), e->getCost());
					for (int n = 0; n < e->getNeighborTriangulationCount(); ++n) {
						Triangulation * t = e->getNeighborTriangulation(n);
						V2f delta = t->centerMass - c;
						g_screen.drawLine(c, c + delta*0.9f);
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
				g_screen.setColor(0xdddddd);
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
			g_screen.setColor(0xff00ff);
			for (int i = 0; i < selectedNode->getEdgeCount(); ++i) {
				const Edge * e = (Edge*)selectedNode->getEdge(i);
#ifndef NO_TEST
				if (!e->isValid()) {
					ABORT
				}
#endif
				e->glDraw();
			}
		}
		const Triangulation * selectedTriangle = NULL;
		for (int i = 0; i < currentTriangles.size(); ++i) {
			const Triangulation* t = &currentTriangles[i];
			if (t->isValidTriangle() && t->polygonContains(specialCursor)) {
				selectedTriangle = t;
			}
		}

		glColor3f(0, 0, .9f);
		for (int i = 0; i < currentNodes.size(); ++i) {
			const VoronoiNode * n = &currentNodes[i];
			if (n->isValidNode()) {
				glDrawCircle(n->getLocation(), 0.2f, false);
				g_screen.printf(n->getLocation(), "%d", i);
			}
		}

		TemplateVector<V2f> edgeBorder;
		for (int i = 0; i < currentTriangles.size(); ++i) {
			if (!currentTriangles[i].isValidTriangle()) continue;

			glDrawCircle(currentTriangles[i].circum.center, .05f, false);

			const Triangulation * t = &currentTriangles[i];

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
					glColor3f(.9f, 0, 0);
					g_screen.drawCircle(edgeCenter, .1f, false);
					glColor3f(0, 0, .9f);
					V2f ray = e->getFace().normal.perp(), end;//-(delta.perp().normal()), end;
					V2f otherPoint;
					a = (VoronoiNode*)e->getNode(0);
					b = (VoronoiNode*)e->getNode(1);
					V2f::closestPointOnLine(a->getLocation(), b->getLocation(), t->centerMass, otherPoint);
					V2f correctDirection = otherPoint - t->centerMass;
					float alignment = V2f::dot(ray, correctDirection);
					if (alignment < 0) { ray *= -1.0f; }
					if (!boundary) {
						end = edgeCenter + ray * 100;
						edgeCenter.glDrawTo(end);
					} else if (boundary->contains(edgeCenter)) {
						
						V2f end, hitNorm;
						float dist;
						if (!boundary->raycast(edgeCenter, ray, dist, end, hitNorm)) {
							end = edgeCenter;
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
				g_screen.setColor(0xff00ff);
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
			g_screen.setColor(selectedNode->isBorderPolygon()?0x0000ff:0xffff00);
			calcedPoly->glDraw(true);
		}
	}

};