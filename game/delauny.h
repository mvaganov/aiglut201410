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
#include "cellspacepartition.h"
#include "mempool.h"
#include "obstacles.h"

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
	class VoronoiFace : public Obstacle {
	public:
		// TODO replace points, center, normal, radius with shape.
		ShapePolygon shape;
		/** the points that determine the ends of the separating edge. in a 3D shape these will be the coordinates of the polygon face between voronoi nodes */
		TemplateVector<V2f> points;
		/** where the separating face is, and it's normal */
		V2f center, normal;
		float radius;

		VoronoiFace();

		void refresh(Edge * e);
		void draw(GLUTRenderingContext * g) const;
		/** positive value on one side, negative on the other. not sure which is which. */
		float sideValue(V2f const & p);

		void gatherOppositePoints(V2f & out_a, V2f & out_b);

		/** allows a shape to ignore shapes that do not share bits in the mask */
		long getMask() const { return Obstacle::VORONOI_EDGE; }
		/** allows a shape to ignore shapes that do not share bits in the mask */
		void setMask(const long a_mask) { int i = 0; i = 1 / i; }
		/** binary-OR the given mask into this mask */
		void ensureMaskOverlaps(const int otherMask) { int i = 0; i = 1 / i; }

		Shape * getShape() { return &shape; }
		const Shape * getShape() const { return &shape; }
	};

	/** the graph edge of the voronoi node, used to calculate the delauny triangulation */
	class Edge : public GraphEdge, public Obstacle {
		/** the triangulations that neighbor this edge */
		TemplateSet<Triangulation*> neighborTri;
		/** the face that was calculated for this edge */
		VoronoiFace voronoiFace;
		ShapeLineP shape;
	public:
		VoronoiFace & getFace();

		int getNeighborTriangulationCount() const;

		Triangulation* getNeighborTriangulation(const int i) const;

		void refreshVoronoiFace();

		void addNeighborTriangulation(Triangulation * nt);
		bool removeNeighborTriangulation(Triangulation * nt);

		void invalidateEdge();
		bool isValid() const;
		Edge();
		Edge(VoronoiNode *a, VoronoiNode *b);
		void draw(GLUTRenderingContext * g) const;
		bool isNeighbor(Triangulation * const t) const;

		/** allows a shape to ignore shapes that do not share bits in the mask */
		long getMask() const { return Obstacle::VORONOI_EDGE; }
		/** allows a shape to ignore shapes that do not share bits in the mask */
		void setMask(const long a_mask) { int i = 0; i = 1 / i; }
		/** binary-OR the given mask into this mask */
		void ensureMaskOverlaps(const int otherMask) { int i = 0; i = 1 / i; }

		Shape * getShape() { return &shape; }
		const Shape * getShape() const { return &shape; }
	};

	class Triangulation : public ObstacleConcrete {
	private:
		/** not a set because the order matters. edges are stored in clockwise fashion. */
		TemplateArray<Edge*> edges;
		/** cached calculation of nodes in this Triangulation */
		TemplateSet<VoronoiNode*> nodeSet;
		long mask;
	public:
		Shape * getShape() { return &circum; }
		const Shape * getShape() const { return &circum; }
		/** circumscription of all edges */
		ShapeCircle circum;
//		Circf circum;
		/** the average point of all of the node locations */
		V2f centerMass;
		/** which node to start drawing from */
		VoronoiNode * startingNode;
		void drawEdges(GLUTRenderingContext * g) const;

		Triangulation & operator=(Triangulation const & t);

		void set(Edge* const * list, const int listCount, Circf const & circumscription);

		Triangulation(Edge* const * list, const int listCount, Circf const & circumscription);

		/** create an invalid triangulation by default */
		Triangulation();

		int getEdgeCount() const;
		Edge * getEdge(const int i) const;

		/**
		 * calculate center of mass, and order edges in clockwise fashion.
		 * @return true if this is a valid Triangulation
		 */
		bool calculate();

		/** @return if this is a valid triangulation */
		bool isValidTriangle() const;

		/** make this an invalid triangulation. put any edges freed (nolonger referenced by any triangles) into the given set */
		void invalidate(MemPool<Edge> & freeEdges);

		/** put all of the nodes from referenced edges into the given set */
		void gatherNodesInto(TemplateSet<VoronoiNode*> & nodes) const;

		const TemplateSet<VoronoiNode*> & getNodeSet() const;

		/** only works for 2D triangulation */
		void gatherPointsClockwise(TemplateArray<V2f> & points) const;

		/** only works for 2D triangulation */
		bool polygonContains(V2f const & p) const;

		/** @return true if this triangulation is using the given edge */
		bool hasEdge(const Edge * const edge) const;

		/** @return true if this triangle referenced the given triangle as a neighbor */
		bool isNeighborsWith(Triangulation * t) const;

		/** get the edge that connects the two given nodes. useful when the list of nodes is known but not the edges */
		int findEdgeIndex(VoronoiNode * n0, VoronoiNode * n1) const;

		/** called after a Triangulation is created, to ensure that it's edges are aware of it */
		void recognizeNeighbors();

		/** check this triangulation's circumscription */
		bool circumscriptionContains(V2f const & p) const;

		/** @return if this Triangulation is related to the given node */
		bool hasNode(const VoronoiNode * node) const;
		/** @return if this Triangulation has every node in the given series */
		bool hasAllNodes(VoronoiNode ** list, const int listSize) const;
		/** @return if this node has exactly the nodes in this given series */
		bool hasExactlyTheseNodes(VoronoiNode ** list, const int listSize) const;

		/** a triangulation is equal if it is brokering between the exact same nodes */
		bool equals(const Triangulation & t) const;
		bool operator==(Triangulation const & t) const;
		bool operator!=(Triangulation const & t) const;

		/** allows a shape to ignore shapes that do not share bits in the mask */
		long getMask() const { return mask; }
		/** allows a shape to ignore shapes that do not share bits in the mask */
		void setMask(const long a_mask) { mask = a_mask; }
		/** binary-OR the given mask into this mask */
		void ensureMaskOverlaps(const int otherMask) { mask |= otherMask; }
	};

	class VoronoiNode : public GraphNode, public Obstacle  {
		ShapePolygon calculatedPolygon;
		ShapePoint point;
		/** identifies if this is on the edge of the voronoi diagram (near the boundaries) */
		Obstacle * atBoundaryOf;
		/** */
		bool valid;
		/** whenever data changes, the node should flag itself as needed recalculation */
		bool needsModelRecalculated;
		long mask;
	public:
		/** allows a shape to ignore shapes that do not share bits in the mask */
		long getMask() const {
			return mask;
		}
		/** allows a shape to ignore shapes that do not share bits in the mask */
		void setMask(const long a_mask) { mask = a_mask; }
		/** binary-OR the given mask into this mask */
		void ensureMaskOverlaps(const int otherMask) { mask |= otherMask; }

		Shape * getShape() { return &point; }
		const Shape * getShape() const { return &point; }
		VoronoiNode();
		VoronoiNode(V2f const & p);

		V2f & getLocation();
		const V2f & getLocation() const;

		/** @return whether the polygon needs to be recalculated */
		bool isDirty();
		/** identify that the polygon needs to be recalculated */
		void setDirty();

		bool isBorderPolygon() const;

		int getEdgeCount() const;
		//		Edge * getEdge(const int i) const { return edges[i]; }
		bool removeEdge(Edge * e);
		//		bool addEdge(Edge * e) { setDirty(); return edges.add(e); }

		void set(V2f const & p);

		bool isValidNode() const;

		void invalidate(DelaunySet & D, TemplateSet<VoronoiNode*> & changedNodes);

		/** calculate the Voronoi polygon*/
		void calculatePolygon(TemplateVector<V2f> & out_points, V2f & out_center) const;

		/** calculate the Voronoi polygon*/
		Polygon2f & getPolygon2f(Obstacle * boundary);

		bool polyhedronContains(V2f const & a_position) const;

		/** @return the edge from this node to the given node */
		Edge * getEdgeTo(VoronoiNode * node) const;

		/** gather a list of every triangulation attached to this node */
		void gatherTriangles(TemplateSet<Triangulation*> & triangles) const;

		void gatherEdges(TemplateSet<Edge*> & edges) const;
	};

	/**
	* @param nodes what nodes to calculate the clockwise order of
	* @param nodeCount how many nodes there are
	* @param out_centerMass the average location of the nodes
	* @param out_inOrder where to put the list re-arranged in order
	*/
	static void orderNodesClockwise(VoronoiNode*const* nodes, const int nodeCount, V2f  & out_centerMass, VoronoiNode** out_inOrder);

	int getNodeCount() const;
	AbstractGraphNode * getNode(const int index);

	void connectNodes(AbstractGraphNode * from, AbstractGraphNode * to);
	AbstractGraphEdge * getEdge(AbstractGraphNode * a, AbstractGraphNode * b, bool createIfNotThere);

	MemPool<VoronoiNode> m_nodePool;
	MemPool<Edge> edgePool;
	MemPool<Triangulation> triangulationPool;

	CellSpacePartition * csp;

	VoronoiNode * selectedNode;
	Obstacle * boundary;

	DelaunySet(Obstacle * boundary);

	~DelaunySet();

	Triangulation * addTriangle(Triangulation const & t);

	void removeTriangle(Triangulation * toRemove);

	void destroyNode(AbstractGraphNode * n);

	void removeNode(VoronoiNode * toRemove);

	void removeNode(VoronoiNode * toRemove, TemplateSet<VoronoiNode*> & changedNodes);

	void moveNode(VoronoiNode * toMove, V2f newPosition);

	Edge * findEdge(VoronoiNode * a, VoronoiNode * b);

	Edge* addEdge(Edge const & e);

	Edge* marshalEdge(VoronoiNode * a, VoronoiNode * b);

	Triangulation * getTriangle(TemplateSet<VoronoiNode*> & nodeCluster);

	/** call this when you know that this is a valid triangulation... */
	Triangulation * marshalTriangulation(TemplateSet<VoronoiNode *> & nodeCluster, Circf circumscription);

	/** @return true if the given node cluster can be triangulated into a circumscription */
	bool triangulateFor(TemplateSet<VoronoiNode*> & nodeCluster, float floatingPointRounding, Circf & out_circ);

	/**
	 * @param node what to create triangulations around
	 * @param createdTriangles if not null, adds all new triangles to this set
	 */
	void createTriangulationsFor(VoronoiNode* node, TemplateSet<Triangulation*> * createdTriangles);

	bool createTriangulationInternal(TemplateSet<VoronoiNode*> & nodeCluster, int startIndex, TemplateSet<Triangulation*> * createdTriangles, Circf whereToLookForNodes);

	void calculateAllTriangles();

	void calculateTriangles(VoronoiNode** list, const int listCount);

	/** @param count how many random points to make in this diagram's boundary */
	void makeRandom(int count);

	void gatherTrianglesContaining(V2f const & point, TemplateSet<Triangulation*> & out_triangleSet);

	TemplateVector<VoronoiNode*> _p;
	TemplateVector<Triangulation*> _t;

	int indexOf(const Triangulation * t, TemplateVector<const Triangulation*> & triangleList);

	bool findDuplicateTriangle(int startIndex, int & a, int & b);

	void dupTest();

	void bigValidationTest();

	AbstractGraphNode * createNode();

	bool addNode(V2f point);

	bool addNode(V2f point, TemplateSet<VoronoiNode*> & changedNodes);

	V2f specialCursor;

	/**
	 * @param location where to look for nodes
	 * @param out_allNodesHere a list to put all of the discovered nodes in. If NULL, the function will stop when the first node is found
	 * @param ignoreList what nodes to ignore
	 * @param ignoreListCount how many nodes to ignore
	 * @return a node found in the given location
	 */
	VoronoiNode * getNodeAt(Circf const & location, TemplateSet<VoronoiNode*> * out_allNodesHere, VoronoiNode** ignoreList, const int ignoreListCount);
	Obstacle * DelaunySet::getNodeAt(Circf const & location, TemplateSet<Obstacle*> & out_allNodesHere);

	//VoronoiNode * getNodePolyhedronContains(V2f const & point);

	//void gatherVoronoi(TemplateVector<VoronoiNode*> & nodes, bool includeBorderPolygons);

	void draw(GLUTRenderingContext * g) const;
};