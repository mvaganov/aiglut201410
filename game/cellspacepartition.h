#pragma once

#include "codegiraffe/templatevectorlist.h"
#include "obstacles.h"
#include "codegiraffe/glutrenderingcontext.h"

class Cell_I : public Obstacle {
public:
	virtual void gatherAtAABB(const AABBf location, TemplateSet<Obstacle*> & result, const long mask) const = 0;
	virtual void gatherAtCircle(const Circf location, TemplateSet<Obstacle*> & result, const long mask) const = 0;
	virtual void gatherAtPosition(const V2f position, TemplateSet<Obstacle*> & result, const long mask) const = 0;
	virtual void gatherAtLayer(TemplateSet<Obstacle*> & result, const long mask) const = 0;
	/** mask parameter not included because each Obstacle already has a mask to use. */
	virtual void gatherCollisions(Obstacle * subject, TemplateSet<Obstacle*> & out_list) const = 0;
	virtual void add(Obstacle* s) = 0;
	virtual bool remove(Obstacle * s) = 0;
	virtual void removeObstaclesMarked(const long mask) = 0;
	virtual void clear() = 0;
	virtual TemplateVectorList<Obstacle*> * getList() = 0;
	virtual ~Cell_I(){}
	virtual bool raycastContainer(Ray const & ray, RaycastHit & out_rh, float maxDistance, bool dontCareAboutObstacle, Obstacle * & out_hit,
		const long mask) const = 0;
		//Obstacle ** ignoreList, int ignoreCount) const = 0;
};

/** the end of a cellspace partition. only holds a list of elements. */
class PartitionCell : public Cell_I {
protected:
	ShapeAABB aabb;
	TemplateVectorList<Obstacle*> list;
public:
	Shape * getShape(){ return &aabb; }
	const Shape * getShape() const { return &aabb; }
	TemplateVectorList<Obstacle*> * getList() { return &list; }
	PartitionCell() {}
	PartitionCell(AABBf area) : aabb(area){ }

	void add(Obstacle * s) { ensureMaskOverlaps(s->getMask()); list.add(s); }
	bool remove(Obstacle * s) { return list.removeDataFast(s); }
	void removeObstaclesMarked(const long mask) {
		for (int i = list.size() - 1; i >= 0; --i) {
			if (list[i]->masksCollide(mask)) {
				list.removeFast(i);
			}
		}
	}
	void clear() { list.clear(); }

	void draw(GLUTRenderingContext * g, bool filled) const {
		g->drawRect(aabb, filled);
		if (filled) return;
		V2f center = aabb.getCenter();
		for (int i = 0; i < list.size(); ++i){
			g->drawLine(center, list[i]->getShape()->getCenter());
		}
	}

	bool raycastContainer(Ray const & ray, RaycastHit & out_rh, float maxDistance, bool dontCareAboutObstacle, Obstacle * & out_obstacle,
		//Obstacle ** ignoreList, int ignoreCount) const {
		const long mask) const {
		RaycastHit best(-1);
		out_obstacle = NULL;
		for (int i = 0; i < list.size(); ++i) {
			//if (TemplateArray<Obstacle*>::indexOf(list[i], ignoreList, 0, ignoreCount) < 0) {
				//bool hit = list[i]->getShape()->raycast(ray, out_rh);
				bool hit = list[i]->raycast(ray, out_rh, mask);
				if (hit && out_rh.distance < maxDistance) {
					if (dontCareAboutObstacle) return true;
					if (best.distance < 0 || out_rh.distance < best.distance) {
						best = out_rh;
						out_obstacle = list[i];
					}
				}
			//}
		}
		if (out_obstacle != NULL) {
			out_rh = best;
		}
		return out_obstacle != NULL;
	}

	void gatherAtAABB(const AABBf location, TemplateSet<Obstacle*> & result, const long mask) const {
		Obstacle* obj;
		for (int i = 0; i < list.size(); ++i){
			obj = list[i];
			//if (obj->getShape()->intersectsAABB(location.min, location.max)) { result.add(obj); }
			if (obj->intersectsAABB(location.min, location.max, mask)) { result.add(obj); }
		}
	}
	void gatherAtCircle(const Circf location, TemplateSet<Obstacle*> & result, const long mask) const {
		Obstacle* obj;
		for (int i = 0; i < list.size(); ++i){
			obj = list[i];
			//if (obj->getShape()->intersectsCircle(location.center, location.radius)) { result.add(obj); }
			if (obj->intersectsCircle(location.center, location.radius, mask)) { result.add(obj); }
		}
	}
	void gatherAtPosition(const V2f position, TemplateSet<Obstacle*> & result, const long mask) const {
		Obstacle* obj;
		for (int i = 0; i < list.size(); ++i){
			obj = list[i];
			//if (obj->getShape()->contains(position)) { result.add(obj); }
			if (obj->contains(position, mask)) { result.add(obj); }
		}
	}
	void gatherAtLayer(TemplateSet<Obstacle*> & result, const long mask) const {
		Obstacle* obj;
		for (int i = 0; i < list.size(); ++i){
			obj = list[i];
			if (obj->masksCollide(mask)) { result.add(obj); }
		}
	}

	void gatherCollisions(Obstacle * subject, TemplateSet<Obstacle*> & out_list) const {
		for (int i = 0; i < list.size(); ++i) {
			//if (subject != list[i] && subject->getShape()->intersects(list[i]->getShape())) {
			if (subject != list[i] && subject->intersects(list[i]->getShape(), list[i]->getMask())) {
				out_list.add(list[i]);
			}
		}
	}
};

/** manages a grid-arranged array of cells (which are either PartitionCell or CellSparePartition objects) */
class CellSpacePartition : public Cell_I {
protected:
public:
	ShapeAABB aabb;
	Shape * getShape(){ return &aabb; }
	const Shape * getShape() const { return &aabb; }
	TemplateVectorList<Obstacle*> * getList() { return NULL; }
	/** the sub-cells of this partition */
	TemplateArray<Cell_I*> cells;
	/** how to divide up this cell into sub-cells */
	V2i gridSize;
	/** how big each sub-cell is */
	V2f cellDimensions;

	CellSpacePartition(AABBf totalArea, V2i gridSize) : aabb(totalArea) {
		this->gridSize.set(gridSize);
		cellDimensions.set(aabb.getWidth() / gridSize.x, aabb.getHeight() / gridSize.y);
		V2f cursor = aabb.min;
		AABBf cursorR;
		int index = 0;
		cells.allocateToSize(gridSize.x*gridSize.y);
		for (int r = 0; r < gridSize.y; ++r)
		{
			for (int c = 0; c < gridSize.x; ++c)
			{
				cursorR.set(cursor, cursor + cellDimensions);
				cells[index] = new PartitionCell(cursorR);
				index++;
				cursor.x += cellDimensions.x;
			}
			cursor.x = aabb.min.x;
			cursor.y += cellDimensions.y;
		}
	}
	~CellSpacePartition(){ for (int i = 0; i < cells.size(); ++i) delete cells[i]; }
private:
	int getIndex(int row, int col) const { return row*gridSize.x + col; }
public:
	void gatherCellListFor(Obstacle * s, TemplateVector<int> & list) const { // TODO is this used? remove if not.
		AABBf location(s->getShape()->getCenter(), s->getShape()->getRadius());
		AABBf rect = getGridIndexeRange(location);
		int index;
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				index = getIndex(r, c);
				if (index < 0 || index >= cells.size()) {
					int i = 0; i = 1 / i;
					continue;
				}
				if (s->getShape()->intersects(cells[index]->getShape())) {
					list.add(index);
				}
			}
		}
	}

	void gatherCellListFor(AABBf location, TemplateVector<int> & list) const { // TODO is this used? remove if not.
		AABBf rect = getGridIndexeRange(location);
		list.ensureCapacity(list.size() + ((int)((rect.getWidth() + 1)*(rect.getHeight() + 1))));
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				list.add(getIndex(r, c));
			}
		}
	}

	void gatherCollisions(Obstacle * subject, TemplateSet<Obstacle*> & out_list) const {
//		printf("\n%s %f,%f:%f", subject->getTypeName(), subject->getCenter().x, subject->getCenter().y, subject->getRadius());
		AABBf location(subject->getShape()->getCenter(), subject->getShape()->getRadius());
		AABBf rect = getGridIndexeRange(location);
		int index;
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				index = getIndex(r, c);
//				printf("%d.%d:%d ",r,c,index);
				if (cells[index]->intersects(subject, subject->getMask())) {//  subject->getShape()->intersects(cells[index]->getShape())) {
					cells[index]->gatherCollisions(subject, out_list);
				}
			}
		}
	}

	// used to test the minimal-cell-testing raycasting algorithm
	void raycastCellList(Ray ray, float maxDistance, TemplateSet<int> & cellIDs) const {
		RaycastHit out_rh;
		V2f point = ray.start;
		if (!aabb.contains(point)) {
			if (!aabb.raycast(ray, out_rh)){
				return;
			}
			point = out_rh.point;
		}
		V2f rowcol = getGridIndex(point);
		// when this line steps on the grid, how does it step?
		V2i step(
			((ray.direction.x > 0) ? 1 : ((ray.direction.x < 0) ? -1 : 0)),
			((ray.direction.y > 0) ? 1 : ((ray.direction.y < 0) ? -1 : 0))
			);
		const int NUM_DIMENSIONS = 2;
		int orderOfMoves[NUM_DIMENSIONS];
		// weight the dimensions stepped
		if (abs(ray.direction.x) > abs(ray.direction.y)) {
			orderOfMoves[0] = 0;
			orderOfMoves[1] = 1;
		} else {
			orderOfMoves[0] = 1;
			orderOfMoves[1] = 0;
		}
		// hwo the line is moving on the grid (which dimension)
		int whichMoveToTake = 0;
		// 0 for X, 1 for Y
		int field;
		bool raycasting = true;
		while (raycasting)
		{
			bool isOutOfBounds = rowcol.x < 0 || rowcol.y < 0 || rowcol.x >= gridSize.x || rowcol.y >= gridSize.y;
			int index = getIndex((int)rowcol.y, (int)rowcol.x);
			if (!isOutOfBounds && (cells[index]->getShape()->contains(ray.start)
			|| (cells[index]->getShape()->raycast(ray, out_rh) && out_rh.distance < maxDistance))) {
				cellIDs.add(index);
				if (step.isZero()) break;
				// next step
				whichMoveToTake = 0;
				field = orderOfMoves[whichMoveToTake];
				rowcol.setField(field, rowcol.getField(field) + step.getField(field));
			} else {
				// clearly that last step didn't work so well. lets try another one
				field = orderOfMoves[whichMoveToTake];
				rowcol.setField(field, rowcol.getField(field) - step.getField(field));
				whichMoveToTake++;
				if (whichMoveToTake >= NUM_DIMENSIONS) break; // if there are no more steps, finished!
				field = orderOfMoves[whichMoveToTake];
				if (step.getField(field) == 0) break; // if the next step is a non-step, finished!
				rowcol.setField(field, rowcol.getField(field) + step.getField(field));
			}
		}
	}

	bool raycastContainer(Ray const & ray, RaycastHit & out_rh, float maxDistance, bool dontCareAboutObstacle, Obstacle * & out_obstacle,
		const long mask) const {
		//Obstacle ** ignoreList, int ignoreCount) const {
		V2f point = ray.start;
		if (!aabb.contains(point)) {
			if (!aabb.raycast(ray, out_rh)){
				return false;
			}
			point = out_rh.point;
		}
		V2f rowcol = getGridIndex(point);
		// when this line steps on the grid, how does it step?
		V2i step(
			((ray.direction.x > 0) ? 1 : ((ray.direction.x < 0) ? -1 : 0)),
			((ray.direction.y > 0) ? 1 : ((ray.direction.y < 0) ? -1 : 0))
			);
		const int NUM_DIMENSIONS = 2;
		int orderOfMoves[NUM_DIMENSIONS];
		// weight the dimensions stepped based on the direction being travelled
		if (abs(ray.direction.x) > abs(ray.direction.y)) {
			orderOfMoves[0] = 0;
			orderOfMoves[1] = 1;
		} else {
			orderOfMoves[0] = 1;
			orderOfMoves[1] = 0;
		}
		// how the line is moving on the grid (which dimension is stepped)
		int whichMoveToTake = 0;
		// 0 for X, 1 for Y
		int field;
		bool raycasting = true, hitInCell = false;
		while (raycasting)
		{
			bool isOutOfBounds = rowcol.x < 0 || rowcol.y < 0 || rowcol.x >= gridSize.x || rowcol.y >= gridSize.y;
			int index = getIndex((int)rowcol.y, (int)rowcol.x);
			//if (!isOutOfBounds && (cells[index]->getShape()->contains(ray.start)
			//	|| (cells[index]->getShape()->raycast(ray, out_rh) && out_rh.distance < maxDistance))) {
			if (!isOutOfBounds && (cells[index]->contains(ray.start, mask)
				|| (cells[index]->raycast(ray, out_rh, mask) && out_rh.distance < maxDistance))) {
				hitInCell = cells[index]->raycastContainer(ray, out_rh, maxDistance, dontCareAboutObstacle, out_obstacle, mask);
				if (hitInCell || step.isZero()) break;
				// next step
				whichMoveToTake = 0;
				field = orderOfMoves[whichMoveToTake];
				rowcol.setField(field, rowcol.getField(field) + step.getField(field));
			} else {
				// clearly that last step didn't work so well. lets try another one
				field = orderOfMoves[whichMoveToTake];
				rowcol.setField(field, rowcol.getField(field) - step.getField(field));
				whichMoveToTake++;
				if (whichMoveToTake >= NUM_DIMENSIONS) break; // if there are no more steps, finished!
				field = orderOfMoves[whichMoveToTake];
				if (step.getField(field) == 0) break; // if the next step is a non-step, finished!
				rowcol.setField(field, rowcol.getField(field) + step.getField(field));
			}
		}
		return hitInCell;
	}

	V2f getGridIndex(V2f position) const {
		V2f rowcol(position);
		rowcol -= aabb.min;
		rowcol /= cellDimensions;
		if (rowcol.x < 0) { rowcol.x = 0; }
		if (rowcol.y < 0) { rowcol.y = 0; }
		if (rowcol.x >= gridSize.x) { rowcol.x = (float)(gridSize.x - 1); }
		if (rowcol.y >= gridSize.y) { rowcol.y = (float)(gridSize.y - 1); }
		rowcol.clampToInt();
		return rowcol;
	}

	AABBf getGridIndexeRange(AABBf location) const {
		//return AABBf(getGridIndex(location.min), getGridIndex(location.max));
		const ShapeAABB * const s = &aabb;
		AABBf rowcol(location);
		V2f totalSize = s->max - s->min;
		rowcol.min -= s->min;
		rowcol.max -= s->min;
		rowcol.min /= cellDimensions;
		rowcol.max /= cellDimensions;
		if (rowcol.min.x < 0)	{ rowcol.min.x = 0; if (rowcol.max.x < rowcol.min.x) { rowcol.max.x = rowcol.min.x; } }
		if (rowcol.min.y < 0)	{ rowcol.min.y = 0; if (rowcol.max.y < rowcol.min.y) { rowcol.max.y = rowcol.min.y; } }
		if (rowcol.max.x >= gridSize.x)	{ rowcol.max.x = (float)(gridSize.x - 1); if (rowcol.min.x > rowcol.max.x) { rowcol.min.x = rowcol.max.x; } }
		if (rowcol.max.y >= gridSize.y)	{ rowcol.max.y = (float)(gridSize.y - 1); if (rowcol.min.y > rowcol.max.y) { rowcol.min.y = rowcol.max.y; } }
		rowcol.clampToInt();
		return rowcol;
	}
#define VALID_cells_index(locationRect, CODE) \
	AABBf _rect = getGridIndexeRange(locationRect);\
	int index;\
	for (int r = (int)_rect.min.y; r <= _rect.max.y; ++r) {\
		for (int c = (int)_rect.min.x; c <= _rect.max.x; ++c) {\
			index = getIndex(r, c);\
			if (index < 0 || index >= cells.size()) { int i = 0; i = 1 / i; }\
			CODE \
		}\
	}

	void gatherAtAABB(const AABBf location, TemplateSet<Obstacle*> & result, const long mask) const {
		VALID_cells_index(location, {
			//if (cells[index]->getShape()->intersectsAABB(location.min, location.max)) {
			if (cells[index]->intersectsAABB(location.min, location.max, mask)) {
				cells[index]->gatherAtAABB(location, result, mask);
			}
		})
	}
	void gatherAtPosition(const V2f position, TemplateSet<Obstacle*> & result, const long mask) const {
		VALID_cells_index(AABBf(position, 0), {
			//if (cells[index]->getShape()->contains(position)) {
			if (cells[index]->contains(position, mask)) {
				cells[index]->gatherAtPosition(position, result, mask);
			}
		})
	}
	void gatherAtCircle(const Circf circle, TemplateSet<Obstacle*> & result, const long mask) const {
		VALID_cells_index(AABBf(circle.center, circle.radius), {
			//if (cells[index]->getShape()->intersectsCircle(circle.center, circle.radius)) {
			if (cells[index]->intersectsCircle(circle.center, circle.radius, mask)) {
				cells[index]->gatherAtCircle(circle, result, mask);
			}
		})
	}

	void gatherAtLayer(TemplateSet<Obstacle*> & result, const long mask) const {
		for (int i = 0; i < cells.size(); ++i){
			if (cells[i]->masksCollide(mask)) { cells[i]->gatherAtLayer(result, mask); }
		}
	}

#undef VALID_cells_index

	void add(Obstacle* s) {
		this->ensureMaskOverlaps(s->getMask());
		AABBf location(s->getShape()->getCenter(), s->getShape()->getRadius());
		AABBf rect = getGridIndexeRange(location);
		int index;
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				index = getIndex(r, c);
				if (index < 0 || index >= cells.size()) {
					printf("B@D TH!NG T#4T H4PP3ND ");
					continue;
				}
				// when adding, don't test the mask, otherwise nothing will ever get added. it isnt until something is added that a mask exists at all.
				if (s->getShape()->intersects(cells[index]->getShape())) {
					cells[index]->add(s);
				}
			}
		}
	}

	bool remove(Obstacle * s) {
		if (!masksCollide(s)) return false;
		AABBf location(s->getShape()->getCenter(), s->getShape()->getRadius());
		AABBf rect = getGridIndexeRange(location);
		int index;
		bool removed = false;
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				index = getIndex(r, c);
				if (index < 0 || index >= cells.size()) {
					printf("B@D TH!NG T#4T H4PP3ND ");
					continue;
				}
				// when adding, don't test the mask, otherwise nothing will ever get added. it isnt until something is added that a mask exists at all.
				if (s->getShape()->intersects(cells[index]->getShape())) {
					removed |= cells[index]->remove(s);
				}
			}
		}
		return removed;
	}

	void clear() {
		this->setMask(0);
		for (int i = 0; i < cells.size(); ++i)
			cells[i]->clear();
	}
	//bool has(AABBf area) const {
	//	return this->getShape()->intersectsAABB(area.min, area.max);
	//}

	void draw(GLUTRenderingContext * g) {
		// draw grid
		V2f cursor = aabb.min;
		V2f line(aabb.getWidth(), 0);
		for (int row = 0; row < gridSize.y-1; ++row) {
			cursor.y += cellDimensions.y;
			cursor.glDrawTo(cursor + line);
		}
		line.set(0, aabb.getHeight());
		cursor.y = aabb.min.y;
		for (int col = 0; col < gridSize.x-1; ++col) {
			cursor.x += cellDimensions.x;
			cursor.glDrawTo(cursor + line);
		}
		aabb.glDraw(false);

		// draw individual cells
		CellSpacePartition * csp;
		for (int i = 0; i < cells.size(); ++i) {
			cursor = cells[i]->getShape()->getCenter();
			g->printf(cursor, "%d", i);
			csp = dynamic_cast<CellSpacePartition*>(cells[i]);
			if (csp != NULL) {
				csp->draw(g);
			}
			else {
				TemplateVectorList<Obstacle*> * list = cells[i]->getList();
				if (list != NULL) {
					V2f p;
					for (int a = 0; a < list->size(); ++a) {
						p = list->get(a)->getShape()->getCenter();
						g->printf(p, "%d", a);
						g->drawLine(cursor, p);
					}
				}
			}
		}
	}

	void removeObstaclesMarked(const long mask) {
		bool mustCheck = mask == MARKED_FOR_DELETE;
		for (int row = 0; row < gridSize.y; ++row) {
			for (int col = 0; col < gridSize.x; ++col) {
				int index = getIndex(row, col);
				if (mustCheck || cells[index]->masksCollide(mask))
					cells[index]->removeObstaclesMarked(mask);
			}
		}
	}
};

/**
 * used when multiple CellSpacePartitions are appropriate. For example, a separate Static-Obstacle vs. Moving-Obstacle layer
 */
class LayeredPartitions : public Cell_I {
	TemplateVector<CellSpacePartition*> mapLayer;
public:
public:
	ShapeAABB aabb;
	Shape * getShape(){ return &aabb; }
	const Shape * getShape() const { return &aabb; }
	TemplateVectorList<Obstacle*> * getList() { return NULL; }

	void draw(GLUTRenderingContext * g, int * layerColor, int colorCount) {
		for (int i = 0; i < mapLayer.size(); ++i) {
			g->setColor(layerColor[i % colorCount]);
			mapLayer[i]->draw(g);
		}
	}

	LayeredPartitions(AABBf totalArea, V2i defaultGridSize) {
		CellSpacePartition * csp = new CellSpacePartition(totalArea, defaultGridSize);
		csp->setMask(Obstacle::EVERYTHING);
		mapLayer.add(csp);
	}
	void addLayer(const long mask) {
		addLayer(mask, mapLayer[0]->aabb, mapLayer[0]->gridSize);
	}
	void addLayer(const long mask, AABBf totalArea, V2i gridSize) {
		// check if exactly-this layer exists already. if it does, ignore this call.
		for (int i = 0; i < mapLayer.size(); ++i) {
			if (mapLayer[i]->getMask() == mask) return;
		}
		// add this layer, and set the mask
		CellSpacePartition * csp = new CellSpacePartition(totalArea, gridSize);
		csp->setMask(Obstacle::EVERYTHING);
		mapLayer.add(csp);
		// go through layer 0, remove all elements that belong to the new layer, and add them to the new layer
		TemplateSet<Obstacle*> list;
		mapLayer[0]->gatherAtLayer(list, mask);
		mapLayer[0]->removeObstaclesMarked(mask);
		for (int i = 0; i < list.size(); ++i)
			csp->add(list[i]);
	}
#define GATHER_AT_(specificFunction) \
	for (int i = 0; i < mapLayer.size(); ++i) { \
		if (mapLayer[i]->masksCollide(mask)) { \
			mapLayer[i]->specificFunction(location, result, mask); \
		} \
	}

	void gatherAtAABB(const AABBf location, TemplateSet<Obstacle*> & result, const long mask) const { GATHER_AT_(gatherAtAABB) }
	void gatherAtCircle(const Circf location, TemplateSet<Obstacle*> & result, const long mask) const { GATHER_AT_(gatherAtCircle) }
	void gatherAtPosition(const V2f location, TemplateSet<Obstacle*> & result, const long mask) const { GATHER_AT_(gatherAtPosition) }
#undef GATHER_AT_
	void gatherAtLayer(TemplateSet<Obstacle*> & result, const long mask) const {
		for (int i = 0; i < mapLayer.size(); ++i){
			if (mapLayer[i]->masksCollide(mask)) { mapLayer[i]->gatherAtLayer(result, mask); }
		}
	}

	/** mask parameter not included because each Obstacle already has a mask to use. */
	void gatherCollisions(Obstacle * subject, TemplateSet<Obstacle*> & out_list) const {
		for (int i = 0; i < mapLayer.size(); ++i) {
			if (mapLayer[i]->masksCollide(subject->getMask())) {
				mapLayer[i]->gatherCollisions(subject, out_list);
			}
		}
	}
	void add(Obstacle* s) {
		bool added = false;
		// check every layer after layer 0 to see if this object fits in there. if it doesnt, add it to layer 0
		for (int i = 1; i < mapLayer.size(); ++i) {
			if (mapLayer[i]->masksCollide(s->getMask())) {
				mapLayer[i]->add(s);
				added = true;
			}
		}
		if (!added) {
			mapLayer[0]->add(s);
		}
	}

	void add(Obstacle** obs, const int size) {
		for (int i = 0; i < size; ++i) {
			add(obs[i]);
		}
	}

	bool remove(Obstacle * s) {
		bool somethingRemoved = false;
		for (int i = 0; i < mapLayer.size(); ++i) {
			if (mapLayer[i]->masksCollide(s->getMask())) {
				mapLayer[i]->remove(s);
				somethingRemoved = false;
			}
		}
		return somethingRemoved;
	}

	void remove(Obstacle** obs, const int size) {
		for (int i = 0; i < size; ++i) {
			remove(obs[i]);
		}
	}

	void removeObstaclesMarked(const long mask) {
		bool mustCheck = mask == MARKED_FOR_DELETE;
		for (int i = 0; i < mapLayer.size(); ++i) {
			if (mustCheck || mapLayer[i]->masksCollide(mask))
				mapLayer[i]->removeObstaclesMarked(mask);
		}
	}
	void clear(const long mask) {
		for (int i = 0; i < mapLayer.size(); ++i) {
			if (mapLayer[i]->masksCollide(mask))
				mapLayer[i]->clear();
		}
	}
	void clear() {
		for (int i = 0; i < mapLayer.size(); ++i) { mapLayer[i]->clear(); }
	}
	~LayeredPartitions(){
		for (int i = 0; i < mapLayer.size(); ++i) {
			delete mapLayer[i];
		}
	}
	bool raycastContainer(Ray const & ray, RaycastHit & out_rh, float maxDistance, bool dontCareAboutObstacle, Obstacle * & out_hit,
		const long mask) const {
		bool raycastHit = false;
		for (int i = 0; i < mapLayer.size(); ++i) {
			raycastHit |= mapLayer[i]->raycastContainer(ray, out_rh, maxDistance, dontCareAboutObstacle, out_hit, mask);
		}
		return raycastHit;
	}
};
