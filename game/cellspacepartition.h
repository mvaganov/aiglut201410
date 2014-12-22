#pragma once

#include "codegiraffe/templatevectorlist.h"
#include "obstacles.h"

class Cell_I : public Obstacle {
public:
	virtual void gatherAt(const RectF location, TemplateSet<Obstacle*> & result) const = 0;
	virtual void gatherAt(const CircF location, TemplateSet<Obstacle*> & result) const = 0;
	virtual void gatherAt(const V2f position, TemplateSet<Obstacle*> & result) const = 0;
	virtual void gatherCollisions(Obstacle * subject, TemplateSet<Obstacle*> & out_list) const = 0;
	virtual void add(Obstacle* s) = 0;
	virtual void clear() = 0;
	virtual TemplateVectorList<Obstacle*> * getList() = 0;
	virtual ~Cell_I(){}
	virtual bool raycastContainer(Ray const & ray, RaycastHit & out_rh, float maxDistance, bool dontCareAboutObstacle, Obstacle * & out_hit,
		Obstacle ** ignoreList, int ignoreCount) const = 0;
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
	PartitionCell(RectF area) : aabb(area){ }

	void add(Obstacle * s) { list.add(s); }
	void clear() { list.clear(); }

	void glDraw(bool filled) const {
		aabb.glDraw(filled);
		if (filled) return;
		V2f center = aabb.getCenter();
		for (int i = 0; i < list.size(); ++i){
			center.glDrawTo(list[i]->getShape()->getCenter());
		}
	}

	bool raycastContainer(Ray const & ray, RaycastHit & out_rh, float maxDistance, bool dontCareAboutObstacle, Obstacle * & out_obstacle,
		Obstacle ** ignoreList, int ignoreCount) const {
		RaycastHit best(-1);
		out_obstacle = NULL;
		for (int i = 0; i < list.size(); ++i) {
			if (TemplateArray<Obstacle*>::indexOf(list[i], ignoreList, 0, ignoreCount) < 0) {
				bool hit = list[i]->getShape()->raycast(ray, out_rh);
				if (hit && out_rh.distance < maxDistance) {
					if (dontCareAboutObstacle) return true;
					if (best.distance < 0 || out_rh.distance < best.distance) {
						best = out_rh;
						out_obstacle = list[i];
					}
				}
			}
		}
		if (out_obstacle != NULL) {
			out_rh = best;
		}
		return out_obstacle != NULL;
	}

	void gatherAt(const RectF location, TemplateSet<Obstacle*> & result) const {
		Obstacle* obj;
		for (int i = 0; i < list.size(); ++i){
			obj = list[i];
			if (obj->getShape()->intersectsAABB(location.min, location.max)) { result.add(obj); }
		}
	}
	void gatherAt(const CircF location, TemplateSet<Obstacle*> & result) const {
		Obstacle* obj;
		for (int i = 0; i < list.size(); ++i){
			obj = list[i];
			if (obj->getShape()->intersectsCircle(location.center, location.radius)) { result.add(obj); }
		}
	}
	void gatherAt(const V2f position, TemplateSet<Obstacle*> & result) const {
		Obstacle* obj;
		for (int i = 0; i < list.size(); ++i){
			obj = list[i];
			if (obj->getShape()->contains(position)) { result.add(obj); }
		}
	}
	void gatherCollisions(Obstacle * subject, TemplateSet<Obstacle*> & out_list) const {
		for (int i = 0; i < list.size(); ++i) {
			if (subject != list[i] && subject->getShape()->intersects(list[i]->getShape())) {
				out_list.add(list[i]);
			}
		}
	}
};

/** manages a grid-arranged array of cells (which are either PartitionCell or CellSparePartition objects) */
class CellSpacePartition : public Cell_I {
protected:
	ShapeAABB aabb;
public:
	Shape * getShape(){ return &aabb; }
	const Shape * getShape() const { return &aabb; }
	TemplateVectorList<Obstacle*> * getList() { return NULL; }
	/** the sub-cells of this partition */
	TemplateArray<Cell_I*> cells;
	/** how to divide up this cell into sub-cells */
	V2i gridSize;
	/** how big each sub-cell is */
	V2f cellDimensions;

	CellSpacePartition(RectF totalArea, V2i gridSize) : aabb(totalArea) {
		this->gridSize.set(gridSize);
		cellDimensions.set(aabb.getWidth() / gridSize.x, aabb.getHeight() / gridSize.y);
		V2f cursor = aabb.min;
		RectF cursorR;
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
	void gatherCellListFor(Obstacle * s, TemplateVector<int> & list) const {
		RectF location(s->getShape()->getCenter(), s->getShape()->getRadius());
		RectF rect = getGridIndexeRange(location);
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

	void gatherCellListFor(RectF location, TemplateVector<int> & list) const {
		RectF rect = getGridIndexeRange(location);
		list.ensureCapacity(list.size() + ((int)((rect.getWidth() + 1)*(rect.getHeight() + 1))));
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				list.add(getIndex(r, c));
			}
		}
	}

	void gatherCollisions(Obstacle * subject, TemplateSet<Obstacle*> & out_list) const {
//		printf("\n%s %f,%f:%f", subject->getTypeName(), subject->getCenter().x, subject->getCenter().y, subject->getRadius());
		RectF location(subject->getShape()->getCenter(), subject->getShape()->getRadius());
		RectF rect = getGridIndexeRange(location);
		int index;
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				index = getIndex(r, c);
//				printf("%d.%d:%d ",r,c,index);
				if (subject->getShape()->intersects(cells[index]->getShape())) {
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
		Obstacle ** ignoreList, int ignoreCount) const {
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
			if (!isOutOfBounds && (cells[index]->getShape()->contains(ray.start)
				|| (cells[index]->getShape()->raycast(ray, out_rh) && out_rh.distance < maxDistance))) {
				hitInCell = cells[index]->raycastContainer(ray, out_rh, maxDistance, dontCareAboutObstacle, out_obstacle, ignoreList, ignoreCount);
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

	RectF getGridIndexeRange(RectF location) const {
		//return RectF(getGridIndex(location.min), getGridIndex(location.max));
		const ShapeAABB * const s = &aabb;
		RectF rowcol(location);
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
	RectF _rect = getGridIndexeRange(locationRect);\
	int index;\
	for (int r = (int)_rect.min.y; r <= _rect.max.y; ++r) {\
		for (int c = (int)_rect.min.x; c <= _rect.max.x; ++c) {\
			index = getIndex(r, c);\
			if (index < 0 || index >= cells.size()) { int i = 0; i = 1 / i; }\
			CODE \
		}\
	}

	void gatherAt(const RectF location, TemplateSet<Obstacle*> & result) const {
		VALID_cells_index(location, {
			if (cells[index]->getShape()->intersectsAABB(location.min, location.max)) {
				cells[index]->gatherAt(location, result);
			}
		})
	}
	void gatherAt(const V2f position, TemplateSet<Obstacle*> & result) const {
		VALID_cells_index(RectF(position, 0), {
			if (cells[index]->getShape()->contains(position)) {
				cells[index]->gatherAt(position, result);
			}
		})
	}
	void gatherAt(const CircF circle, TemplateSet<Obstacle*> & result) const {
		VALID_cells_index(RectF(circle.center, circle.radius), {
			if (cells[index]->getShape()->intersectsCircle(circle.center, circle.radius)) {
				cells[index]->gatherAt(circle, result);
			}
		})
	}
#undef VALID_cells_index

	void add(Obstacle* s) {
		RectF location(s->getShape()->getCenter(), s->getShape()->getRadius());
		RectF rect = getGridIndexeRange(location);
		int index;
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				index = getIndex(r, c);
				if (index < 0 || index >= cells.size()) {
					printf("B@D TH!NG T#4T H4PP3ND ");
					continue;
				}
				if (s->getShape()->intersects(cells[index]->getShape())) {
					cells[index]->add(s);
				}
			}
		}
	}
	void clear() {
		for (int i = 0; i < cells.size(); ++i)
			cells[i]->clear();
	}
	bool has(RectF area) const {
		return this->getShape()->intersectsAABB(area.min, area.max);
	}

	void glDraw() {
		V2f cursor = aabb.min;
		V2f line(aabb.getWidth(), 0);
		for (int row = 0; row < gridSize.y-1; ++row) {
			cursor.y += cellDimensions.y;
			cursor.glDrawTo(cursor + line);
		}
		line.set(0, aabb.getHeight());
		for (int col = 0; col < gridSize.x-1; ++col) {
			cursor.x += cellDimensions.x;
			cursor.glDrawTo(cursor + line);
		}
		aabb.glDraw(false);
	}
};
