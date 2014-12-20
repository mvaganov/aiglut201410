#include "codegiraffe/templatevectorlist.h"
#include "obstacles.h"

// TYPE should implement Shape
template<typename TYPE>
class Cell_I : public AABBObject {
public:
	Cell_I(RectF const & r) : AABBObject(r) {}
	virtual void gatherAt(const RectF location, TemplateSet<TYPE*> & result) const = 0;
	virtual void gatherAt(const CircF location, TemplateSet<TYPE*> & result) const = 0;
	virtual void gatherAt(const V2f position, TemplateSet<TYPE*> & result) const = 0;
	virtual void gatherCollisions(TYPE * subject, TemplateSet<TYPE*> & out_list) const = 0;
	virtual void add(TYPE* s) = 0;
	virtual void clear() = 0;
	virtual TemplateVectorList<TYPE*> * getList() = 0;
	virtual ~Cell_I(){}
	virtual bool raycastContainer(V2f start, V2f direction, float maxDistance, bool dontCareAboutObstacle, TYPE * & out_obstacle,
		float & out_dist, V2f & out_point, V2f & out_normal, TYPE ** ignoreList, int ignoreCount) const = 0;
};

template<typename TYPE>
class PartitionCell : public Cell_I<TYPE> {
protected:
	TemplateVectorList<TYPE*> list;
public:
	TemplateVectorList<TYPE*> * getList() { return &list; }
	PartitionCell() : Cell_I(RectF()) {}
	PartitionCell(RectF area) : Cell_I(area){ }

	void add(TYPE * entity) { list.add(entity); }
	void clear() { list.clear(); }

	void glDraw(bool filled) const {
		AABBObject::glDraw(filled);
		if (filled) return;
		V2f center = getCenter();
		for (int i = 0; i < list.size(); ++i){
			center.glDrawTo(list[i]->getCenter());
		}
	}

	bool raycastContainer(V2f start, V2f direction, float maxDistance, bool dontCareAboutObstacle, TYPE * & out_obstacle,
		float & out_dist, V2f & out_point, V2f & out_normal, TYPE ** ignoreList, int ignoreCount) const {
		float bestDist = -1;
		V2f bestPoint, bestNormal;
		out_obstacle = NULL;
		for (int i = 0; i < list.size(); ++i) {
			if (TemplateArray<TYPE*>::indexOf(list[i], ignoreList, 0, ignoreCount) < 0) {
				bool hit = list[i]->raycast(start, direction, out_dist, out_point, out_normal);
				if (hit && out_dist < maxDistance) {
					if (dontCareAboutObstacle) return true;
					if (bestDist < 0 || out_dist < bestDist) {
						bestDist = out_dist;
						bestPoint = out_point;
						bestNormal = out_normal;
						out_obstacle = list[i];
					}
				}
			}
		}
		if (out_obstacle != NULL) {
			out_dist = bestDist;
			out_point = bestPoint;
			out_normal = bestNormal;
		}
		return out_obstacle != NULL;
	}

	void gatherAt(const RectF location, TemplateSet<TYPE*> & result) const {
		TYPE* obj;
		for (int i = 0; i < list.size(); ++i){
			obj = list[i];
			if (obj->intersectsAABB(location.min, location.max)) { result.add(obj); }
		}
	}
	void gatherAt(const CircF location, TemplateSet<TYPE*> & result) const {
		TYPE* obj;
		for (int i = 0; i < list.size(); ++i){
			obj = list[i];
			if (obj->intersectsCircle(location.center, location.radius)) { result.add(obj); }
		}
	}
	void gatherAt(const V2f position, TemplateSet<TYPE*> & result) const {
		TYPE* obj;
		for (int i = 0; i < list.size(); ++i){
			obj = list[i];
			if (obj->contains(position)) { result.add(obj); }
		}
	}
	void gatherCollisions(TYPE * subject, TemplateSet<TYPE*> & out_list) const {
		for (int i = 0; i < list.size(); ++i) {
			if (subject != list[i] && subject->intersects(list[i])) {
				out_list.add(list[i]);
			}
		}
	}
};

template<typename TYPE>
class CellSpacePartition : public Cell_I<TYPE> {
public:
	TemplateVectorList<TYPE*> * getList() { return NULL; }
	TemplateArray<Cell_I<TYPE>*> cells;
	V2i gridSize;
	V2f cellDimensions;
	CellSpacePartition(RectF totalArea, V2i gridSize) :Cell_I(totalArea) {
		this->gridSize.set(gridSize);
		cellDimensions.set(getWidth() / gridSize.x, getHeight() / gridSize.y);
		V2f cursor = this->min;
		RectF cursorR;
		int index = 0;
		cells.allocateToSize(gridSize.x*gridSize.y);
		for (int r = 0; r < gridSize.y; ++r)
		{
			for (int c = 0; c < gridSize.x; ++c)
			{
				cursorR.set(cursor, cursor + cellDimensions);
				cells[index] = new PartitionCell<TYPE>(cursorR);
				index++;
				cursor.x += cellDimensions.x;
			}
			cursor.x = this->min.x;
			cursor.y += cellDimensions.y;
		}
	}
	~CellSpacePartition(){ for (int i = 0; i < cells.size(); ++i) delete cells[i]; }
private:
	int getIndex(int row, int col) const { return row*gridSize.x + col; }
public:
	void gatherCellListFor(TYPE * s, TemplateVector<int> & list) const {
		RectF location(s->getCenter(), s->getRadius());
		RectF rect = getGridIndexeRange(location);
		int index;
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				index = getIndex(r, c);
				if (index < 0 || index >= cells.size()) {
					int i = 0; i = 1 / i;
					continue;
				}
				if (s->intersects(cells[index])) {
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

	void gatherCollisions(TYPE * subject, TemplateSet<TYPE*> & out_list) const {
//		printf("\n%s %f,%f:%f", subject->getTypeName(), subject->getCenter().x, subject->getCenter().y, subject->getRadius());
		RectF location(subject->getCenter(), subject->getRadius());
		RectF rect = getGridIndexeRange(location);
		int index;
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				index = getIndex(r, c);
//				printf("%d.%d:%d ",r,c,index);
				if (subject->intersectsAABB(cells[index]->min, cells[index]->max)) {
					cells[index]->gatherCollisions(subject, out_list);
				}
			}
		}
	}

	// used to test the minimal-cell-testing raycasting algorithm
	void raycastCellList(V2f start, V2f direction, float maxDistance, TemplateSet<int> & cellIDs) const {
		V2f out_point, out_normal;
		float out_dist;
		V2f point = start;
		if (!contains(point)) {
			if (!raycast(start, direction, out_dist, point, out_normal)){
				return;
			}
		}
		V2f rowcol = getGridIndex(point);
		// when this line steps on the grid, how does it step?
		V2i step(
			((direction.x > 0) ? 1 : ((direction.x < 0) ? -1 : 0)),
			((direction.y > 0) ? 1 : ((direction.y < 0) ? -1 : 0))
			);
		const int NUM_DIMENSIONS = 2;
		int orderOfMoves[NUM_DIMENSIONS];
		// weight the dimensions stepped
		if (abs(direction.x) > abs(direction.y)) {
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
			if (!isOutOfBounds && (cells[index]->contains(start) || (cells[index]->raycast(start, direction, out_dist, out_point, out_normal) && out_dist < maxDistance))) {
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

	bool raycastContainer(V2f start, V2f direction, float maxDistance, bool dontCareAboutObstacle, TYPE * & out_obstacle,
		float & out_dist, V2f & out_point, V2f & out_normal, TYPE ** ignoreList, int ignoreCount) const {
		V2f point = start;
		if (!contains(point)) {
			if (!raycast(start, direction, out_dist, point, out_normal)){
				return false;
			}
		}
		V2f rowcol = getGridIndex(point);
		// when this line steps on the grid, how does it step?
		V2i step(
			((direction.x > 0) ? 1 : ((direction.x < 0) ? -1 : 0)),
			((direction.y > 0) ? 1 : ((direction.y < 0) ? -1 : 0))
			);
		const int NUM_DIMENSIONS = 2;
		int orderOfMoves[NUM_DIMENSIONS];
		// weight the dimensions stepped based on the direction being travelled
		if (abs(direction.x) > abs(direction.y)) {
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
			if (!isOutOfBounds && (cells[index]->contains(start) || (cells[index]->raycast(start, direction, out_dist, out_point, out_normal) && out_dist < maxDistance))) {
				hitInCell = cells[index]->raycastContainer(start, direction, maxDistance, dontCareAboutObstacle, out_obstacle,
					out_dist, out_point, out_normal, ignoreList, ignoreCount);
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
		rowcol -= this->min;
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
		RectF rowcol(location);
		V2f totalSize = max - min;
		rowcol.min -= this->min;
		rowcol.max -= this->min;
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

	void gatherAt(const RectF location, TemplateSet<TYPE*> & result) const {
		VALID_cells_index(location, {
			if (cells[index]->intersectsAABB(location.min, location.max)) {
				cells[index]->gatherAt(location, result);
			}
		})
	}
	void gatherAt(const V2f position, TemplateSet<TYPE*> & result) const {
		//return gatherAt(RectF(position, 0), result);
		VALID_cells_index(RectF(position, 0), {
			if (cells[index]->contains(position)) {
				cells[index]->gatherAt(position, result);
			}
		})
	}
	void gatherAt(const CircF circle, TemplateSet<TYPE*> & result) const {
		//return gatherAt(RectF(position, 0), result);
		VALID_cells_index(RectF(circle.center, circle.radius), {
			if (cells[index]->intersectsCircle(circle.center, circle.radius)) {
				cells[index]->gatherAt(circle, result);
			}
		})
	}
#undef VALID_cells_index

	void add(TYPE* entity) {
		RectF location(entity->getCenter(), entity->getRadius());
		RectF rect = getGridIndexeRange(location);
		int index;
		for (int r = (int)rect.min.y; r <= rect.max.y; ++r) {
			for (int c = (int)rect.min.x; c <= rect.max.x; ++c) {
				index = getIndex(r, c);
				if (index < 0 || index >= cells.size()) {
					printf("B@D TH!NG T#4T H4PP3ND ");
					continue;
				}
				if (entity->intersects(cells[index])) {
					cells[index]->add(entity);
				}
			}
		}
	}
	void clear() {
		for (int i = 0; i < cells.size(); ++i)
			cells[i]->clear();
	}
	bool has(RectF area) const {
		return this->intersectsAABB(area.min, area.max);
	}

	void glDraw() {
		for (int i = 0; i < cells.size(); ++i) {
			cells[i]->glDraw(false);
		}
	}
};
