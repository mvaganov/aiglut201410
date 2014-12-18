#pragma once
#include <GL/freeglut.h>
#include "v2.h"
#include "templatevector.h"
#include "templateset.h"
#include "box.h"
#include "cone.h"

/** convex polygons only please! TODO implement a structure to deal with concave polygons, by splitting convex up into multiple concave */
template <typename TYPE>
class Polygon2 {
	/** if there was a circle that this polygon were instribed in, what would it's radius be */
	TYPE radius;
	/** where is the center of the circle that this polygon is inscribed in */
	V2<TYPE> center;
public:

	/** offsets relative to the origin, expected in CW order */
	TemplateVector< V2<TYPE> > points;
	struct pair {
		int a, b;
		pair(int a, int b) :a(a), b(b){}
		pair() :a(-1), b(-1){}
		bool operator==(pair const & pair) const { return (a == pair.a && b == pair.b) || (a == pair.b && b == pair.a); }
		int min() const { return (a < b) ? a : b; }
		int max() const { return (a > b) ? a : b; }
		bool operator<(pair const & pair) const {
			int m = min(), pairm = pair.min();
			if (m == pairm) { return max() < pair.max(); }
			return m < pairm;
		}
		bool operator>(pair const & pair) const {
			int m = min(), pairm = pair.min();
			if (m == pairm) { return max() > pair.max(); }
			return m > pairm;
		}
		bool contains(const int value) const { return a == value || b == value; }
	};
	/** how the points connect. If this is size 0 then the points are considered a line loop */
	TemplateSet<pair> connectionPairs;
	/** where the center of this polygon is */
	V2<TYPE> origin;
	TYPE rotation;

	Polygon2() : radius(0) {}
	Polygon2(V2<TYPE> const & origin, V2<TYPE> const & rotation) : origin(origin), rotation(rotation) : radius(0) {}

	bool isValid() const { return radius > 0; }
	bool isLineStripPolygon() const { return connectionPairs.size() == 0; }

	int getLineCount() const { return isLineStripPolygon() ? points.size() : connectionPairs.size(); }
	/** @return the line in relation to the origin (before translation to the origin and rotation) */
	void gatherLineRelative(const int index, V2<TYPE> & out_a, V2<TYPE> & out_b) const {
		if (isLineStripPolygon()) {
			out_a = points[index];
			out_b = points[(index + 1) % points.size()];
		} else {
			out_a = points[connectionPairs[index].a];
			out_b = points[connectionPairs[index].b];
		}
	}
	
	void gatherLine(const int index, V2<TYPE> & out_a, V2<TYPE> & out_b) const {
		if (isLineStripPolygon()) {
			out_a = points[index];
			out_b = points[(index + 1) % points.size()];
		} else {
			out_a = points[connectionPairs[index].a];
			out_b = points[connectionPairs[index].b];
		}
		out_a.rotate(rotation);
		out_b.rotate(rotation);
		out_a += origin;
		out_b += origin;
	}

	void setPoints(const V2<TYPE> * const & list, const int count) {
		points.setSize(count);
		for (int i = 0; i < points.size(); ++i) {
			points[i] = list[i];
			center += points[i];
		}
		center /= (TYPE)points.size();
		for (int i = 0; i < points.size(); ++i) {
			TYPE mag = (center - points[i]).magnitude();
			if (mag > radius) {
				radius = mag;
			}
		}
	}

	Polygon2(V2<TYPE> const & origin, TYPE const & rotation, V2<TYPE> * const & list, const int count)
		:radius(0), origin(origin), rotation(rotation) {
		setPoints(list, count);
		// TODO implement the inscribed-circle center
		if (V2<TYPE>::isPolyCW(points.getRawList(), points.size())) {
			points.reverse();
		}
		if (!V2<TYPE>::isPolyCCW(points.getRawList(), points.size())) {
			printf("polygon not ordered correctly!\n");
			radius = -1;
			return;
		}
	}

	Polygon2(Box<TYPE> box):radius(0),origin(box.center),rotation(box.rotation.piRadians()) {
		V2<TYPE> boxP[4];
		box.writePointsRelative(boxP, 4);
		setPoints(boxP, 4);
	}
	Polygon2(V2<TYPE> const & origin, TYPE const & rotation, const V2<TYPE> * const & list, const int count, const pair * const & pairs, const int pairsCount) {
		set(origin, rotation, list, count, pairs, pairsCount, true);
	}

	void set(V2<TYPE> const & origin, TYPE const & rotation, const V2<TYPE> * const & list, const int count, const pair * const & pairs, const int pairsCount, bool attemptClockwiseOptimization) {
		this->origin = origin;
		this->rotation = rotation;
		bool clockwiseOptimize = attemptClockwiseOptimization;
		if (clockwiseOptimize) {
			// try to make this non-clockwise polygon into a clockwise polygon.
			//for (int i = 0; i < pairsCount; ++i) { printf("%d %d, ", pairs[i].a, pairs[i].b); } printf("\n");
			// every point needs to be referenced by the pairs twice to make the loop (to and from)
			int numIncompletes = 0, start = -1, end = -1;
			TemplateVector<int> ignoredPoints;
			for (int i = 0; i < count; ++i) {
				int timesReferenced = 0;
				for (int p = 0; p < pairsCount; ++p) {
					if (pairs[p].contains(i)) {
						timesReferenced++;
						if (timesReferenced > 2) break;
					}
				}
				//printf("%d:%d ", i, timesReferenced);
				if(timesReferenced == 0) {
					ignoredPoints.add(i);
				} else if (timesReferenced == 1) {
					if (start == -1) start = i;
					else end = i;
					numIncompletes++;
					if(numIncompletes > 2) {
						clockwiseOptimize = false;
						break;
					}
				} else if (timesReferenced > 2) {
					clockwiseOptimize = false;
					break;
				}
			}
			if(clockwiseOptimize && ignoredPoints.size() > 0) {
				//printf("removing unused points\n");
				for(int i = ignoredPoints.size()-1; i >= 0; --i) {
					points.remove(ignoredPoints[i]);
					if (ignoredPoints[i] < start) start--;
					if (ignoredPoints[i] < end) end--;
				}
			}
			//printf("\n");
			// every point needs to be clockwise to the next one
			if (clockwiseOptimize) {
				// if it's true, the points can be re-arranged to be clockwise, and the pairs can be discarded
				connectionPairs.clear();
				V2f startP, endP;
				// TODO if incompletes is 2 and start and end are valid, make connection pairs as a strip from start to end (or end to start, whichever is actually CW)
				if (numIncompletes == 2 && start != -1 && end != -1) {
					startP = list[start];
					endP = list[end];
				}
				//printf("%s\n", (numIncompletes==0)?"clockwise!":"close enough to clockwise");
				setPoints(list, count);
				calculatePolygonCW(points.getRawList(), points.size(), center);
				if (numIncompletes == 2 && start != -1 && end != -1) {
					//for (int i = 0; i < points.size(); ++i) {
					//	if (points[i] == startP) printf("start@%d ", i);
					//	if (points[i] == endP) printf("end@%d ", i);
					//}
					int index;
					for (index = 0; index < points.size(); ++index) {
						if (points[index] == endP || points[index] == startP) {
							if (points[index + 1] == endP || points[index + 1] == startP) {
								index++;
							}
							break;
						}
					}
					for (int i = 0; i < points.size() - 1; ++i) {
						int a = (index + i) % points.size(), b = (index + i + 1) % points.size();
					//	printf("%d %d, ", a, b);
						connectionPairs.add(pair(a, b));
					}
					//printf("\n");
				}
			}
		}
		if (!clockwiseOptimize) {
			setPoints(list, count);
			connectionPairs.clear();
			connectionPairs.ensureCapacity(pairsCount);
			for (int i = 0; i < pairsCount; ++i) {
				if (pairs[i].a < points.size() && pairs[i].b < points.size()) {
					connectionPairs.add(pairs[i]);
				}
				else {
					int i = 0; i = 1 / i; // OOB...
				}
				connecitonPairsAreSafe();
			}
			connecitonPairsAreSafe();
		}
		connecitonPairsAreSafe();
	}

	bool connecitonPairsAreSafe() const {
		if (this->isLineStripPolygon()) return true;
		for (int i = 0; i < connectionPairs.size(); ++i) {
			if (connectionPairs[i].a >= points.size() || connectionPairs[i].b >= points.size()) {
				int x = 0; x = 1 / x;
				return false;
			}
		}
		return true;
	}

	V2<TYPE> getCenter() const { return origin + center; }
	float getRadius() const { return radius; }

	bool contains(V2<TYPE> const & p) const {
		V2<TYPE> relativeP = p - origin, a, b, lineP, lineC, lineDelta;
		relativeP.rotate(-rotation);
		bool knownToBeClockwise = this->isLineStripPolygon();
		bool isClockwise = knownToBeClockwise;
		TYPE pSign, cSign;
		for (int i = 0; i < getLineCount(); ++i) {
			gatherLineRelative(i, a, b);
			lineDelta = b - a;
			lineP = relativeP - a;
			if (!knownToBeClockwise) {
				lineC = center - a;
				cSign = lineC.sign(lineDelta);
				isClockwise = (cSign > 0);
			}
			// if this point is not on the same side of this line as the center
			pSign = lineP.sign(lineDelta);
			if ((pSign > 0) != isClockwise) {
				return false;
			}
		}
		return true;
	}

	bool raycast(V2<TYPE> const & rayStart, V2<TYPE> const & rayDirection,
		TYPE & out_dist, V2<TYPE> & out_point, V2<TYPE> & out_normal) const {
		V2<TYPE> relativeRayStart = rayStart - origin;
		relativeRayStart.rotate(-rotation);
		V2<TYPE> relativeRayDirection = rayDirection.rotated(-rotation);
		bool hitSomething = false;
		// test if the ray would hit a sphere that contains this polygon
		if (V2<TYPE>::rayCrossesCircle(relativeRayStart, relativeRayDirection, center, radius)) {
			V2<TYPE> a, b, delta, point, centerDelta;
			TYPE line, ray, dist, bestDist = -1;
			for (int i = 0; i < this->getLineCount(); ++i) {
				gatherLineRelative(i, a, b);
				bool hitPossible = V2<TYPE>::rayIntersection(a, b, relativeRayStart, relativeRayStart + relativeRayDirection, line, ray);
				if (hitPossible && ray >= 0 && line >= 0 && line < 1) {
					delta = b - a;
					point = a + (delta * line);
					dist = (relativeRayStart - point).magnitude();
					if (bestDist < 0 || dist < bestDist) {
						bestDist = dist;
						out_point = point;
						out_normal = delta.perp().normal();
						out_dist = (relativeRayDirection * ray).magnitude();
						// make the normal face away from the center
						centerDelta = out_point - center;
						if (V2<TYPE>::dot(out_normal, centerDelta) < 0) {
							out_normal *= -1;
						}
						hitSomething = true;
					}
				}
			}
		}
		if (hitSomething) {
			out_point.rotate(rotation);
			out_point += origin;
			out_normal.rotate(rotation);
		}
		return hitSomething;
	}

	V2<TYPE> getClosestPointOnEdge(const V2<TYPE> point, V2<TYPE> & out_normal) const {
		TYPE closestDist = -1, d;
		V2<TYPE> relativeP = point - origin, a, b;
		relativeP.rotate(-rotation);
		V2<TYPE> closestPoint, p, delta;
		int closestLineIndex = -1;
		bool isOnLine;
		// get the closest point in the list
		closestPoint = points.get(V2<TYPE>::indexOfClosest(relativeP, points.getRawListConst(), points.size()));
		delta = relativeP - closestPoint;
		out_normal = delta.normal();
		closestDist = delta.magnitude();
		// check the closest point on each line
		int linesToCheck = getLineCount();
		for (int i = 0; i < linesToCheck; ++i) {
			gatherLineRelative(i, a, b);
			isOnLine = V2<TYPE>::closestPointOnLine(a, b, relativeP, p);
			if (isOnLine) {
				delta = p - relativeP;
				d = delta.magnitude();
				// keep the closest solution
				if (d < closestDist) {
					closestDist = d;
					closestLineIndex = i;
					closestPoint = p;
				}
			}
		}
		// if there are no solutions, give the closest point in the list
		if (closestLineIndex >= 0) {
			delta = points.get((closestLineIndex + 1) % points.size()) - points.get(closestLineIndex);
			out_normal = delta.normal().perp();
		}
		closestPoint.rotate(rotation);
		out_normal.rotate(rotation);
		return closestPoint + origin;
	}

	void glDraw(bool filled) const {
		glPushMatrix();
		origin.glTranslate();
		TYPE degrees = rotation * 180 / (TYPE)V_PI;
		glRotatef(degrees, 0, 0, 1);
		if (isLineStripPolygon()) {
			//glDrawCircle(center, radius, false);
			if (!filled) {
				glBegin(GL_LINE_LOOP);
			} else {
				glBegin(GL_TRIANGLE_FAN);
				center.glVertex();
			}
			for (int i = 0; i < points.size(); ++i) {
				points[i].glVertex();
			}
			if (points.size() > 0) points[0].glVertex();
			glEnd();

			//if (filled) {
			//	glColor3f(0, 0, 0);
			//	glBegin(GL_LINES);
			//	for (int i = 0; i < points.size(); ++i) {
			//		center.glVertex();
			//		points[i].glVertex();
			//	}
			//	glEnd();
			//}
		} else {
			connecitonPairsAreSafe();
			if (!filled) {
				glBegin(GL_LINES);
			} else {
				glBegin(GL_TRIANGLES);
			}
			int numa, numb;
			for (int i = 0; i < connectionPairs.size(); ++i) {
				numa = connectionPairs[i].a;
				numb = connectionPairs[i].b;
				if (filled) center.glVertex();
				points[numa].glVertex();
				points[numb].glVertex();
			}
			glEnd();
		}
		glPopMatrix();
		//glBegin(GL_LINES);
		//V2<TYPE> a, b;
		//for (int i = 0; i < getLineCount(); ++i) {
		//	gatherLine(i, a, b);
		//	a.glDrawTo(b);
		//}
		////printf("\n");
		//glEnd();
	}

	bool intersectsCircle(const V2<TYPE> & circCenter, const TYPE circRadius) const {
		V2<TYPE> relativeP = circCenter - origin;
		relativeP.rotate(-rotation);
		V2<TYPE> delta = center - relativeP, a, b;
		TYPE dist = delta.magnitude();
		// check if the circle approximation would collide.
		if (dist < radius + circRadius) {
			// check if the circle's center is inside the approximation, and if it is inside the polygon
			if (dist < radius && relativeP.isInsidePolyCW(points.getRawListConst(), points.size()))
				return true;
			for (int i = 0; i < points.size(); ++i) {
				delta = points.get(i) - relativeP;
				dist = delta.magnitude();
				// check if any of the points are inside
				if (dist < circRadius) return true;
			}
			for (int i = 0; i < getLineCount(); ++i) {
				gatherLineRelative(i, a, b);
				// check if a line crosses the circle
				if (i > 0 && V2<TYPE>::lineCrossesCircle(a, b, relativeP, circRadius, delta)) {
					return true;
				}
			}
		}
		return false;
	}

	bool intersectsPolygon(Polygon2<TYPE> const & poly) const {
		V2<TYPE> delta = (origin + center.rotated(rotation)) - (poly.origin + poly.center.rotated(poly.rotation));
		TYPE dist = delta.magnitude();
		// check if circle approximations collide
		if (dist < radius + poly.radius) {
			V2<TYPE> a, b, c, d;
			// check if any of his points are in my space, or mine are in his space (only need to check one)
			gatherLine(0, a, b);
			poly.gatherLine(0, c, d);
			if (contains(c) || poly.contains(a)) return true;
			// check if any lines cross
			for (int mine = 0; mine < getLineCount(); ++mine) {
				gatherLine(mine, a, b);
				for (int his = 0; his < poly.getLineCount(); ++his) {
					poly.gatherLine(his, c, d);
					if (V2<TYPE>::lineIntersection(a, b, c, d)) return true;
				}
			}
		}
		return false;
	}

	bool intersectsBox(Box<TYPE> box) const {
		// treat the box like another polygon
		Polygon2<TYPE> boxPoly(box);// .center, box.rotation.piRadians(), boxP, 4);
		return intersectsPolygon(boxPoly);
	}

	bool intersectsAABB(AABB<TYPE> const & box) const {
		return intersectsBox(Box<TYPE>(box));
	}

	bool intersectsAABB(V2<TYPE> const & min, V2<TYPE> const & max) const {
		return intersectsAABB(AABB<TYPE>(min,max));
	}

	bool intersectsCone(Cone<TYPE> const & cone) const {
		V2<TYPE> delta = (origin + center.rotated(rotation)) - cone.origin;
		TYPE dist = delta.magnitude();
		// check if the circle approximations meet
		if (dist < cone.radius + radius) {
			if (contains(cone.origin)) return true;
			V2<TYPE> a, b;
			for (int i = 0; i < getLineCount(); ++i) {
				gatherLine(i, a, b);
				if (cone.contains(a) || cone.contains(b) || cone.intersectsLine(a, b)) return true;
			}
		}
		return false;
	}

	/**
	* @param io_points __OUT __IN a list of points to re-arrange in order (clockwise)
	* @param numPoints how many points there are
	* @param out_center where the average location of these points is
	*/
	static void calculatePolygonCW(V2<TYPE> * io_points, const int numPoints, V2<TYPE>  & out_center) {
		if (numPoints == 0) return;
		// calculate where to start the polygon, and where to center it too.
		int start = 0;
		out_center = io_points[0];
		for (int i = 1; i < numPoints; ++i) {
			if (io_points[i].y > io_points[start].y) { start = i; }
			out_center += io_points[i];
		}
		out_center /= (TYPE)numPoints;
		// calculate the angle of the starting point
		V2<TYPE>  currentAngle, testAngle, *startAngle = &io_points[start];
		currentAngle = (io_points[start] - out_center).normal();
		// calculate the angles of all of the points relative to the starting point
		TYPE * angles = new TYPE[numPoints];
		angles[start] = 0;
		for (int i = 0; i < numPoints; ++i) {
			if (i != start) {
				testAngle = (io_points[i] - out_center).normal();
				angles[i] = currentAngle.piRadians(testAngle);
				if (angles[i] < 0) angles[i] += (TYPE)V_2PI; // no going backwards
			}
		}
		// sort the list of points according to it's angle
		quicksort<TYPE, V2f>(angles, io_points, 0, numPoints - 1);
		// the angle list has served it's purpose.
		delete[] angles;
	}

	template<typename COMPARE, typename PARALLEL>
	static void quicksort(COMPARE * comparorList, PARALLEL * parallelList, int first, int last) {
		int i = first - 1, j = last;
		COMPARE v = comparorList[last], tempf;
		if (last <= first) return;
		PARALLEL tempv;
		do {
			while (comparorList[++i] < v);
			while (v < comparorList[--j]) if (j == first) break;
			if (i >= j) break;
			//swap(i, j);
			tempf = comparorList[i]; comparorList[i] = comparorList[j]; comparorList[j] = tempf;
			tempv = parallelList[i]; parallelList[i] = parallelList[j]; parallelList[j] = tempv;
		} while (true);
		//swap(i, last);
		tempf = comparorList[i]; comparorList[i] = comparorList[last]; comparorList[last] = tempf;
		tempv = parallelList[i]; parallelList[i] = parallelList[last]; parallelList[last] = tempv;
		quicksort(comparorList, parallelList, first, i - 1);
		quicksort(comparorList, parallelList, i + 1, last);
	}


};

typedef Polygon2<float> Polygon2f;
