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
		center /= (float)points.size();
		for (int i = 0; i < points.size(); ++i) {
			float mag = (center - points[i]).magnitude();
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

	Polygon2(V2<TYPE> const & origin, TYPE const & rotation, const V2<TYPE> * const & list, const int count, const pair * const & pairs, const int pairsCount)
		:radius(0), origin(origin), rotation(rotation) {
		setPoints(list, count);
		for (int i = 0; i < pairsCount; ++i) {
			connectionPairs.add(pairs[i]);
		}
		//for (int i = 0; i < connectionPairs.size(); ++i) {
		//	printf("%d %d\n", connectionPairs[i].min(), connectionPairs[i].max());
		//}
	}

	V2<TYPE> getCenter() const { return origin + center; }

	bool contains(V2<TYPE> const & p) const {
		V2<TYPE> relativeP = p - origin, a, b, lineP, lineC, lineDelta;
		relativeP.rotate(-rotation);
		bool knownToBeClockwise = this->isLineStripPolygon();
		bool isClockwise = knownToBeClockwise;
		float pSign, cSign;
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
			float line, ray, dist, bestDist = -1;
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
		float degrees = rotation * 180 / (float)V_PI;
		glRotatef(degrees, 0, 0, 1);
		if (isLineStripPolygon()) {
			glDrawCircle(center, radius, false);
			if (!filled) {
				glBegin(GL_LINE_LOOP);
			}
			else {
				glBegin(GL_TRIANGLE_FAN);
				center.glVertex();
			}
			for (int i = 0; i < points.size(); ++i) {
				points[i].glVertex();
			}
			glEnd();
		} else {
			glBegin(GL_LINES);
			for (int i = 0; i < connectionPairs.size(); ++i) {
				points[connectionPairs[i].a].glVertex();
				points[connectionPairs[i].b].glVertex();
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

	bool intersectCircle(const V2<TYPE> & circCenter, const TYPE circRadius) const {
		V2<TYPE> relativeP = circCenter - origin;
		relativeP.rotate(-rotation);
		V2<TYPE> delta = center - relativeP, a, b;
		float dist = delta.magnitude();
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
		float dist = delta.magnitude();
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

	bool intersectsCone(Cone<TYPE> const & cone) const {
		V2<TYPE> delta = (origin + center.rotated(rotation)) - cone.origin;
		float dist = delta.magnitude();
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
};

typedef Polygon2<float> Polygon2f;
