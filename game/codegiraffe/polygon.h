#pragma once
#include <GL/freeglut.h>
#include "v2.h"
#include "templatevector.h"
#include "box.h"
#include "cone.h"

/** convex polygons only please! TODO implement a structure to deal with concave polygons, by splitting convex up into multiple concave */
template <typename TYPE>
class Polygon2 {
	/** if there was a circle that this polygon were instribed in, what would it's radius be */
	TYPE radius;
	/** where is the center of the circle that this polygon is inscribed in */
	V2<TYPE> center;

	int hitSide;
	V2<TYPE> outpoint;
public:
	/** offsets relative to the origin, expected in CW order */
	TemplateVector< V2<TYPE> > points;
	V2<TYPE> origin;
	TYPE rotation;

	Polygon2() : radius(0), hitside(-1) {}
	Polygon2(V2<TYPE> const & origin, V2<TYPE> const & rotation) : origin(origin), rotation(rotation) : radius(0), hitside(-1) {}

	bool isValid() { return radius > 0; }

	Polygon2(V2<TYPE> const & origin, TYPE const & rotation, V2<TYPE> * const & list, const int count)
		:radius(0), origin(origin), rotation(rotation) {
		points.add(count, list);
		// TODO implement the inscribed-circle center
		for (int i = 0; i < points.size(); ++i) {
			center += points[i];
		}
		center /= (float)points.size();
		if (V2f::isPolyCW(points.getRawList(), points.size())) {
			points.reverse();
		}
		if (!V2f::isPolyCCW(points.getRawList(), points.size())) {
			printf("polygon not ordered correctly!\n");
			radius = -1;
			return;
		}
		for (int i = 0; i < points.size(); ++i) {
			float mag = (center - points[i]).magnitude();
			if (mag > radius) {
				radius = mag;
			}
		}
	}

	V2<TYPE> getCenter() const { return origin + center; }

	bool contains(V2<TYPE> const & p) const {
		V2<TYPE> relativeP = p - origin;
		relativeP.rotate(-rotation);
		return contains(relativeP, points.getRawListConst(), points.size());
	}

	static bool contains(V2<TYPE> const & p, const V2<TYPE> * const & points, const int numPoints) {
		for (int i = 0; i < numPoints; ++i) {
			V2<TYPE> a = points[i];
			V2<TYPE> b = points[(i + 1) % numPoints];
			if (p.isCW(a, b))
				return false;
		}
		return true;
	}

	bool raycast(V2<TYPE> const & rayStart, V2<TYPE> const & rayDirection,
		TYPE & out_dist, V2<TYPE> & out_point, V2<TYPE> & out_normal) const {
		V2<TYPE> relativeRayStart = rayStart - origin;
		relativeRayStart.rotate(-rotation);
		V2<TYPE> relativeRayDirection = rayDirection.rotated(-rotation);
		// test if the ray would hit a sphere that contains this polygon
		if (V2<TYPE>::rayCrossesCircle(relativeRayStart, relativeRayDirection, center, radius)) {
			// test the ray against every line
			int whichSegmentWasHit = V2<TYPE>::raycastPolygon(
				relativeRayStart, relativeRayDirection, points.getRawListConst(), points.size(),
				out_dist, out_point, out_normal);
			if (whichSegmentWasHit >= 0) {
				out_point.rotate(rotation);
				out_point += origin;
				out_normal.rotate(rotation);
			}
			return whichSegmentWasHit >= 0;
		}
		return false;
	}

	V2<TYPE> getClosestPointOnEdge(const V2<TYPE> point, V2<TYPE> & out_normal) const {
		TYPE closestDist = -1, d;
		V2<TYPE> relativeP = point - origin;
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
		for (int i = 0; i < points.size(); ++i) {
			isOnLine = V2<TYPE>::closestPointOnLine(points.get(i), points.get((i + 1) % points.size()), relativeP, p);
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
		glDrawCircle(center, radius, false);
		if (!filled) {
			glBegin(GL_LINE_LOOP);
		}
		else {
			glBegin(GL_TRIANGLE_FAN);
		}
		for (int i = 0; i < points.size(); ++i) {
			points.get(i).glVertex();
		}
		glEnd();
		if (hitSide >= 0) {
			glColor3f(1, 0, 0);
			glBegin(GL_LINES);
			points[hitSide].glVertex();
			points[(hitSide + 1) % points.size()].glVertex();
			glEnd();
			glDrawCircle(outpoint, .05f, true);
		}
		glPopMatrix();
	}

	bool intersectCircle(const V2<TYPE> & circCenter, const TYPE circRadius) const {
		V2<TYPE> relativeP = circCenter - origin;
		relativeP.rotate(-rotation);
		V2<TYPE> delta = center - relativeP;
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
				// if a pair of points are outside the circle, check if the line between them crosses the circle
				if (i > 0 && V2<TYPE>::lineCrossesCircle(points[i - 1], points[i], relativeP, circRadius, delta)) {
					return true;
				}
			}
			// check the last pair
			if (V2<TYPE>::lineCrossesCircle(points[points.size() - 1], points[0], relativeP, circRadius, delta)) {
				return true;
			}
		}
		return false;
	}

	bool intersectsPolygon(Polygon2<TYPE> const & poly) const {
		V2<TYPE> delta = (origin + center.rotated(rotation)) - (poly.origin + poly.center.rotated(poly.rotation));
		float dist = delta.magnitude();
		// check if circle approximations collide
		if (dist < radius + poly.radius) {
			// make the points relative to this polygon
			TemplateArray<V2<TYPE> > relativePoints(poly.points.size(), poly.points.getRawListConst());
			for (int i = 0; i < relativePoints.size(); ++i) {
				relativePoints[i].rotate(poly.rotation);
				relativePoints[i] += poly.origin;
				relativePoints[i] -= origin;
				relativePoints[i].rotate(-rotation);
			}
			// check intersection
			return intersectPolygons(points.getRawListConst(), points.size(),
				relativePoints.getRawListConst(), relativePoints.size());
		}
		return false;
	}

	static bool intersectPolygons(const V2<TYPE> * const & pointsA, const int numPointsA,
		const V2<TYPE> * const & pointsB, const int numPointsB) {
		// if a point from B is in A, or vice-versa, then there is clear intersection
		if (contains(pointsB[0], pointsA, numPointsA) || contains(pointsA[0], pointsB, numPointsB)) {
			return true;
		}
		// check the lines from B colliding with A
		for (int a = 0; a < numPointsA; ++a) {
			for (int b = 0; b < numPointsB; ++b) {
				if (V2<TYPE>::lineIntersection(
					pointsA[a], pointsA[(a + 1) % numPointsA],
					pointsB[b], pointsB[(b + 1) % numPointsB]))
				return true;
			}
		}
		return false;
	}

	bool intersectsBox(Box<TYPE> box) const {
		// treat the box like another polygon
		V2<TYPE> boxP[4];
		box.writePoints(boxP, 4);
		for (int i = 0; i < 4; ++i) {
			boxP[i] -= origin;
			boxP[i].rotate(-rotation);
		}
		return intersectPolygons(points.getRawListConst(), points.size(), boxP, 4);
	}

	bool intersectsCone(Cone<TYPE> const & cone) const {
		V2<TYPE> delta = (origin + center.rotated(rotation)) - cone.origin;
		float dist = delta.magnitude();
		// check if the circle approximations meet
		if (dist < cone.radius + radius) {
			// make the points relative to absolute space
			TemplateArray<V2<TYPE> > relativePoints(points.size(), points.getRawListConst());
			for (int i = 0; i < relativePoints.size(); ++i) {
				relativePoints[i].rotate(rotation);
				relativePoints[i] += origin;
			}
			// use the cone's intersectsPolygonLineList method
			return cone.intersectsPolygonLineList(relativePoints.getRawListConst(), relativePoints.size());
		}
		return false;
	}
};

typedef Polygon2<float> Polygon2f;
