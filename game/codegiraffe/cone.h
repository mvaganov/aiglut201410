#pragma once
#include <GL/freeglut.h>
#include "codegiraffe/v2.h"
#include "codegiraffe/circle.h"

template <typename TYPE>
class Cone {
public:
	V2<TYPE> origin;
	TYPE radius;
	TYPE startAngle;
	TYPE endAngle;

	Cone(V2<TYPE> origin, TYPE rad, TYPE startAngle, TYPE endAngle)
		:origin(origin), radius(rad), startAngle(startAngle), endAngle(endAngle)
	{}

	V2<TYPE> getCenter() const { 
		TYPE totalAngle = endAngle - startAngle;
		TYPE middleAngle = totalAngle / 2 + startAngle;
		V2<TYPE> midRad = V2<TYPE>(middleAngle);
		TYPE rad = radius / 2;
		if (totalAngle > V_PI) {
			rad -= (float)((totalAngle - V_PI) / V_PI * (radius / 2));
		}
		return midRad * rad;
	}

	static bool isInAngle(V2<TYPE> testPoint, V2<TYPE> origin, TYPE startAngle, TYPE endAngle) {
		V2<TYPE> relativeP = p - origin;
		V2<TYPE> startRad = V2<TYPE>(startAngle);
		V2<TYPE> endRad = V2<TYPE>(endAngle);
		return isInAngle(relativeP, startRad, endRad);
	}

	/**
	 * @param relativeP point relative to the origin
	 * @param startRad line from the origin that is CCW from the endRad
	 * @param endRad line from the origin that is CW from the startRad
	 */
	static bool isInAngle(V2<TYPE> relativeP, V2<TYPE> startRad, V2<TYPE> endRad) {
		bool cwStart = relativeP.sign(startRad) <= 0, cwEnd = relativeP.sign(endRad) <= 0;
		return (cwStart && !cwEnd);
	}

	bool contains(V2<TYPE> const & p) const { 
		V2<TYPE> relativeP = p - origin;
		V2<TYPE> startRad = V2<TYPE>(startAngle);
		V2<TYPE> endRad = V2<TYPE>(endAngle);
		if (isInAngle(relativeP, startRad, endRad)) {
			TYPE distance = relativeP.magnitude();
			return distance < radius;
		}
		return false;
	}

	int indexOfContained(const V2<TYPE> * points, const int numPoints, const int startIndex) const {
		V2<TYPE> startRad = V2<TYPE>(startAngle);
		V2<TYPE> endRad = V2<TYPE>(endAngle);
		for (int i = startIndex; i < numPoints; ++i) {
			V2<TYPE> relativeP = points[i] - origin;
			if (isInAngle(relativeP, startRad, endRad)) {
				TYPE distance = relativeP.magnitude();
				if (distance < radius)
					return i;
			}
		}
		return -1;
	}

	static int indexOfSmallestPositive(const TYPE * list, const int listSize) {
		int closest = -1;
		for (int i = 0; i < listSize; ++i) {
			if (list[i] >= 0 && ((closest == -1) || (list[i] < list[closest])))
				closest = i;
		}
		return closest;
	}

	bool raycast(V2<TYPE> const & rayStart, V2<TYPE> const & rayDirection,
		TYPE & out_dist, V2<TYPE> & out_point, V2<TYPE> & out_normal) const {
		CircF circ(origin, radius);
		V2<TYPE> arcPoint, arcNormal;
		TYPE distances[3] = { -1, -1, -1 };
		const int numDistances = sizeof(distances) / sizeof(distances[0]);
		bool hit = circ.raycast(rayStart, rayDirection, distances[0], arcPoint, arcNormal);
		V2<TYPE> relativeP = arcPoint - origin;
		V2<TYPE> startRad = V2<TYPE>(startAngle) * radius;
		V2<TYPE> endRad = V2<TYPE>(endAngle) * radius;
		if(!hit || !isInAngle(relativeP, startRad, endRad)){
			distances[0] = -1;
		}
		TYPE ray, line;
		hit = V2<TYPE>::rayIntersection(rayStart, rayStart+rayDirection, origin, origin + startRad, ray, line);
		V2<TYPE> startRad_point, startRad_normal;
		if(hit && ray > 0 && line > 0 && line < 1) {
			startRad_point = startRad * line + origin;
			startRad_normal = -V2<TYPE>(startAngle).perp();
			distances[1] = (rayDirection * ray).magnitude();
		} else distances[1] = -1;
		hit = V2<TYPE>::rayIntersection(rayStart, rayStart+rayDirection, origin, origin + endRad, ray, line);
		V2<TYPE> endRad_point, endRad_normal;
		if (hit && ray > 0 && line > 0 && line < 1) {
			endRad_point = endRad * line + origin;
			endRad_normal = V2<TYPE>(endAngle).perp();
			distances[2] = (rayDirection * ray).magnitude();
		} else distances[2] = -1;
		int closest = indexOfSmallestPositive(distances, numDistances);
		if (closest != -1) {
			out_dist = distances[closest];
			switch (closest) {
			case 0:	out_point = arcPoint;	out_normal = arcNormal;	break;
			case 1:	out_point = startRad_point;	out_normal = startRad_normal;	break;
			case 2:	out_point = endRad_point;	out_normal = endRad_normal;	break;
			}
		}
		return closest != -1;
	}

	V2<TYPE> getClosestPointOnEdge(const V2<TYPE> point, V2<TYPE> & out_normal) const {
		CircF circ(origin, radius);
		V2<TYPE> arcPoint = circ.getClosestPointOnEdge(point, out_normal);
		V2<TYPE> relativeP = arcPoint - origin;
		V2<TYPE> startRad = V2<TYPE>(startAngle) * radius,endRad = V2<TYPE>(endAngle) * radius;
		V2<TYPE> startRadCorner = origin + startRad;
		V2<TYPE> endRadCorner = origin + endRad;
		V2<TYPE> onStart, onEnd, p;
		bool arcPointValid = isInAngle(relativeP, startRad, endRad);
		TYPE line, ray;
		bool startLineHit = V2<TYPE>::rayIntersection(origin, startRadCorner, point, point + startRad.perp(), line, ray);
		if (startLineHit && line <= 0 || line >= 1) {
			startLineHit = false;
		} else { onStart = origin + startRad * line; }
		bool endLineHit = V2<TYPE>::rayIntersection(origin, endRadCorner, point, point + endRad.perp(), line, ray);
		if (endLineHit && line <= 0 || line >= 1) {
			endLineHit = false;
		} else { onEnd = origin + endRad * line; }
		TYPE distances[] = {
			(point - origin).magnitude(), // origin corner
			(point - startRadCorner).magnitude(), // start corner
			(point - endRadCorner).magnitude(), // end corner
			arcPointValid ? (point - arcPoint).magnitude() : distances[0], // arc line
			startLineHit ? (point - onStart).magnitude() : distances[1], // start line
			endLineHit ? (point - onEnd).magnitude() : distances[2], // end line
		};
		const int numDistances = sizeof(distances) / sizeof(distances[0]);
		int smallest = indexOfSmallestPositive(distances, numDistances);
		switch (smallest) {
		case 0:	out_normal = (point - origin).normal();	p = origin;	break;
		case 1:	out_normal = (point - startRadCorner).normal();	p = startRadCorner;	break;
		case 2:	out_normal = (point - endRadCorner).normal();	p = endRadCorner;	break;
		case 3: p = arcPoint; break;
		case 4: out_normal = -V2<TYPE>(startAngle).perp();	p = onStart;	break;
		case 5: out_normal = V2<TYPE>(endAngle).perp();	p = onEnd;	break;
		}
		return p;
	}
	void glDraw(bool filled) const {
		//glDrawCircle(origin, radius, false);
		V2<TYPE> start(startAngle);
		V2<TYPE> points[32];
		const int numPoints = sizeof(points)/sizeof(points[0]);
		TYPE anglePerPoint = (endAngle-startAngle) / (numPoints-1);
		V2<TYPE>::arc(start, V2<TYPE>(anglePerPoint), points, numPoints);
		V2<TYPE> startRad = start * radius;
		V2<TYPE> endRad = V2<TYPE>(endAngle) * radius;
		V2<TYPE> zero = V2<TYPE>::ZERO();
		glPushMatrix();
		glTranslatef(origin.x, origin.y, 0);
		if(!filled) {
			zero.glDrawTo(startRad);
			zero.glDrawTo(endRad);
			glScalef(radius, radius, 1);
			glBegin(GL_LINE_STRIP);
		} else {
			glBegin(GL_TRIANGLE_FAN);
			zero.glVertex();
		}
		for(int i = 0; i < numPoints; ++i) {
			points[i].glVertex();
		}
		glEnd();
		glPopMatrix();
	}

	bool intersectCircle(const V2<TYPE> & circCenter, const TYPE circRadius) const {
		V2<TYPE> delta = circCenter - origin;
		// check the easy stuff: against the 3 corners
		float dist = delta.magnitude();
		if (dist < circRadius)return true;
		V2<TYPE> start(startAngle);
		V2<TYPE> end(endAngle);
		V2<TYPE> startRad = start * radius;
		V2<TYPE> endRad = end * radius;
		V2<TYPE> relativeP = circCenter - origin;
		// if it is in the circle range, check the harder stuff
		if (dist < circRadius + radius) {
			if (isInAngle(relativeP, startRad, endRad)) {
				return true;
			}
			V2<TYPE> prpStart = start.perp(), prpEnd = end.perp();
			V2<TYPE> points[] = {
				circCenter + prpStart * circRadius,
				circCenter - prpStart * circRadius,
				circCenter + prpEnd * circRadius,
				circCenter - prpEnd * circRadius,
			};
			const int numPoints = sizeof(points) / sizeof(points[0]);
			for (int i = 0; i < numPoints; ++i) {
				if (isInAngle(points[i] - origin, startRad, endRad))
					return true;
			}
		}
		dist = (relativeP - startRad).magnitude();
		if (dist < circRadius) return true;
		dist = (relativeP - endRad).magnitude();
		if (dist < circRadius) return true;
		return false;
	}

	bool intersectsBox(Box<TYPE> box) const {
		// check radius intersect, to see if collision is possible
		TYPE boxRad = box.size.magnitude() / 2;
		V2<TYPE> delta = box.center - origin;
		TYPE dist = delta.magnitude();
		if (dist < radius + boxRad) {
			if (box.contains(origin)) return true;
			// check if box corners are in the cone
			const int numPoints = 4;
			V2<TYPE> points[numPoints];
			box.writePoints(points, numPoints);
			if (indexOfContained(points, numPoints, 0) >= 0) {
				return true;
			}
			// check if cone lines cross box lines
			V2<TYPE> startP = V2f(startAngle) * radius + origin;
			if (box.contains(startP)) return true;
			V2<TYPE> endP = V2f(endAngle) * radius + origin;
			if (box.contains(endP)) return true;
			TYPE boxLine, coneLine;
			for (int i = 0; i < numPoints; ++i) {
				if (V2<TYPE>::rayIntersection(origin, startP, points[i], points[(i + 1) % numPoints], coneLine, boxLine)
					&& coneLine > 0 && coneLine < 1 && boxLine > 0 && boxLine < 1) {
					return true;
				}
			}
		}
		return false;
	}

	bool intersectsCone(Cone<TYPE> otherCone) const {
		// check radius intersect
		V2<TYPE> delta = otherCone.origin - origin;
		TYPE dist = delta.magnitude();
		if (dist < radius + otherCone.radius) {
			// check circle collision against each other, if there is at least one collision there
			// check corners inside each other
			// check crossing lines
		}
		return false;
	}
};

typedef Cone<float> ConeF;