#pragma once

#include <math.h>

// PI constants redefined to reduce dependency on weird _USE_MATH_DEFINES issues
const double V_PI = 3.14159265358979323846;
const double V_2PI = V_PI * 2;
const double V_HALFPI = V_PI / 2;
const double V_QUARTERPI = V_PI / 4;

/**
 * A two dimensional vector class, which includes methods to help with 2D math, including rotation, and some 2D geometry testing.
 * Intended to be used with OpenGL, and more specifically, freeglut http://freeglut.sourceforge.net/docs/api.php
 * Written to be useful, readable, and informative for students studying 2D Vector math.
 *
 * @author mvaganov@hotmail.com October 2014
 */
template<typename TYPE>
struct V2 {
	/** position in 2 Dimensional space */
	TYPE x, y;
	/** how many dimensions a V2<TYPE> keeps track of */
	static const int NUM_DIMENSIONS = 2;
	static const int DIMENSION_X = 0;
	static const int DIMENSION_Y = 1;

	/** access X */
	const TYPE getX() const {return x;}
	/** access Y */
	const TYPE getY() const {return y;}
	/** mutate X */
	void setX(const TYPE a_value) { x = a_value; }
	/** mutate Y */
	void setY(const TYPE a_value) { y = a_value; }
	/** mutate X */
	void addX(const TYPE a_value) { x += a_value; }
	/** mutate Y */
	void addY(const TYPE a_value) { y += a_value; }
	/** mutate in bulk */
	void add(const TYPE a_x, const TYPE a_y) { x += a_x; y += a_y; }

	/** @return array containing points of this vector ([V2::DIMENSION_X], [V2::DIMENSION_Y]) */
	const TYPE * getDimensions() const { return &x; }
	/** @return getDimensions()[a_dimensionField] */
	const TYPE getField(const int a_dimensionField) const { return (&x)[a_dimensionField]; }
	/** getDimensions()[a_dimensionField]=a_value; */
	void setField(const int a_dimensionField, const TYPE a_value) { (&x)[a_dimensionField]=a_value; }
	/** @param a_axis which axis to flip around {DIMENSION_X, DIMENSION_Y} */
	void flipAxis(const int & a_axis){ (&x)[a_axis]*=-1; }
	/** @return a version of this flipped around the given axis */
	V2<TYPE> flippedAxis(int a_axis){ V2<TYPE> p(*this); p.flipAxis(a_axis); return p;}

	/** resets the value of this vector */
	void set(const TYPE a_x, const TYPE a_y){x = a_x;	y = a_y;}
	/** copy another vector */
	void set(const V2<TYPE> & a_v2d){set(a_v2d.x, a_v2d.y);}
	/** make this a cos/sin of the given degree radian */
	void set(const TYPE a_piRadians){x=cos(a_piRadians);y=sin(a_piRadians);}

	/** default constructor */
	V2():x(0),y(0){}
	/** complete constructor */
	V2(TYPE a_x, TYPE a_y):x(a_x),y(a_y){}
	/** de-serialization constructor */
	V2(TYPE * a_twoValues):x(a_twoValues[0]),y(a_twoValues[1]){}
	/** turns a pi-radians angle into a vector */
	explicit V2(TYPE a_piRadians):x(cos(a_piRadians)), y(sin(a_piRadians)){}
	/** copy constructor */
	V2(const V2<TYPE> & v):x(v.x),y(v.y){}

	/** sets x and y to zero */
	void zero(){ x=y=0; }

	/**
	 * declares a "global" variable in a function, which is OK in a template! 
	 * (can't declare globals in headers otherwise)
	 */
	static const V2<TYPE> & ZERO()			{static V2<TYPE> ZERO(0,0);        return ZERO;}
	static const V2<TYPE> & ZERO_DEGREES()	{static V2<TYPE> ZERODEGREES(1,0);	return ZERODEGREES;}

	/**
	 * @return true if both x and y are zero
	 * @note function is not allowed to modify members in V2
	 */
	bool isZero() const {return x == 0 && y == 0;}

	/** @return if both x and y are less than the given x and y */
	bool isLessThan(const V2<TYPE> & p) const { return x < p.x && y < p.y; }
	/** @return if both x and y are greaterthan or equl to the given x and y */
	bool isGreaterThanOrEqualTo(const V2<TYPE> & p) const { return x >= p.x && y >=p.y; }

	/** @return the squared length of the vector (avoids sqrt) "quadrance" */
	TYPE magnitudeSq() const { return (TYPE)(x*x+y*y); }

	/** @return the length of the vector (uses sqrt) */
	TYPE magnitude() const { return sqrt((TYPE)magnitudeSq()); }

	/** @return dot product of this and v2. how "aligned" are these two? 0 = perp, +1 = same dir, -1 = opposite dir */
	TYPE dot(const V2<TYPE> & v2) const { return (TYPE)(x * v2.x + y * v2.y); }

	/** @return dot product of v1 and v2. how "aligned" are these two? 0 = perp, +1 = same dir, -1 = opposite dir */
	static TYPE dot(const V2<TYPE> & v1, const V2<TYPE> & v2) { return (TYPE)(v1.x * v2.x + v1.y * v2.y); }

	/**
	 * @return positive if position p is clockwise of this vector 
	 * (assume Y points down, X to right)
	 */
	TYPE sign(const V2<TYPE> & p) const { return (x*p.y)-(y*p.x); }

	/** @return if this point is clockwise of line a->b */
	bool isCW(V2<TYPE> const & a, V2<TYPE> const & b) const {
		return difference(a).sign(b.difference(a)) <= 0;
	}
	/** @return if this point is counter-clockwise of line a->b */
	bool isCCW(V2<TYPE> const & a, V2<TYPE> const & b) const {
		return difference(a).sign(b.difference(a)) >= 0;
	}

	/** @return true if this point is inside the given triangle */
	bool isInsideTriangle(V2<TYPE> const & a, V2<TYPE> const & b, V2<TYPE> const & c) const {
		TYPE signab = difference(a).sign(b.difference(a)), 
				signbc = difference(b).sign(c.difference(b)), 
				signac = difference(c).sign(a.difference(c));
		return(((signab>=0) == (signbc>=0)) && ((signbc>=0) == (signac>=0)))
			||(((signab<=0) == (signbc<=0)) && ((signbc<=0) == (signac<=0)));
	}

	/** triangle points are in clock-wise order */
	bool isInsideTriangleCW(V2<TYPE> const & a, V2<TYPE> const & b, V2<TYPE> const & c) const {
		return isCW(a,b) && isCW(b,c) && isCW(c,a);
	}

	/** triangle points are in counter-clock-wise order */
	bool isInsideTriangleCCW(V2<TYPE> const & a, V2<TYPE> const & b, V2<TYPE> const & c) const {
		return isCCW(a,b) && isCCW(b,c) && isCCW(c,a);
	}

	/** @return if the given polygon has all of it's points in clock-wise order. false if CCW or not convex */
	static bool isPolyCW(V2<TYPE> * const & points, int points_length) {
		for(int i = 0; i < points_length; ++i) {
			if(!points[(i+2)%points_length].isCW(points[i], points[(i+1)%points_length]))
				return false;
		}
		return true;
	}

	/** @return if the given polygon has all of it's points in counter-clock-wise order. false if CW or not convex */
	static bool isPolyCCW(V2<TYPE> * const & points, int points_length) {
		for(int i = 0; i < points_length; ++i) {
			if(!points[(i+2)%points_length].isCCW(points[i], points[(i+1)%points_length]))
				return false;
		}
		return true;
	}

	/** @return if this point is inside the given convex polygon, who's points are in clock-wise order */
	bool isInsidePolyCW(V2<TYPE> * const & points, int points_length) const {
		for(int i = 0; i < points_length; ++i) {
			if(isCW(points[i], points[(i+1)%points_length]))
				return false;
		}
		return true;
	}

	/** @return if this point is inside the given convex polygon, who's points are in counter-clock-wise order */
	bool isInsidePolyCCW(V2<TYPE> * const & points, int points_length) const {
		for(int i = 0; i < points_length; ++i){
			if(isCCW(points[i], points[(i+1)%points_length]))
				return false;
		}
		return true;
	}

	/** @return true if this point is between points A and B */
	bool isBetween(V2<TYPE> const & a, V2<TYPE> const & b)const{
		V2<TYPE> bracket = b.difference(a);
		bracket.setPerp();
		TYPE signa = difference(a).sign(bracket), 
				signb = difference(b).sign(bracket); 
		return ((signa >= 0) ^ (signb >= 0));
	}
	
	/** @return the vector perpendicular to this one */
	V2<TYPE> perp() const { return V2(-y, x); }

	/** make this vector perpendicular to itself */
	void setPerp() { TYPE t = x; x = -y; y = t; }

	/** @param a_length what to set this vector's magnitude to */
	void setMagnitude(const TYPE a_length) {
		TYPE m = magnitude();
		// "x = (x * a_length) / m" is more precise than "x *= (a_length/m)"
		x = (x * a_length) / m;
		y = (y * a_length) / m;
	}

	/**
	 * @param a_maxLength x and y adjusts so that the magnitude does not exceed this
	 */
	void truncate(const TYPE a_maxLength) {
		TYPE max2 = a_maxLength*a_maxLength;
		TYPE mm = magnitudeSq();
		// do check in quadrance space, to save a sqrt operation
		if(mm > max2) {
			TYPE m = sqrt(mm);
			// "x = (x * a_maxLength) / l" is more precise than "x *= (a_maxLength/l)"
			x = (x * a_maxLength) / m;
			y = (y * a_maxLength) / m;
		}
	}

	/** @return the quadrance (distance squared) between this vector and v */
	TYPE distanceSq(const V2<TYPE> & v) const {
		TYPE dx = x-v.x, dy = y-v.y;
		return dx*dx+dy*dy;
	}

	/** @return the pythagorean distance between this vector and v */
	TYPE distance(const V2<TYPE> & v) const{ return sqrt(distanceSq(v)); }

	static TYPE distance(const V2<TYPE> & a, const V2<TYPE> & b) { return a.distance(b); }

	static TYPE distanceSq(const V2<TYPE> & a, const V2<TYPE> & b) { return a.distanceSq(b); }

	/** @return the manhattan distance between this vector and v */
	TYPE distanceManhattan(const V2<TYPE> & v) const { return abs(v.x-x)+abs(v.y-y); }

	/** @return the vector that is the reverse of this vector, same as multiplying by -1 */
	V2<TYPE> getReverse() const { return V2<TYPE>(-x, -y); }

	/** make this direction vector go in the opposite direction */
	void setReverse() { *this *= -1; }

	/** @return a new V2<TYPE> that is the sum(+) of this V2<TYPE> and v */
	V2<TYPE> sum(const V2<TYPE> & v) const { return V2<TYPE>(x+v.x, y+v.y); }

	/** @return a new V2<TYPE> that is the difference(-) of this V2<TYPE> and v */
	V2<TYPE> difference(const V2<TYPE> & v) const { return V2<TYPE>(x-v.x, y-v.y); }

	/** @return a new V2<TYPE> that is the product(*) of this V2<TYPE> and v */
	V2<TYPE> product(const V2<TYPE> & v) const { return V2<TYPE>(x*v.x, y*v.y); }

	/** @return a new V2<TYPE> that is the product(*) of this V2<TYPE> and v */
	V2<TYPE> product(const TYPE a_value) const { return V2<TYPE>(x*a_value, y*a_value); }

	/** @return a new V2<TYPE> that is the quotient of this V2<TYPE> and the given V2*/
	V2<TYPE> quotient(const TYPE a_value) const { return V2<TYPE>(x/a_value, y/a_value); }

	/** @return a new V2<TYPE> that is the quotient of this V2<TYPE> and the given V2<TYPE> */
	V2<TYPE> quotient(const V2<TYPE> & v) const { return V2<TYPE>(x/v.x, y/v.y); }

	/** @return this V2<TYPE> after adding v */
	V2<TYPE> & add(const V2<TYPE> & v) {
		x += v.x;
		y += v.y;
		return *this;
	}

	/** @return this V2<TYPE> after subtracting v */
	V2<TYPE> & subtract(const V2<TYPE> & v) {
		x -= v.x;
		y -= v.y;
		return *this;
	}

	/** @return this V2<TYPE> after multiplying v */
	V2<TYPE> & multiply(const TYPE a_value) {
		x *= a_value;
		y *= a_value;
		return *this;
	}

	/** @return this V2<TYPE> after multiplying v */
	V2<TYPE> & multiply(const V2<TYPE> & v) {
		x *= v.x;
		y *= v.y;
		return *this;
	}

	/** @return this V2<TYPE> after dividing v */
	V2<TYPE> & divide(const TYPE a_value) {
		x /= a_value;
		y /= a_value;
		return *this;
	}

	/** @return this V2<TYPE> after dividing v */
	V2<TYPE> & divide(const V2<TYPE> & v) {
		x /= v.x;
		y /= v.y;
		return *this;
	}

	/** @return if this V2<TYPE> is euqal to v */
	bool isEqual(const V2<TYPE> & v) const { return (x == v.x && y == v.y); }

	/** @return true if this point is within a_radius from a_point */
	bool isWithin(const TYPE a_radius, const V2<TYPE> & a_point) const {
		TYPE rr = a_radius*a_radius;
		return distanceSq(a_point) <= rr; // calculate in quadrance space to avoid sqrt
	}

	/** @return if this point is between the Axis-Aligned Bounding Box (rectangle) inscribed by the given corners */
	bool isBetweenAABB(const V2<TYPE> & a, const V2<TYPE> & b) {
		return((a.x <= b.x)?(x>=a.x&&x<=b.x):(x>=b.x&&x<=a.x))
			&&((a.y <= b.y)?(y>=a.y&&y<=b.y):(y>=b.y&&y<=a.y));
	}

	/** forces vector to have a length of 1 */
	V2<TYPE> & normalize() {
		divide(magnitude());
		return *this;
	}

	/** a normalized verson of this vector */
	V2<TYPE> normal() const {
		V2<TYPE> norm(*this);
		return norm.normalize();
	}

	/** @return radians between these normalized vector is */
	TYPE piRadians(const V2<TYPE> & v) const { return V_COS(dot(v)); }

	/** @return radians that this normalized vector is, using ZERO_DEGREES() as the starting point */
	TYPE piRadians() const { return (y<0?-1:1)*piRadians(ZERO_DEGREES()); }
	
	/** @return how many degrees (standard 360 degree scale) this is */
	TYPE degrees() const { return ((y > 0)?-1:1)*acos(x)*180/V_PI; }

	/** @param a_normal cos(theta), sin(theta) as x,y values */
	void rotate(const V2<TYPE> & a_normal) {
		TYPE len = magnitude();	// remember length data
		// normalize()	// turn vector into simple angle, lose length data
		divide(len);	// same as normalize, but avoids calculating length again
		// calculate turn in-place within a new data structure
		V2<TYPE> turned(
			// x_ = x*cos(theta) - y*sin(theta)
			x*a_normal.x - y*a_normal.y, 
			// y_ = x*sin(theta) + y*cos(theta)
			x*a_normal.y + y*a_normal.x);
		// memory copy of structure
		*this = turned;
		// put length data back into normalized vector
		multiply(len);
	}

	/**
	 * @param a_normal cos(theta),sin(theta) as x,y values 
	 * NOTE: only use if this is a unit vector!
	 */
	void rotateUnitVectors(const V2<TYPE> a_normal) {
		// x_ = x*cos(theta) - y*sin(theta)
		TYPE x0 = x*a_normal.x - y*a_normal.y; 
		// y_ = x*sin(theta) + y*cos(theta)
		TYPE y0 = x*a_normal.y + y*a_normal.x;
		x = x0;
		y = y0;
	}

	/** @param a_degreePiRadians in piRadians */
	void rotate(const TYPE a_degreePiRadians) {
		rotate(V2<TYPE>(cos(a_degreePiRadians), sin(a_degreePiRadians)));
	}

	/**
	 * Used to undo toLocalSpace transformations
	 *
	 * @param a_xAxis rotate this point's x axis to match this vector
	 * @param a_yAxis rotate this point's y axis to match this vector
	 */
	V2<TYPE> toWorldSpace(V2<TYPE> const & a_xAxis, V2<TYPE> const & a_yAxis,
		const V2<TYPE> & a_pos) const
	{
		return V2<TYPE>(
			(a_xAxis.x*x) + (a_yAxis.x*y) + (a_pos.x),
			(a_xAxis.y*x) + (a_yAxis.y*y) + (a_pos.y));
	}

	/**
	 * Used to undo toLocalSpace transformations
	 *
	 * @param a_xAxis rotate this point's x axis to match this vector
	 * @param a_yAxis rotate this point's y axis to match this vector
	 */
	V2<TYPE> toWorldSpace(V2<TYPE> const & a_xAxis, V2<TYPE> const & a_yAxis,
		const V2<TYPE> & a_pos, const V2<TYPE> & a_scale) const
	{
		return V2<TYPE>(
			(a_scale.x*a_xAxis.x*x) + (a_scale.y*a_yAxis.x*y) + (a_pos.x),
			(a_scale.x*a_xAxis.y*x) + (a_scale.y*a_yAxis.y*y) + (a_pos.y));
	}

	/**
	 * Used to undo toLocalSpace transformations (angle only, for unit vectors)
	 *
	 * @param a_xAxis rotate this point's x axis to match this vector
	 * @param a_yAxis rotate this point's y axis to match this vector
	 */
	V2<TYPE> toWorldSpace(V2<TYPE> const & a_xAxis, V2<TYPE> const & a_yAxis) const
	{
		return V2<TYPE>((a_xAxis.x*x) + (a_yAxis.x*y), (a_xAxis.y*x) + (a_yAxis.y*y));
	}

	/**
	 * Used when it is more convenient to rotate one point to match a system
	 * rather than rotate a system to match a point. For example, when testing
	 * collision of a rotated Axis-Aligned Bounding-Box
	 *
	 * @param a_xVector vector of translated x dimension (V2f(1,0).rotate(rotationAngle))
	 * @param a_yVector vector of translated y dimension (V2f(0,1).rotate(rotationAngle))
	 * @param a_origin origin to translate this point in relation to (center of the system)
	 */
	V2<TYPE> toLocalSpace(V2<TYPE> const & a_xVector, V2<TYPE> const & a_yVector, V2<TYPE> const & a_origin) const
	{
		TYPE tx = -a_origin.dot(a_xVector);
		TYPE ty = -a_origin.dot(a_yVector);
		return V2<TYPE>(
			(a_xVector.x*x) + (a_xVector.y*y) + (tx),
			(a_yVector.x*x) + (a_yVector.y*y) + (ty));
	}

	static void convertToLocalSpace(const V2<TYPE> * points, const int numPoints, V2<TYPE> const & a_xVector, 
		V2<TYPE> const & a_yVector, V2<TYPE> const & a_origin) {
		TYPE tx = -a_origin.dot(a_xVector);
		TYPE ty = -a_origin.dot(a_yVector);
		for (int i = 0; i < numPoints; ++i) {
			points[i] = V2<TYPE>(
				(a_xVector.x*points[i].x) + (a_xVector.y*points[i].y) + (tx),
				(a_yVector.x*points[i].x) + (a_yVector.y*points[i].y) + (ty))
		}
	}

	/**
	 * organizes a list of 2D vectors into a circular curve (arc)
	 *
	 * @param a_startVector what to start at (normalized vector)
	 * @param a_angle the angle to increment the arc by V2(piRadians), where piRadians is totalArcPiRadians/a_arcPoints
	 * @param a_list the list of 2D vectors to map along the arc
	 * @param a_arcs the number of elements in a_list
	 */
	static void arc(V2<TYPE> const & a_startVector, V2<TYPE> const & a_angle, V2<TYPE> * const & a_list, int const & a_arcPoints) {
		TYPE len = a_startVector.magnitude();	// remember magnitude
		a_list[0] = a_startVector;				// copy starting point for calculations
		a_list[0].divide(len);					// normalize starting point
		V2<TYPE> * lastPoint = &a_list[0];	// faster memory reference than a_list[i-1]
		// calculate all points in the arc as normals (faster: no division)
		for (int i = 1; i < a_arcPoints; ++i) {
			// calculate rotation in next point
			(lastPoint+1)->set(
				// x_ = x*cos(theta) - y*sin(theta)
				lastPoint->x*a_angle.x - lastPoint->y*a_angle.y, 
				// y_ = x*sin(theta) + y*cos(theta)
				lastPoint->x*a_angle.y + lastPoint->y*a_angle.x);
			++lastPoint;
		}
		if(len != 1) {
			// put length data back into normalized vector
			for (int i = 0; i < a_arcPoints; ++i) {
				// embarassingly parallel
				a_list[i].multiply(len);
			}
		}
	}

	/**
	 * ensures wraps this V2's x/y values around the given Axis-Aligned Bounding Box range (like a torroid)
	 * @param a_min the minimum x/y
	 * @param a_max the maximum x/y
	 */
	void wrapAround(V2<TYPE> const & a_min, V2<TYPE> const & a_max) {
		TYPE width = a_max.x - a_min.x;
		TYPE height= a_max.y - a_min.y;
		while(x < a_min.x){	x += width;	}
		while(x > a_max.x){	x -= width;	}
		while(y < a_min.y){	y +=height;	}
		while(y > a_max.y){	y -=height;	}
	}

	/** @return the position half-way between line a->b */
	static V2<TYPE> between(V2<TYPE> const & a, V2<TYPE> const & b) {
		V2<TYPE> average = b.sum(a);
		average.divide(2.0);
		return average;
	}

	/**
	 * @param percentage 0 is a, 1 is b, .5 is between(a,b). any numeric value should work.
	 * @return an interpolated point between points a and b.
	 */
	static V2<TYPE> lerp(V2<TYPE> const & a, V2<TYPE> const & b, const float percentage) {
		V2<TYPE> delta = b.difference(a);
		TYPE len = delta.magnitude();
		delta.setMagnitude(len * percentage);
		return a.sum(delta);
	}

	/**
	 * @param A,B line 1
	 * @param C,D line 2
	 * @param point __OUT to the intersection of line AB and CD
	 * @param dist __OUT the distance along line AB to the intersection
	 * @return true if intersection occurs between the lines
	 */
	static bool lineIntersection(const V2<TYPE> & A, const V2<TYPE> & B, 
								const V2<TYPE> & C, const V2<TYPE> & D, 
								TYPE & dist, V2<TYPE> & point)
	{
		TYPE r, s;
		if (rayIntersection(A, B, C, D, r, s)) {
			V2<TYPE> ray = (B - A) * r;
			dist = ray.magnitude();
			point = A + ray;
			return ( (r > 0) && (r < 1) && (s > 0) && (s < 1) );
		}
		return false;
	}

	/**
	* @param A,B line 1
	* @param C,D line 2
	* @param out_alongAB __OUT how far along AB (percentage) the collision happened
	* @param out_alongCD __OUT how far along CD (percentage) the collision happened
	* @return false if parallel, true otherwise
	*/
	static bool rayIntersection(const V2<TYPE> & A, const V2<TYPE> & B,
		const V2<TYPE> & C, const V2<TYPE> & D, 
		TYPE & out_alongAB, TYPE & out_alongCD)
	{
		TYPE rTop = (A.y - C.y)*(D.x - C.x) - (A.x - C.x)*(D.y - C.y);
		TYPE rBot = (B.x - A.x)*(D.y - C.y) - (B.y - A.y)*(D.x - C.x);
		TYPE sTop = (A.y - C.y)*(B.x - A.x) - (A.x - C.x)*(B.y - A.y);
		TYPE sBot = (B.x - A.x)*(D.y - C.y) - (B.y - A.y)*(D.x - C.x);
		if ((rBot == 0) || (sBot == 0)) {
			//lines are parallel
			return false;
		}
		out_alongAB = rTop / rBot; // how far along AB the collision happened
		out_alongCD = sTop / sBot; // how far along CD the collision happened
		return true;
	}

	/**
	 * Does ray-casting against one polygon by detecting line collision
	 *
	 * @param rayStart
	 * @param rayDirection must be a unit vector
	 * @param polygonLoop a list of points
	 * @param polygonLoopCount how many points
	 * @param out_dist how far along the ray the collision happened
	 * @param out_point where the collision happened
	 * @param out_normal the normal of the collision line
	 * 
	 * @return -1 if no collision happened, othewise, the loop segment where the collision happened
	 */
	static int raycastPolygon(V2<TYPE> const & rayStart, V2<TYPE> const & rayDirection,
		const V2<TYPE> * polygonLoop, const int polygonLoopCount, float & out_dist,
		V2<TYPE> & out_point, V2<TYPE> & out_normal) {
		V2<TYPE> rayEnd = rayStart + rayDirection;
		float surfacePercent, closestDistance, closestSurfacePercent;
		int closestEdge = - 1;
		for (int i = 0; i < polygonLoopCount; ++i) {
			int j = (i + 1) % polygonLoopCount;
			if (rayIntersection(rayStart, rayEnd, polygonLoop[i], polygonLoop[j], out_dist, surfacePercent)) {
				if (out_dist > 0 && surfacePercent > 0 && surfacePercent < 1
				&& (closestEdge == -1 || out_dist < closestDistance)) {
					closestEdge = i;
					closestDistance = out_dist;
					closestSurfacePercent = surfacePercent;
				}
			}
		}
		if (closestEdge != -1) {
			int other = (closestEdge + 1) % polygonLoopCount;
			V2<TYPE> delta = (polygonLoop[other] - polygonLoop[closestEdge]);
			out_point = polygonLoop[closestEdge] + delta * closestSurfacePercent;
			out_normal = delta.perp().normal();
		}
		return closestEdge;
	}

	/**
	 * @param a_out_closestPoint will be closest point to a_point on line AB
	 * @return true if a_out_closestPoint is actually on line AB
	 */
	static bool closestPointOnLine(const V2<TYPE> & A, const V2<TYPE> & B, 
		const V2<TYPE> & a_point, V2<TYPE> & a_out_closestPoint)
	{
		V2<TYPE> r = B.difference(A).perp();
		V2<TYPE> d = r.product(2).sum(a_point);
		V2<TYPE> c = a_point.difference(r);
		TYPE dist;
		bool intersected = lineIntersection(A, B, c, d, dist, a_out_closestPoint);
		return intersected;
	}

	/** @return if circle (a_point,a_radius) crosses line (A,B) */
	static bool lineCrossesCircle(const V2<TYPE> & A, const V2<TYPE> & B, 
		const V2<TYPE> & a_point, const TYPE a_radius, V2<TYPE> & a_out_closePoint)
	{
		bool connectionOnLine = closestPointOnLine(A, B, a_point, a_out_closePoint);
		return(connectionOnLine && a_out_closePoint.distance(a_point) <= a_radius)
			|| a_point.distance(A) < a_radius
			|| a_point.distance(B) < a_radius;
	}

	/**
	* Does ray-casting against one circle by detecting line collision
	*
	* @param rayStart
	* @param rayDirection must be a unit vector
	* @param center center of the circle
	* @param radius radius of the circle
	* @param out_dist how far along the ray the collision happened
	* @param out_point where the collision happened
	* @param out_normal the normal of the collision line
	*
	* @return -1 if no collision happened, othewise, the loop segment where the collision happened
	*/
	static bool raycastCircle(V2<TYPE> const & rayStart, V2<TYPE> const & rayDirection,
		V2<TYPE> const center, const TYPE radius, float & out_dist,
		V2<TYPE> & out_point, V2<TYPE> & out_normal) {
		V2<TYPE> radiusDirection = rayDirection.perp();
		TYPE radiusNeededToHit;
		if (rayIntersection(
			rayStart, rayStart + rayDirection,
			center, center + radiusDirection, out_dist, radiusNeededToHit)) {
			if (out_dist > 0 && abs(radiusNeededToHit) <= radius) {
				V2<TYPE> linesMeet = rayStart + rayDirection * out_dist;
				// pythagorean theorum to get missing triangle side
				float b = sqrt(radius*radius - radiusNeededToHit*radiusNeededToHit);
				V2<TYPE> hitDelta = rayDirection * b;
				// find which of the two solutions is closer to the ra start
				V2<TYPE> p1 = linesMeet + hitDelta;
				float d1 = V2<TYPE>::distance(rayStart, p1);
				V2<TYPE> p2 = linesMeet - hitDelta;
				float d2 = V2<TYPE>::distance(rayStart, p2);
				if (d1 < d2) {
					out_point = p1;
					out_dist = d1;
				} else {
					out_point = p2;
					out_dist = d2;
				}
				out_normal = (out_point - center).normal();
				return true;
			}
		}
		return false;
	}

	/**
	* @param velocity the vector to be reflected
	* @param normal Vector2 that is normal to the tangent being reflected (not the wall, but the wall's normal)
	* @param out_reflected will be set to the reflected velocity
	* @return the dot-product, which is useful to determine if the velocity
	* is going against the normal (positive value, expected if going into the wall) or with the
	* normal (negative value, probably an error)
	*/
	static TYPE calculateReflection(V2<TYPE> const & velocity, V2<TYPE> const & normal,
		V2<TYPE> & out_reflected){
		TYPE velocityDotProduct = V2<TYPE>::dot(normal, velocity);
		out_reflected.set(velocity.x - 2 * velocityDotProduct * normal.x,
			velocity.y - 2 * velocityDotProduct * normal.y);
		return velocityDotProduct;
	}

	// overloaded operators
	V2<TYPE> operator-()const{return this->getReverse();}
	V2<TYPE> operator+(V2<TYPE> const & rhs)const{return sum(rhs);}
	V2<TYPE> operator-(V2<TYPE> const & rhs)const{return difference(rhs);}
	V2<TYPE> operator*(V2<TYPE> const & rhs)const{ return product(rhs); }
	V2<TYPE> operator/(V2<TYPE> const & rhs)const{ return quotient(rhs); }
	V2<TYPE> operator*(TYPE const rhs)const{ return product(rhs); }
	V2<TYPE> operator/(TYPE const rhs)const{ return quotient(rhs); }
	bool operator==(V2<TYPE> const & rhs) const{return isEqual(rhs);}
	bool operator!=(V2<TYPE> const & rhs) const{return !isEqual(rhs);}
	V2<TYPE> operator+=(V2<TYPE> const & rhs){ return add(rhs); }
	V2<TYPE> operator-=(V2<TYPE> const & rhs){ return subtract(rhs); }
	V2<TYPE> operator*=(V2<TYPE> const & rhs){ return multiply(rhs); }
	V2<TYPE> operator/=(V2<TYPE> const & rhs){ return divide(rhs); }
	V2<TYPE> operator*=(TYPE const rhs){ return multiply(rhs); }
	V2<TYPE> operator/=(TYPE const rhs){ return divide(rhs); }

#ifdef __GL_H__
// OpenGL specific functions

	/** calls glVertex2fv on this data structure */
	void glVertex() const;

	/**
	 * translates the open GL rendering context.
	 * @note: Dont forget to push and pop the matrix!
	 */
	void glTranslate() const { glTranslatef((GLfloat)x, (GLfloat)y, (GLfloat)0); }

	/**
	 * scale the open GL rendering context.
	 * @note: Dont forget to push and pop the matrix!
	 */
	void glScale() const { glScalef((GLfloat)x, (GLfloat)y, (GLfloat)1); }

	/**
	 * rotate the open GL rendering context.
	 * @note: Dont forget to push and pop the matrix!
	 */
	void glRotate() const { glRotatef(-degrees(), 0, 0, 1); }

	/**
	 * draw an OpenGL line from (0,0) to (V2.x, V2.y)
	 * @return true if something was drawn
	 */
	void glDraw() const {
		glBegin(GL_LINES);
		ZERO().glVertex();
		glVertex();
		glEnd();
	}

	bool glDrawTo(const V2<TYPE> & a_next) const {
		glBegin(GL_LINES);
		glVertex();
		a_next.glVertex();
		glEnd();
		return true;
	}

	/**
	 * @param GLenum_mode draws using this GL enum. eg: GL_LINES, GL_LINE_LOOP, GL_LINE_STRIP
	 */
	static void glVertexList(V2<TYPE> * a_list, const int & a_count) {
		for(int i = 0; i < a_count; ++i) {
			a_list[i].glVertex();
		}
	}
#endif
};

// V2 using float, useful for all kinds of games and simulations
typedef V2<float> V2f;
// V2 using double, useful for high-precision simulations
typedef V2<double> V2d;
// V2 using int, useful for grids, tables, row/col coordinates
typedef V2<int> V2i;
// V2 using short, useful for limited grids, tables, row/col coordinates
typedef V2<short> V2s;

#ifdef __GL_H__
// specialized template implementations. inlined so that the linker wont complain about redefinition
inline void V2f::glVertex() const { glVertex2fv((float*)this); }
inline void V2d::glVertex() const { glVertex2dv((double*)this); }
inline void V2i::glVertex() const { glVertex2iv((int*)this); }
inline void V2s::glVertex() const { glVertex2sv((short*)this); }
#endif
