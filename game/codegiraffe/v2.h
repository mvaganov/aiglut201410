#pragma once
#include <GL/freeglut.h>
#include <math.h>
#include "m.h"
#include "random.h"

// TODO move polygon-related code into Polygon2

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

	/** const access X */
	TYPE getX() const {return x;}
	/** const access Y */
	TYPE getY() const {return y;}
	/** mutate X */
	void setX(const TYPE a_value) { x = a_value; }
	/** mutate Y */
	void setY(const TYPE a_value) { y = a_value; }
	/** mutate in bulk */
	void add(const TYPE dx, const TYPE dy) { x += dx; y += dy; }

	/** @return thisAsAnArray[a_dimensionField] */
	const TYPE getField(const int a_dimensionField) const { return (&x)[a_dimensionField]; }
	/** getDimensions()[a_dimensionField]=a_value; */
	void setField(const int a_dimensionField, const TYPE a_value) { (&x)[a_dimensionField]=a_value; }
	/** @param a_axis which axis to flip around {DIMENSION_X, DIMENSION_Y} */
	void flipAxis(const int a_axis){ (&x)[a_axis]*=-1; }
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
	explicit V2(TYPE * a_twoValues):x(a_twoValues[0]),y(a_twoValues[1]){}
	/** turns a pi-radians angle into a vector. explicit so that scalars are not all-of-the-sudden turned into vectors */
	explicit V2(TYPE a_piRadians):x((TYPE)cos(a_piRadians)), y((TYPE)sin(a_piRadians)){}
	/** copy constructor */
	V2(const V2<TYPE> & v):x(v.x),y(v.y){}

	/** sets x and y to zero */
	void setZero() { x=y=0; }

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
	TYPE magnitude() const { return (TYPE)sqrt((TYPE)magnitudeSq()); }

	/** @return dot product of v1 and v2. how "aligned" are these two? 0 = perp, +1 = same dir, -1 = opposite dir */
	static TYPE dot(V2<TYPE> const & v1, V2<TYPE> const & v2) { return (TYPE)(v1.x * v2.x + v1.y * v2.y); }

	/**
	 * @return positive if position p is clockwise of this vector, 0 if the same. same as dot product of a vector perpendicular to the given vector
	 * (assume Y points down, X to right)
	 */
	TYPE sign(const V2<TYPE> & p) const { return (x*p.y)-(y*p.x); }

	/** @return if this point is clockwise of line a->b (assume Y points down, X to right) */
	bool isCW(V2<TYPE> const & a, V2<TYPE> const & b) const {
		return (*this - a).sign(b - a) <= 0;
	}
	/** @return if this point is counter-clockwise of line a->b (assume Y points down, X to right) */
	bool isCCW(V2<TYPE> const & a, V2<TYPE> const & b) const {
		return (*this - a).sign(b - a) >= 0;
	}

	/** triangle points are in clock-wise order (assume Y points down, X to right) */
	bool isInsideTriangleCW(V2<TYPE> const & a, V2<TYPE> const & b, V2<TYPE> const & c) const {
		return isCW(a,b) && isCW(b,c) && isCW(c,a);
	}

	/** triangle points are in counter-clock-wise order (assume Y points down, X to right) */
	bool isInsideTriangleCCW(V2<TYPE> const & a, V2<TYPE> const & b, V2<TYPE> const & c) const {
		return isCCW(a,b) && isCCW(b,c) && isCCW(c,a);
	}

	/** @return true if this point is inside the given triangle */
	bool isInsideTriangle(V2<TYPE> const & a, V2<TYPE> const & b, V2<TYPE> const & c) const {
		TYPE signab = (*this - a).sign(b - a),
			signbc = (*this - b).sign(c - b),
			signac = (*this - c).sign(a - c);
		return(((signab >= 0) == (signbc >= 0)) && ((signbc >= 0) == (signac >= 0)))
			|| (((signab <= 0) == (signbc <= 0)) && ((signbc <= 0) == (signac <= 0)));
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
	bool isInsidePolyCW(const V2<TYPE> * const & points, const int points_length) const {
		for(int i = 0; i < points_length; ++i) {
			if(isCW(points[i], points[(i+1)%points_length]))
				return false;
		}
		return true;
	}

	/** @return if this point is inside the given convex polygon, who's points are in counter-clock-wise order */
	bool isInsidePolyCCW(const V2<TYPE> * const & points, const int points_length) const {
		for(int i = 0; i < points_length; ++i){
			if(isCCW(points[i], points[(i+1)%points_length]))
				return false;
		}
		return true;
	}
	
	/** @return the vector perpendicular to this one */
	V2<TYPE> perp() const { return V2(-y, x); }

	/** make this vector perpendicular to itself */
	void setPerp() { TYPE t = x; x = -y; y = t; }

	/** @return true if this point is between points A and B */
	bool isBetween(V2<TYPE> const & a, V2<TYPE> const & b)const{
		V2<TYPE> bracket = b - a;
		bracket.setPerp();
		TYPE signa = (*this - a).sign(bracket),
			signb = (*this - b).sign(bracket);
		return ((signa >= 0) ^ (signb >= 0));
	}

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
		TYPE m = magnitude();
		if (m > a_maxLength) {
			// "x = (x * a_maxLength) / l" is more precise than "x *= (a_maxLength/l)"
			x = (x * a_maxLength) / m;
			y = (y * a_maxLength) / m;
		}
	}

	static TYPE distance(const V2<TYPE> & a, const V2<TYPE> & b) { return (a - b).magnitude(); }

	TYPE distance(const V2<TYPE> & v) const { return (*this - v).magnitude(); }

	static TYPE distanceSq(const V2<TYPE> & a, const V2<TYPE> & b) { return (a - b).magnitudeSq(); }

	/** @return the manhattan distance between this vector and v */
	static TYPE distanceManhattan(const V2<TYPE> & a, const V2<TYPE> & b) { return abs(b.x - a.x) + abs(b.y - a.y); }

	/** @return the vector that is the reverse of this vector, same as multiplying by -1 */
	V2<TYPE> getReverse() const { return V2<TYPE>(-x, -y); }

	/** make this direction vector go in the opposite direction */
	void setReverse() { *this *= -1; }

	/** @return a new V2<TYPE> that is the negative value of this V2<TYPE> */
	V2<TYPE> operator-() const { return V2<TYPE>(-x, -y); }
	/** @return a new V2<TYPE> that is the sum(+) of this V2<TYPE> and v */
	V2<TYPE> operator+(V2<TYPE> const & v) const { return V2<TYPE>(x + v.x, y + v.y); }
	/** @return a new V2<TYPE> that is the difference(-) of this V2<TYPE> and v */
	V2<TYPE> operator-(V2<TYPE> const & v) const { return V2<TYPE>(x - v.x, y - v.y); }
	/** @return a new V2<TYPE> that is the product(*) of this V2<TYPE> and v */
	V2<TYPE> operator*(V2<TYPE> const & v) const { return V2<TYPE>(x * v.x, y * v.y); }
	/** @return a new V2<TYPE> that is the product(*) of this V2<TYPE> and v */
	V2<TYPE> operator*(TYPE const a_scalar) const { return V2<TYPE>(x * a_scalar, y * a_scalar); }
	/** @return a new V2<TYPE> that is the quotient of this V2<TYPE> and the given scalar */
	V2<TYPE> operator/(TYPE const a_scalar) const { return V2<TYPE>(x / a_scalar, y / a_scalar); }
	/** @return a new V2<TYPE> that is the quotient of this V2<TYPE> and the given vector */
	V2<TYPE> operator/(V2<TYPE> const & v) const { return V2<TYPE>(x / v.x, y / v.y); }
	/** @return this V2<TYPE> after adding v */
	V2<TYPE> & operator+=(V2<TYPE> const & v){ x += v.x; y += v.y; return *this; }
	/** @return this V2<TYPE> after subtracting v */
	V2<TYPE> & operator-=(V2<TYPE> const & v){ x -= v.x; y -= v.y; return *this; }
	/** @return this V2<TYPE> after multiplying v */
	V2<TYPE> & operator*=(V2<TYPE> const & v){ x *= v.x; y *= v.y; return *this; }
	/** @return this V2<TYPE> after multiplying the given scalar */
	V2<TYPE> & operator*=(TYPE const a_scalar){ x *= a_scalar; y *= a_scalar; return *this; }
	/** @return this V2<TYPE> after dividing v */
	V2<TYPE> & operator/=(TYPE const a_scalar){ x /= a_scalar; y /= a_scalar; return *this; }
	/** @return this V2<TYPE> after dividing the given scalar */
	V2<TYPE> & operator/=(V2<TYPE> const & v){ x /= v.x; y /= v.y; return *this; }
	/** @return if this V2<TYPE> is equal to v */
	bool operator==(V2<TYPE> const & v) const { return (x == v.x && y == v.y); }
	/** @return if this V2<TYPE> is not equal to v */
	bool operator!=(V2<TYPE> const & v) const { return (x != v.x || y != v.y); }

	V2<TYPE> sum(V2<TYPE> const & v) const { return operator+(v); }
	V2<TYPE> difference(const V2<TYPE> & v) const { return operator-(v); }
	V2<TYPE> product(const V2<TYPE> & v) const { return operator*(v); }
	V2<TYPE> product(const TYPE a_scalar) const { return operator*(a_scalar); }
	V2<TYPE> quotient(const TYPE a_scalar) const { return operator/(a_scalar); }
	V2<TYPE> quotient(const V2<TYPE> & v) const { return operator/(v); }
	V2<TYPE> & add(const V2<TYPE> & v) { return operator+=(v);  }
	V2<TYPE> & subtract(const V2<TYPE> & v) { return operator-=(v);  }
	V2<TYPE> & multiply(const TYPE a_scalar) { return operator*=(a_scalar); }
	V2<TYPE> & multiply(const V2<TYPE> & v) { return operator*=(v); }
	V2<TYPE> & divide(const TYPE a_scalar) { return operator/=(a_scalar); }
	V2<TYPE> & divide(const V2<TYPE> & v) { return operator/=(v); }
	bool isEqual(const V2<TYPE> & v) const { return operator==(v); }

	/** @return true if this point is within a_radius from a_point */
	bool isWithin(const TYPE a_radius, const V2<TYPE> & a_point) const {
//		TYPE rr = a_radius*a_radius;
//		return distanceSq(a_point) <= rr; // calculate in quadrance space to avoid sqrt
		return V2<TYPE>::distance(*this, a_point) < a_radius;
	}

	/** @return if this point is between the Axis-Aligned Bounding Box (rectangle) inscribed by the given corners */
	bool isBetweenAABB(V2<TYPE> const & a, V2<TYPE> const & b) const {
		return((a.x <= b.x)?(x>=a.x&&x<=b.x):(x>=b.x&&x<=a.x))
			&&((a.y <= b.y)?(y>=a.y&&y<=b.y):(y>=b.y&&y<=a.y));
	}

	/** forces vector to have a length of 1 */
	V2<TYPE> & normalize() {
		divide(magnitude());
		return *this;
	}

	/** forces vector to have a length of 1 if the vector is not a zero vector */
	V2<TYPE> & normalizeIfNotZero() {
		if (!isZero()) divide(magnitude());
		return *this;
	}

	/** a normalized verson of this vector */
	V2<TYPE> normal() const {
		V2<TYPE> norm(*this);
		return norm.normalize();
	}

	/** @return radians between these normalized vector is (accounts for negative rotation) */
	TYPE piRadians(V2<TYPE> const & v) const { return ((sign(v) < 0) ? 1 : -1) * (TYPE)acos(dot(*this, v)); }

	/** @return radians that this normalized vector is, using ZERO_DEGREES() as the starting point */
	TYPE piRadians() const { return ((y > 0) ? 1 : -1) * (TYPE)acos(x); }
	
	/** @return how many degrees (standard 360 degree scale) this is */
	TYPE degrees() const { return ((y > 0) ? 1 : -1) * (TYPE)acos(x) * 180 / (TYPE)V_PI; }

	/** @param a_normal cos(theta), sin(theta) as x,y values */
	void rotate(const V2<TYPE> & a_normal) {
		if (isZero()) return;
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
		this->operator*=(len);
	}

	/**
	 * @param a_normal cos(theta),sin(theta) as x,y values 
	 * NOTE: only use if this is a unit vector!
	 * NOTE: also, interestingly, this is the same as ComplexNumber multiplication 
	 * (real and imaginary parts, like "4 + 5i"). x is the real part, y is imaginary
	 */
	void rotateUnitVectors(const V2<TYPE> a_normal) {
		// x_ = x*cos(theta) - y*sin(theta)
		TYPE x0 = x*a_normal.x - y*a_normal.y; 
		// y_ = x*sin(theta) + y*cos(theta)
		TYPE y0 = x*a_normal.y + y*a_normal.x;
		set(x0, y0);
	}

	/**
	 * modifies the vector
	 * @param a_degreePiRadians in piRadians
	 */
	void rotate(const TYPE a_degreePiRadians) {
		rotate(V2<TYPE>((TYPE)cos(a_degreePiRadians), (TYPE)sin(a_degreePiRadians)));
	}

	/**
	 * does not modify the calling vector, it returns a vector just like this one after a rotation
	 * @param a_degreePiRadians in piRadians
	 */
	V2<TYPE> rotated(const TYPE a_degreePiRadians) const {
		V2<TYPE> result(*this);
		if (!result.isZero()) result.rotate(V2<TYPE>((TYPE)cos(a_degreePiRadians), (TYPE)sin(a_degreePiRadians)));
		return result;
	}

	/**
	 * does not modify the calling vector, it returns a vector just like this one after a rotation
	 * @param a_normal cos(theta), sin(theta) as x,y values
	 */
	V2<TYPE> rotated(const V2<TYPE> & a_normal) const {
		V2<TYPE> result(*this);
		result.rotate(a_normal);
		return result;
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
	 * Used to undo toLocalSpace transformations (angle only, for unit vectors most likely)
	 *
	 * @param a_xAxis rotate this point's x axis to match this vector
	 * @param a_yAxis rotate this point's y axis to match this vector
	 */
	V2<TYPE> toWorldSpace(V2<TYPE> const & a_xAxis, V2<TYPE> const & a_yAxis) const {
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
	V2<TYPE> toLocalSpace(V2<TYPE> const & a_xVector, V2<TYPE> const & a_yVector, V2<TYPE> const & a_origin) const {
		TYPE tx = dot(-a_origin, a_xVector);
		TYPE ty = dot(-a_origin, a_yVector);
		return V2<TYPE>(
			(a_xVector.x*x) + (a_xVector.y*y) + (tx),
			(a_yVector.x*x) + (a_yVector.y*y) + (ty));
	}

	static void convertToLocalSpace(const V2<TYPE> * points, const int numPoints, V2<TYPE> const & a_xVector, 
		V2<TYPE> const & a_yVector, V2<TYPE> const & a_origin) {
		TYPE tx = dot(-a_origin, a_xVector);
		TYPE ty = dot(-a_origin, a_yVector);
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
				a_list[i] *= (len);
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

	void clampToInt() { x = (TYPE)((int)x); y = (TYPE)((int)y); }
	void roundToInt() { x = (TYPE)((int)(x + 0.5)); y = (TYPE)((int)(y + 0.5)); }

	/** @return the position half-way between line a->b */
	static V2<TYPE> between(V2<TYPE> const & a, V2<TYPE> const & b) { return (b + a) / 2; }

	/**
	 * @param percentage 0 is a, 1 is b, .5 is between(a,b). any numeric value should work.
	 * @return an interpolated point between points a and b.
	 */
	static V2<TYPE> lerp(V2<TYPE> const & a, V2<TYPE> const & b, const TYPE percentage) {
		V2<TYPE> delta = b - a;
		delta *= percentage;
		return a + delta;
	}

	/**
	 * @param a
	 * @param b
	 * @param c
	 * @param out_circumcenter __OUT the circumcenter, if there is one
	 * @param out_radius __OUT the radius of the inscribed circle
	 * @return if there is a circle that this triangle is circumscribed into
	 */
	static bool circumcenter(V2<TYPE> const & a, V2<TYPE> const & b, V2<TYPE> const & c,
		V2<TYPE> & out_circumcenter, TYPE & out_radius) {
		V2<TYPE> ab = b - a, bc = c - b;
		V2<TYPE> midAB = a + ab / 2, midBC = b + bc / 2;
		V2<TYPE> bisectorAB = ab.perp(), bisectorBC = bc.perp();
		float rayMidAB, rayMidBC;
		if (V2<TYPE>::rayIntersection(midAB, midAB + bisectorAB, midBC, midBC + bisectorBC, rayMidAB, rayMidBC)) {
			out_circumcenter = midAB + bisectorAB * rayMidAB;
			out_radius = (out_circumcenter - a).magnitude();
			return true;
		}
		return false;
	}

	/**
	 * @param A start of line AB
	 * @param B end of line AB
	 * @param C start of line CD
	 * @param D end of line CD
	 * @param point __OUT to the intersection of line AB and CD
	 * @param dist __OUT the distance along line AB to the intersection
	 * @return true if intersection occurs between the lines
	 */
	static bool lineIntersection(const V2<TYPE> & A, const V2<TYPE> & B, 
								const V2<TYPE> & C, const V2<TYPE> & D, 
								TYPE & dist, V2<TYPE> & point)
	{
		TYPE ab, cd;
		if (rayIntersection(A, B, C, D, ab, cd)) {
			V2<TYPE> ray = (B - A) * ab;
			dist = ray.magnitude();
			point = A + ray;
			return ( (ab > 0) && (ab < 1) && (cd > 0) && (cd < 1) );
		}
		return false;
	}

	/**
	* @param A start of line AB
	* @param B end of line AB
	* @param C start of line CD
	* @param D end of line CD
	* @return true if intersection occurs between the lines, and on the lines
	*/
	static bool lineIntersection(const V2<TYPE> & A, const V2<TYPE> & B, const V2<TYPE> & C, const V2<TYPE> & D)
	{
		TYPE ab, cd;
		if (rayIntersection(A, B, C, D, ab, cd)) {
			return ((ab > 0) && (ab < 1) && (cd > 0) && (cd < 1));
		}
		return false;
	}

	/**
	* @param A start of line AB
	* @param B end of line AB
	* @param C start of line CD
	* @param D end of line CD
	* @param out_alongAB __OUT how far along AB (percentage) the collision happened
	* @param out_alongCD __OUT how far along CD (percentage) the collision happened
	* @return false if parallel, true otherwise
	*/
	static bool rayIntersection(const V2<TYPE> & A, const V2<TYPE> & B, const V2<TYPE> & C, const V2<TYPE> & D, 
		TYPE & out_alongAB, TYPE & out_alongCD)
	{
		V2<TYPE> deltaAB = B - A, deltaCD = D - C, deltaOrigins = A - C;
		TYPE alignmentDifference = deltaAB.sign(deltaCD); // checks the dot product of perpendicular angles
		if (alignmentDifference == 0) { return false; } // fail if the lines are parallel
		out_alongAB = deltaCD.sign(deltaOrigins) / alignmentDifference; // how far along AB the collision happened
		out_alongCD = deltaAB.sign(deltaOrigins) / alignmentDifference; // how far along CD the collision happened
		return true;
		//TYPE abTop = (A.y - C.y)*(D.x - C.x) - (A.x - C.x)*(D.y - C.y);
		//TYPE abBot = (B.x - A.x)*(D.y - C.y) - (B.y - A.y)*(D.x - C.x);
		//TYPE cdTop = (A.y - C.y)*(B.x - A.x) - (A.x - C.x)*(B.y - A.y);
		//TYPE cdBot = (B.x - A.x)*(D.y - C.y) - (B.y - A.y)*(D.x - C.x);
		//if ((abBot == 0) || (cdBot == 0)) {
		//	return false; // lines are parallel
		//}
		//out_alongAB = abTop / abBot; // how far along AB the collision happened
		//out_alongCD = cdTop / cdBot; // how far along CD the collision happened
		//return true;
	}

	/**
	 * Does ray-casting against one polygon by detecting line collision TODO move this into Polygon, so that Ray and RaycastHit can be used
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
				if (out_dist >= 0 && surfacePercent >= 0 && surfacePercent <= 1
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
			out_dist = closestDistance;
			out_normal = delta.perp().normal();
		}
		return closestEdge;
	}

	/**
	 * @param a_out_closestPoint will be closest point to a_point on line AB
	 * @return true if a_out_closestPoint is actually on line AB, false if parallel or not on line
	 */
	static bool closestPointOnLine(const V2<TYPE> & A, const V2<TYPE> & B, 
		const V2<TYPE> & a_point, V2<TYPE> & a_out_closestPoint)
	{
		V2<TYPE> perp = (B - A).perp();
		float percentOfAB, percentOfCD;
		bool intersected = rayIntersection(A, B, a_point, a_point + perp, percentOfAB, percentOfCD);
		if (!intersected)
			return false;
		intersected = percentOfAB > 0 && percentOfAB < 1; // true only if point is on the given line
		a_out_closestPoint = A + ((B - A) * percentOfAB);
		return intersected;
	}

	/** @return if circle (a_point,a_radius) crosses the given ray */
	static bool rayCrossesCircle(const V2<TYPE> & rayStart, const V2<TYPE> & rayDirection,
		const V2<TYPE> & a_center, const TYPE a_radius)
	{
		V2<TYPE> radiusDirection = rayDirection.perp();
		radiusDirection.normalize();
		float rayDist, radiusDist;
		rayIntersection(rayStart, rayDirection, a_center, a_center + radiusDirection, rayDist, radiusDist);
		return (radiusDist <= a_radius) && rayDist > 0;
	}
	/** @return if circle (a_point,a_radius) crosses line (A,B) */
	static bool lineCrossesCircle(const V2<TYPE> & A, const V2<TYPE> & B, 
		const V2<TYPE> & a_center, const TYPE a_radius, V2<TYPE> & a_out_closePoint)
	{
		V2<TYPE> radiusDirection = (B - A).perp();
		radiusDirection.normalize();
		float lineDist, radiusDist;
		rayIntersection(A, B, a_center, a_center + radiusDirection, lineDist, radiusDist);
		a_out_closePoint = radiusDirection * radiusDist + a_center;
		return (abs(radiusDist) <= a_radius) && (lineDist > 0 && lineDist < 1);
	}

	/**
	 * @param point the point you are looking for neighbors to
	 * @param list multiple points
	 * @param listSize how many points are in list
	 * @return the index (from list) of the point that is closest to point
	 */
	static int indexOfClosest(V2<TYPE> const & point, const V2<TYPE> * const & list, const int listSize) {
		if (listSize < 0) return -1;
		int closestIndex = 0;
		V2<TYPE> delta = list[0] - point;
		TYPE mag, closestSoFar = delta.magnitude();
		for (int i = 1; i < listSize; ++i) {
			delta = list[i] - point;
			mag = delta.magnitude();
			if (mag < closestSoFar) {
				closestSoFar = mag;
				closestIndex = i;
			}
		}
		return closestIndex;
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
	* @return false if no collision happened
	*/
	static bool raycastCircle(V2<TYPE> const & rayStart, V2<TYPE> const & rayDirection,
		V2<TYPE> const center, const TYPE radius, float & out_dist,
		V2<TYPE> & out_point, V2<TYPE> & out_normal) {
		V2f p1, p2;
		if (raycastCircle(rayStart, rayDirection, center, radius, p1, p2)) {
			V2<TYPE> dp1 = (p1 - rayStart), dp2 = (p2 - rayStart);
			// check if the points are in the correct direction
			bool p1good = V2<TYPE>::dot(rayDirection, dp1) > 0;
			bool p2good = V2<TYPE>::dot(rayDirection, dp2) > 0;
			float dist1 = dp1.magnitude(), dist2 = dp2.magnitude();
			// if both are good, invalidate the further one
			if (p1good && p2good) {
				if (dist1 < dist2) {
					p2good = false;
				} else {
					p1good = false;
				}
			}
			// if at least one is good...
			if (p1good || p2good) {
				// give the one that is good
				if (p1good) {
					out_dist = dist1;
					out_point = p1;
				} else {
					out_dist = dist2;
					out_point = p2;
				}
				out_normal = (out_point - center).normal();
				return true;
			}
		}
		return false;
	}

	/**
	* Does ray-casting against one circle by detecting line collision TODO move to circle, so that Ray and RaycastHit can be used
	*
	* @param rayStart
	* @param rayDirection must be a unit vector
	* @param center center of the circle
	* @param radius radius of the circle
	* @param out_p1 one of the points on the line and the circle
	* @param out_p2 the other point on the line and the circle (may be the same as out_p1)
	*
	* @return false if no collision happened
	*/
	static bool raycastCircle(V2<TYPE> const & rayStart, V2<TYPE> const & rayDirection,
		V2<TYPE> const center, const TYPE radius, V2<TYPE> & out_p1, V2<TYPE> & out_p2) {
		V2f rayNorm = rayDirection.normal();
		if (rayStart == center) {
			out_p1 = out_p2 = center + rayNorm * radius;
			return true;
		}
		V2<TYPE> radiusDirection = rayNorm.perp();
		TYPE radHit, rayHit;
		if (rayIntersection(rayStart, rayStart + rayNorm, center, center + radiusDirection, rayHit, radHit)) {
			if (abs(radHit) <= radius) {
				V2<TYPE> intersection = center + radiusDirection * radHit;
				// pythagorean theorum to get missing triangle side, where radius is the hypontinuse
				float side = (float)sqrt(radius*radius - radHit*radHit);
				V2<TYPE> hitDelta = rayDirection * side;
				// find which of the two solutions is closer to the ra start
				out_p1 = intersection + hitDelta;
				out_p2 = intersection - hitDelta;
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

	static V2<TYPE> randomUnitVector() {
		V2<TYPE> r;
		do {
			r.x = Random::PRNGf(-1, 1);
			r.y = Random::PRNGf(-1, 1);
		} while (r.isZero());
		return r.normal();
	}

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
	void glRotate() const { glRotatef(degrees(), 0, 0, 1); }

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

	void glDrawTo(const V2<TYPE> & a_next) const {
		glBegin(GL_LINES);
		glVertex();
		a_next.glVertex();
		glEnd();
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

template<typename TYPE>
struct RaycastHit_ {
	V2<TYPE> point, normal;
	TYPE distance;
	RaycastHit_() {}
	RaycastHit_(TYPE const distance):distance(distance) {}
	RaycastHit_(V2<TYPE> const & point, V2<TYPE> const & normal, const TYPE distance) :point(point), normal(normal), distance(distance) {}
	void set(V2<TYPE> const & point, V2<TYPE> const & normal, TYPE const & distance) {
		this->point = point;
		this->normal = normal;
		this->distance = distance;
	}

	void setFromLine(V2<TYPE> const & from, V2<TYPE> const & to) {
		point = to;
		// done backwards here because normal is a surface normal that ends the ray from->to.
		normal = from - to;
		distance = normal.magnitude();
		if (distance != 0)
			normal /= distance;
	}
};

typedef RaycastHit_<float> RaycastHit;

template<typename TYPE>
struct Ray_ {
	V2<TYPE> start, direction;
	Ray_(V2<TYPE> const & start, V2<TYPE> const & direction) :start(start), direction(direction){}
	Ray_(){}

	/**
	 * @param ray another ray for this ray to intersect with
	 * @param out_distance how far this ray needs to travel to reach the given ray as a percentage of the ray.direction length
	 * @return false if there is no collision possible (the rays are parallel
	 */
	bool intersection(Ray_<TYPE> const & ray, TYPE & out_distance) const {
		V2<TYPE> deltaOrigins = start - ray.start;
		TYPE alignmentDifference = this->direction.sign(ray.direction); // checks the dot product of perpendicular angles
		if (alignmentDifference == 0) { return false; } // fail if the lines are parallel
		out_distance = ray.direction.sign(deltaOrigins) / alignmentDifference; // how far along AB the collision happened
		return true;
	}

	/**
	* @param ray another ray for this ray to intersect with
	* @param out_distance how far this ray needs to travel to reach the given ray as a percentage of the ray.direction length
	* @return false if there is no collision possible (the rays are parallel
	*/
	bool intersection(Ray_<TYPE> const & ray, TYPE & out_distance, TYPE & out_rayDistance) const {
		V2<TYPE> deltaOrigins = start - ray.start;
		TYPE alignmentDifference = this->direction.sign(ray.direction); // checks the dot product of perpendicular angles
		if (alignmentDifference == 0) { return false; } // fail if the lines are parallel
		out_distance = ray.direction.sign(deltaOrigins) / alignmentDifference; // how far along AB the collision happened
		out_rayDistance = direction.sign(deltaOrigins) / alignmentDifference; // how far along CD the collision happened
		return true;
	}

	void glDraw() const { start.glDrawTo(start + direction); }

	TYPE distanceFrom(V2<TYPE> const & p) const {
		V2<TYPE> perp = direction.perp();
		TYPE dist;
		intersection(Ray(p, perp), dist);
		return dist;
	}
};

typedef Ray_<float> Ray;
