#pragma once

#include "Vector_.h"

template<typename Scalar>
class Angle2D_ {
	Vector_<_2D<Scalar> > startVector, endVector;
	Scalar startAngle, endAngle;
public:
	Vector_<_2D<Scalar> > origin;

	Scalar GetAngle() const { return endAngle - startAngle; }
	Scalar GetStartAngle() const { return startAngle; }
	Scalar GetEndAngle() const { return endAngle; }

	Angle2D_(const Vector_<_2D<Scalar> >& origin, const Scalar startAngle, const Scalar endAngle) 
		: startVector(UnitVector(startAngle)), endVector(UnitVector(endAngle)), startAngle(startAngle), endAngle(endAngle), origin(origin) {}

	/** @return if the given point p is inside of this angle */
	bool IsInAngle(const Vector_<_2D<Scalar> >& p) const {
		Vector_<_2D<Scalar> > relativeP = p - origin;
		Scalar totalAngle = GetAngle();
		if (totalAngle < M_PI) {
			// if the angle is less than PI, relativeP should be (CW to start and CCW to end)
			return Sign(startVector, relativeP) <= 0 && Sign(endVector, relativeP) >= 0;
		} else {
			// if the angle is greater than PI, relativeP should NOT be (CCW to start and CW to end)
			return !(Sign(startVector, relativeP) >= 0 && Sign(endVector, relativeP) <= 0);
		}
	}

	static Vector_<_2D<Scalar> > UnitVector(const Scalar piRadians) { return Vector_<_2D<Scalar> >(cos(piRadians), sin(piRadians)); }

	static const Vector_<_2D<Scalar> >& ZERO() { Vector_<_2D<Scalar> > zero_degrees(1, 0); return zero_degrees; }

	/** @return positive if position p is clockwise of vector v, 0 if the same. same as dot product of a vector perpendicular to the given vector */
	static Scalar Sign(const Vector_<_2D<Scalar> >& v, const Vector_<_2D<Scalar> >& p) { return (v.x*p.y) - (v.y*p.x); }

	/** @return if this point p is clockwise (or directly on top of) of line a->b */
	static bool IsCW(const Vector_<_2D<Scalar> >& a, const Vector_<_2D<Scalar> >& b, const Vector_<_2D<Scalar> >& p) { return Sign(b - a, p - a) <= 0; }

	/** @return if this point p is counter-clockwise (or directly on top of) of line a->b */
	static bool IsCCW(const Vector_<_2D<Scalar> >& a, const Vector_<_2D<Scalar> >& b, const Vector_<_2D<Scalar> >& p) { return Sign(b - a, p - a) >= 0; }

	/** @return radians between normalized vectors Vector_::ZERO->a and Vector_::ZERO->b (accounts for negative rotation) */
	static Scalar PiRadians(const Vector_<_2D<Scalar> >& a, const Vector_<_2D<Scalar> >& b) { return ((Sign(a, b) < 0) ? 1 : -1) * acos(Dot(a, b)); }

	/** @return radians that this normalized vector v is, using Angle2D::DEGREES() as the starting point */
	static Scalar PiRadians(const Vector_<_2D<Scalar> >& v) { return ((v.y > 0) ? 1 : -1) * acos(v.x); }

	/** @return how many degrees (standard 360 degree scale) this is */
	static Scalar Degrees(const Vector_<_2D<Scalar> >& v) { return PiRadians(v) * 180 / (Scalar)M_PI; }

	/**
	* @param p [in,out] the unit-vector to rotate (about the origin)
	* @param a_normal cos(theta),sin(theta) as x,y values
	* NOTE: only use if this is a unit vector!
	*/
	static void RotateUnitVector(const Vector_<_2D<Scalar> >& unitVector, const Vector_<_2D<Scalar> >& a_normal) {
		Set(x*a_normal.x - y*a_normal.y, x*a_normal.y + y*a_normal.x);
	}

	/**
	* @param p [in,out] the point to rotate (about the origin)
	* @param a_normal cos(theta), sin(theta) as x,y values
	*/
	static void Rotate(const Vector_<_2D<Scalar> >& p, const Vector_<_2D<Scalar> >& a_normal) {
		if (p.IsZero()) return;
		Scalar len = p.Magnitude();	// remember length data
		p /= len;	// same as normalize, but avoids calculating length again
		RotateUnitVector(p, a_normal);
		// put length data back into normalized vector
		p *= (len);
	}

	/**
	* @param v [in,out] a vector to be rotated (around Vector_::ZERO)
	* @param a_degreePiRadians in piRadians
	*/
	static void Rotate(Vector_<_2D<Scalar> >& v, const Scalar a_degreePiRadians) { Rotate(v, UnitVector(a_degreePiRadians)); }

	/**
	* @param unitVector [in,out] a unit vector to be rotated (around Vector_::ZERO)
	* @param a_degreePiRadians in piRadians
	*/
	static void RotateUnitVector(Vector_<_2D<Scalar> >& unitVector, const Scalar a_degreePiRadians) { RotateUnitVector(v, UnitVector(a_degreePiRadians)); }

	/**
	* @param v returns a vector just like this one, after a rotation
	* @param a_degreePiRadians in piRadians
	*/
	static Vector_<_2D<Scalar> > Rotated(const Vector_<_2D<Scalar> >& v, const Scalar a_degreePiRadians) {
		Vector_<_2D<Scalar> > result(*this);
		if (!result.IsZero()) Rotate(result, UnitVector(a_degreePiRadians));
		return result;
	}

	/**
	* @param v returns a vector just like this one, after a rotation
	* @param a_normal cos(theta), sin(theta) as x,y values
	*/
	static Vector_<_2D<Scalar> > Rotated(const Vector_<_2D<Scalar> >& v, const Vector_<_2D<Scalar> >& a_normal) {
		Vector_<_2D<Scalar> > result(*this);
		Rotate(result, a_normal);
		return result;
	}

};

typedef Angle2D_<float> Angle2f;