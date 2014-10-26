#include "v2.h"
#include "aabb.h"

/**
 * Handles data & methods related to screen size and cartesian-vs-pixel calculations
 * Intended to be used with OpenGL, and more specifically, freeglut http://freeglut.sourceforge.net/docs/api.php
 *
 * @author mvaganov@hotmail.com October 2014
 */
struct GLUTRenderingContext {
	/** how many pixels wide and tall the screen is */
	V2f dimensionPixels;
	RectF dimensionCartesian;

	/**
	 * @param a_dimensionPixels pixel width/height
	 * @param a_min cartesian minimum
	 * @param a_max cartesian maximum
	 */
	GLUTRenderingContext(const V2f & a_dimensionPixels, 
		const V2f & a_min, const V2f & a_max)
		:dimensionPixels(a_dimensionPixels), dimensionCartesian(a_min, a_max)
	{}

	/** sets up the openGL screen according to current variables */
	void glSetupScreen();

	/** sets up the cartesian plane for openGL/GLUT */
	void ortho2D();

	/** 
	 * @param a_position x/y in pixels (like from a mouse press) 
	 * @return cartesian coordinate equivalent
	 */
	V2f convertPixelsToCartesian(const V2f & a_position) const;

	/**
	 * resizes the cartesian plane based on how the window has been resized
	 */
	void reshape(int a_widthPixels, int a_heightPixels);

	/** change what the view area is looking at. add the given delta to the min/max coordinates */
	void scroll(V2f const delta);

	/** change what the view area is looking at. called when drag scrolling begins */
	void scrollStart(V2f const mouseClick);

	/** change what the view area is looking at. called when drag scrolling begins */
	void scrollStop(V2f const mouseClick);

	/** change what the view area is looking at. called as drag scrolling continues (just pass the mouse position) */
	void scrollDrag(V2f const draggedMousePosition);

	/** zoom in to the given percentage of the current view, centered at the given location */
	void zoom(float percentageChange, V2f const centeredHere);

	/**
	 * draw the cartesian plane 
	 * @param a_dashSize how big to make the unit marks (sample value: 0.1)
	 */
	void glDraw(float a_dashSize);

	/** @param a_string to draw at the given location */
	void GLUTRenderingContext::glDrawString(const V2f position, const char * a_string);

	/**
	* @param font one of these: {
	*            GLUT_BITMAP_9_BY_15,
	*            GLUT_BITMAP_8_BY_13,
	*            GLUT_BITMAP_TIMES_ROMAN_10,
	*            GLUT_BITMAP_TIMES_ROMAN_24,
	*            GLUT_BITMAP_HELVETICA_10,
	*            GLUT_BITMAP_HELVETICA_12,
	*            GLUT_BITMAP_HELVETICA_18
	* }
	*/
	void GLUTRenderingContext::glDrawString(const V2f position, const char * a_string, void * font);
};
