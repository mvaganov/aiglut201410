#include <GL/freeglut.h>
#include "glutrenderingcontext.h"

#include <stdio.h>

/** sets up the openGL screen according to current variables */
void GLUTRenderingContext::glSetupScreen() {
	// re-define the window size
	glutInitWindowSize((int)dimensionPixels.x,(int)dimensionPixels.y);
	// reset The Current Viewport
	glViewport(0,0,(int)dimensionPixels.x,(int)dimensionPixels.y);
}

/** sets up the cartesian plane for openGL/GLUT */
void GLUTRenderingContext::ortho2D() {
	// reset the view matrix
	glLoadIdentity();
	// set the view matrix to match this context
	gluOrtho2D(
		dimensionCartesian.getMinX(), dimensionCartesian.getMaxX(), 
		dimensionCartesian.getMinY(), dimensionCartesian.getMaxY());
}

/** 
 * @param a_position x/y in pixels (like from a mouse press) 
 * @return cartesian coordinate equivalent
 */
V2f GLUTRenderingContext::convertPixelsToCartesian(const V2f & a_position) const {
	V2f dim = dimensionCartesian.getDimension();
	// because the y axis is inverted...
	V2f pos(a_position.x, dimensionPixels.y-a_position.y);
	V2f cartesian = pos.product(dim);
	cartesian.divide(dimensionPixels);
	cartesian += dimensionCartesian.getMin();
	return cartesian;
}

/**
 * resizes the cartesian plane based on how the window has been resized
 */
void GLUTRenderingContext::reshape(int a_widthPixels, int a_heightPixels) {
	// check to see which rectangle edge moved
	int dLft = 0, dBot = 0, dRgt = 0, dTop = 0;
	static int init_x = glutGet(GLUT_INIT_WINDOW_X), init_y = glutGet(GLUT_INIT_WINDOW_Y);
	int pos_x = glutGet(GLUT_WINDOW_X), pos_y = glutGet(GLUT_WINDOW_Y);
	dLft = pos_x - init_x, dTop = pos_y - init_y; // check if the upper-left (minX,maxY) corner moved
	init_x = pos_x;	init_y = pos_y; // remember the old position
	// figure out percentage change of new size vs old size (in pixels)
	V2f newDimensionPixels((float)a_widthPixels, (float)a_heightPixels);
	dRgt = (int)(newDimensionPixels.x - dimensionPixels.x) + dLft;
	dBot = (int)(newDimensionPixels.y - dimensionPixels.y) + dTop;
	//printf("top: %d, left: %d, right: %d, bottom:%d\n", dTop, dLft, dRgt, dBot);
	V2f changePercent = newDimensionPixels / dimensionPixels;
	// use the new size (dimensionPixels) to setup the screen
	dimensionPixels = newDimensionPixels;
	glSetupScreen();
	// find how much the cartesian system must change based on the pixel system
	V2f newDimensionCarte(dimensionCartesian.getDimension() * changePercent);
	// calculate the net cartesian coordinate change
	V2f dDim = newDimensionCarte - dimensionCartesian.getDimension();
	float dMinX = (dDim.x == 0) ? 0 : -dDim.x * dLft / (dLft + dRgt);
	float dMaxX = (dDim.x == 0) ? 0 : dDim.x * dRgt / (dLft + dRgt);
	float dMinY = (dDim.y == 0) ? 0 : -dDim.y * dBot / (dTop + dBot);
	float dMaxY = (dDim.y == 0) ? 0 : dDim.y * dTop / (dTop + dBot);
	// adjust existing cartesian coordinate system based on calculated deltas
	dimensionCartesian.add(dMinX, dMinY, dMaxX, dMaxY);
	// reset the display
	ortho2D();
}

void GLUTRenderingContext::scroll(V2f const delta) {
	dimensionCartesian.add(delta.x, delta.y, delta.x, delta.y);
	// reset the display
	ortho2D();
}

/** this variable is used only in this compilation unit, to help keep track of drag-scrolling */
static V2f lastMousePosition;
bool scrolling = false;

/** change what the view area is looking at. called when drag scrolling begins */
void GLUTRenderingContext::scrollStart(V2f const mouseClick) {
	lastMousePosition = mouseClick;
	scrolling = true;
}

/** change what the view area is looking at. called when drag scrolling begins */
void GLUTRenderingContext::scrollStop(V2f const mouseClick) {
	lastMousePosition = mouseClick;
	scrolling = false;
}

/**
 * change what the view area is looking at. called as drag scrolling continues (just pass the mouse position)
 * will only scroll if scrollStart(V2f const) was called first
 */
void GLUTRenderingContext::scrollDrag(V2f const draggedMousePosition) {
	if (scrolling) {
		V2f delta = draggedMousePosition - lastMousePosition;
		scroll(-delta);
		lastMousePosition = draggedMousePosition - delta;
	}
}

/** zoom in to the given percentage of the current view, centered at the given location */
void GLUTRenderingContext::zoom(float percentageChange, V2f const centeredHere) {
	V2f dimension = dimensionCartesian.getDimension();
	V2f percentOfScreenWhereMouseIs = centeredHere / dimension;
	dimensionCartesian.m_min *= percentageChange;
	dimensionCartesian.m_max *= percentageChange;
	V2f whereMouseIsAfterZoom = percentOfScreenWhereMouseIs * dimensionCartesian.getDimension();
	V2f delta = centeredHere - whereMouseIsAfterZoom;
	scroll(delta);
	ortho2D();
}

/**
 * draw the cartesian plane 
 * @param a_dashSize how big to make the unit marks (sample value: 0.1)
 */
void GLUTRenderingContext::glDraw(float a_dashSize) {
	glBegin(GL_LINES);
	glVertex2f(0, dimensionCartesian.getMaxY());	V2f::ZERO().glVertex();
	glVertex2f(0, dimensionCartesian.getMinY());	V2f::ZERO().glVertex();
	glVertex2f(dimensionCartesian.getMinX(), 0);	V2f::ZERO().glVertex();
	glVertex2f(dimensionCartesian.getMaxX(), 0);	V2f::ZERO().glVertex();
	if(a_dashSize > 0) {
		int i;
		for(i = 1; i < dimensionCartesian.getMaxY(); i++)	{glVertex2f((GLfloat)-a_dashSize, (GLfloat)i);	glVertex2f((GLfloat)a_dashSize, (GLfloat)i);}
		for(i =-1; i > dimensionCartesian.getMinY(); i--)	{glVertex2f((GLfloat)-a_dashSize, (GLfloat)i);	glVertex2f((GLfloat)a_dashSize, (GLfloat)i);}
		for(i = 1; i < dimensionCartesian.getMaxX(); i++)	{glVertex2f((GLfloat)i, (GLfloat)-a_dashSize);	glVertex2f((GLfloat)i, (GLfloat)a_dashSize);}
		for(i =-1; i > dimensionCartesian.getMinX(); i--)	{glVertex2f((GLfloat)i, (GLfloat)-a_dashSize);	glVertex2f((GLfloat)i, (GLfloat)a_dashSize);}
	}
	glEnd();
}

/** @param a_string to draw at the given location */
void GLUTRenderingContext::glDrawString(const V2f position, const char * a_string) {
	return glDrawString(position, a_string, GLUT_BITMAP_HELVETICA_10);
}

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
void GLUTRenderingContext::glDrawString(const V2f position, const char * a_string, void * font) {
	float lineHeight = 30;
	switch ((int)font) {
	case (int)GLUT_BITMAP_9_BY_15:        lineHeight = 15;	break;
	case (int)GLUT_BITMAP_8_BY_13:        lineHeight = 13;	break;
	case (int)GLUT_BITMAP_TIMES_ROMAN_10: lineHeight = 10;	break;
	case (int)GLUT_BITMAP_TIMES_ROMAN_24: lineHeight = 24;	break;
	case (int)GLUT_BITMAP_HELVETICA_10:   lineHeight = 10;	break;
	case (int)GLUT_BITMAP_HELVETICA_12:   lineHeight = 12;	break;
	case (int)GLUT_BITMAP_HELVETICA_18:   lineHeight = 18;	break;
	}
	const int distanceBetweenLines = 3;
	lineHeight += distanceBetweenLines;

	// determine how to scale letter advancement based on how far moving 1x1 changes the 2D raster position
	V2f pixelsPerUnit = dimensionPixels / dimensionCartesian.getDimension();
	//printf("%f %f\n", pixelsPerUnit.x, pixelsPerUnit.y);
	lineHeight /= pixelsPerUnit.y;
	int linesPrinted = 0, index = 0;
	do {
		glRasterPos2f(position.x, position.y - linesPrinted*lineHeight);
		char c = a_string[index];
		while (c) {
			if (c != '\n' && c != '\r') {
				glutBitmapCharacter(font, a_string[index]);
			}
			index++;
			c = a_string[index];
			if (c == '\n') {
				++linesPrinted;
				break;
			}
		}
	} while (a_string[index]);
}