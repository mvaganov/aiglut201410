#include <GL/freeglut.h>
#include "glutrenderingcontext.h"

#include <stdio.h> // for printf

/**
* @param a_dimensionPixels pixel width/height
* @param a_min cartesian minimum
* @param a_max cartesian maximum
*/
GLUTRenderingContext::GLUTRenderingContext(const V2f & a_dimensionPixels,
	const V2f & a_min, const V2f & a_max)
	:dimensionPixels(a_dimensionPixels), dimensionCartesian(a_min, a_max), glutFont(GLUT_BITMAP_HELVETICA_10)
{}

void GLUTRenderingContext::setColor(long fourByteColor) {
	glColor3ubv((GLubyte*)&fourByteColor);
}
void GLUTRenderingContext::drawLine(const float aX, const float aY, const float bX, const float bY) {
	glBegin(GL_LINES);
	glVertex2f(aX, aY);
	glVertex2f(bX, bY);
	glEnd();
}

void GLUTRenderingContext::drawLine(V2f const & a, V2f const & b) {
	glBegin(GL_LINES);
	a.glVertex();
	b.glVertex();
	glEnd();
}
void GLUTRenderingContext::drawCircle(Circf const & c, bool filled) {
	c.glDraw(filled);
}
void GLUTRenderingContext::drawCircle(V2f const & center, const float radius, const bool filled) {
	glDrawCircle(center, radius, filled);
}
void GLUTRenderingContext::drawRect(AABBf const & r, const bool filled) {
	r.glDraw(filled);
}
void GLUTRenderingContext::drawRect(V2f const & min, V2f const & max, const bool filled) {
	AABBf(min, max).glDraw(filled);
}
void GLUTRenderingContext::drawBox(Boxf const & b, const bool filled) {
	b.glDraw(filled);
}
void GLUTRenderingContext::drawBox(V2f const & center, V2f const & size, const float rotationInRadians, const bool filled) {
	Boxf(center, size, rotationInRadians).glDraw(filled);
}

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
	dimensionCartesian.min *= percentageChange;
	dimensionCartesian.max *= percentageChange;
	V2f whereMouseIsAfterZoom = percentOfScreenWhereMouseIs * dimensionCartesian.getDimension();
	V2f delta = centeredHere - whereMouseIsAfterZoom;
	scroll(delta);
	ortho2D();
}

/**
 * draw the cartesian plane 
 * @param a_gridLineSize how big to make the unit marks (sample value: 0.1f).
 */
void GLUTRenderingContext::drawPlanarAxis(float const a_gridLineSize) {
	glBegin(GL_LINES);
	V2f(0, dimensionCartesian.min.y).glVertex();
	V2f(0, dimensionCartesian.max.y).glVertex();
	V2f(dimensionCartesian.min.x, 0).glVertex();
	V2f(dimensionCartesian.max.x, 0).glVertex();
	if (a_gridLineSize != 0) {
		int i, min, max;
		min = (int)(dimensionCartesian.min.x - 1);
		max = (int)(dimensionCartesian.max.x + 1);
		for (i = min; i < max; ++i)	{
			V2f((float)i, a_gridLineSize).glVertex();
			V2f((float)i, -a_gridLineSize).glVertex();
		}
		min = (int)(dimensionCartesian.min.y - 1);
		max = (int)(dimensionCartesian.max.y + 1);
		for (i = min; i < max; ++i)	{
			V2f(-a_gridLineSize, (float)i).glVertex();
			V2f(a_gridLineSize, (float)i).glVertex();
		}
	}
	glEnd();
}

void GLUTRenderingContext::drawGrid(V2f const gridScale) {
	glBegin(GL_LINES);
	float divs, extra, i;
	float minx = (dimensionCartesian.min.x);
	float maxx = (dimensionCartesian.max.x);
	float miny = (dimensionCartesian.min.y);
	float maxy = (dimensionCartesian.max.y);
	divs = minx / gridScale.x;
	extra = divs - (int)divs;
	for (i = minx - extra * gridScale.x; i < maxx; i += gridScale.x)	drawLine(i, miny, i, maxy);
	divs = miny / gridScale.y;
	extra = divs - (int)divs;
	for (i = miny - extra * gridScale.y; i < maxy; i += gridScale.y)	drawLine(minx, i, maxx, i);
	glEnd();
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
void GLUTRenderingContext::print(const V2f position, const char * a_string) {
	float lineHeight = 30;
	switch ((int)glutFont) {
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
				glutBitmapCharacter(glutFont, a_string[index]);
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

int GLUTRenderingContext::printf(const V2f position, const char *fmt, ...)
{
	int bufferSize = 20, charsPrinted;
	char * buffer;
	do {
		va_list valist;
		va_start(valist, fmt);
#if __STDC_WANT_SECURE_LIB__
		charsPrinted = _vscprintf(fmt, valist);
		bufferSize = charsPrinted + 2;
		buffer = (char*)malloc(bufferSize);
		vsnprintf_s(buffer, bufferSize, _TRUNCATE, fmt, valist);
#else
		buffer = (char*)malloc(bufferSize);
		charsPrinted = vsnprintf(buffer, bufferSize, fmt, valist);
		if (charsPrinted >= bufferSize) {
			free(buffer);
			bufferSize *= 2;
		}
#endif
		va_end(valist);
	} while (charsPrinted >= bufferSize);
	this->print(position, buffer);
	free(buffer);
	return charsPrinted;
}
