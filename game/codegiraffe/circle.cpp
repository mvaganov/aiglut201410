#include <GL/freeglut.h>
#include "v2.h"

bool g_circleArcCalclated = false;
/** a circle will be drawn with this many outer segments. reducing this value may improve performance */
const int g_circlePoints = 32;
V2f g_circle[g_circlePoints];

void calculateCircleArc() {
	V2f::arc(V2f::ZERO_DEGREES(), V2f((float)V_2PI / g_circlePoints), g_circle, g_circlePoints);
}

bool glDrawCircle(const float a_x, const float a_y, const float a_radius, bool filled) {
	if (!g_circleArcCalclated) {
		calculateCircleArc();
		g_circleArcCalclated = true;
	}
	// move the matrix so that drawing a unit circle is actually exactly what we want
	glPushMatrix();
	glTranslatef(a_x, a_y, 0);
	glScalef(a_radius, a_radius, a_radius);
	// draw the circle
	if (filled) V2f::ZERO().glVertex(); // a triangle fan needs the center point first
	glBegin(filled?GL_TRIANGLE_FAN:GL_LINE_LOOP);
	V2f::glVertexList(g_circle, g_circlePoints);
	glEnd();
	glPopMatrix(); // move the matrix back
	return true;
}