#pragma once

#include "obstacles.h"
#include "codegiraffe/glutrenderingcontext.h"

class Agent {
public:
	/** where is this? */
	CircleObject body;
	/** where is this going? */
	V2f velocity;
	/** where is this facing? (cos, sin) of the angle this is facing */
	V2f direction;
	/** a 4 byte color value */
	long color;

	Agent(CircF circle):body(circle),direction(V2f::ZERO_DEGREES()){}

	void update(int a_ms) {
		body.center += velocity * (float)a_ms / 1000.0f;
	}
	void draw(GLUTRenderingContext & g_screen) {
		g_screen.setColor(color);
		body.glDraw(false);
		g_screen.drawLine(body.center, body.center + direction * body.radius);
	}
};