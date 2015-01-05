// properties->"VC++ Directories"->"Include Directories"->add "$(ProjectDir)\freeglut\include"
// properties->"VC++ Directories"->"Library Directories"->add "$(ProjectDir)\freeglut\lib"
// properties->"Debugging"->"Environment"->set to "PATH=%PATH%;$(ProjectDir)\freeglut\bin"

#include <stdio.h>
#include <GL/freeglut.h>
#include "codegiraffe/glutrenderingcontext.h"
#include "game.h"
#include "bullet.h"

/** the screen pixels, cartesian plane min, cartesian plane max */
GLUTRenderingContext g_screen(V2f(600, 600), V2f(-10.0, -10.0), V2f(10.0, 10.0));
/** whether or not to keep running the game loop */
bool g_mouseButtonsDown[3];
const int g_numMouseButtons = sizeof(g_mouseButtonsDown) / sizeof(g_mouseButtonsDown[0]);
Game g_game;

/**
 * @param key what keyboard key was pressed (ASCII codes)
 * @param x/y where the mouse was at the time
 * @note: the following may be helpful:<code>
 int state=glutGetModifiers();
 if (state & GLUT_ACTIVE_SHIFT)	// shift is pressed</code>
 */
void keyboard(unsigned char key, int x, int y) {
	g_game.key(g_screen.convertPixelsToCartesian(V2f((float)x, (float)y)), key, GLUT_DOWN);
	glutPostRedisplay();
}

/** @param a_width/a_height new dimensions of the resized window */
void reshape(int a_width, int a_height) {
	// have the graphics context calculate the changes...
	g_screen.reshape(a_width, a_height);
	glutPostRedisplay();
}

/**
* @param button {GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, GLUT_RIGHT_BUTTON, 3 (mousewheel up), 4 (mousewheel down)}
* @param state {GLUT_DOWN, GLUT_UP}
* @param x/y the coordinate of where the mouse is
*/
void mouse(int button, int state, int x, int y) {
	g_mouseButtonsDown[button] = (state == GLUT_DOWN) ? true : false;
	g_game.mouse(g_screen.convertPixelsToCartesian(V2f((float)x, (float)y)), button, state);
	glutPostRedisplay();
}

/** @param x/y the coordinate of where the mouse is */
void passiveMotion(int x, int y) {
	g_game.mouse(g_screen.convertPixelsToCartesian(V2f((float)x, (float)y)), -1, -1);
}

/** @param x/y the coordinate of where the mouse is being dragged */
void draggedMotion(int x, int y) {
	V2f drag = g_screen.convertPixelsToCartesian(V2f((float)x, (float)y));
	for (int i = 0; i < g_numMouseButtons; ++i) {
		if (g_mouseButtonsDown[i]) {
			g_game.mouse(drag, i, 2);
		}
	}
//	glutPostRedisplay();
}

/** how many milliseconds have passed since the last time update was called */
void update(int a_ms) {
	g_game.update(a_ms);
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);	//Clear the screen

	// insert stuff to draw here
	g_screen.setColor(0xeeeeee); // very light gray
	g_screen.drawGrid(V2f(4, 4));
	g_screen.setColor(0x888888); // gray
	g_screen.drawPlanarAxis(-0.2f);			// draw the cartisian plane

	g_game.draw(&g_screen);

	glFlush();						// print everything to the (off?) screen buffer
	glutSwapBuffers();				// swap the draw buffers
}

void init(const char * windowName) {
	// setup the screen using the graphics context
	g_screen.glSetupScreen();
	g_game.g_screen = &g_screen;
	// define the drawing mode
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	//Create our window
	glutCreateWindow(windowName);
	// setup the call-back functions
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutPassiveMotionFunc(passiveMotion);
	glutMotionFunc(draggedMotion);
	glutMouseFunc(mouse);
	glClearColor(1, 1, 1, 0); // set default backgorund color to white
}


#include <crtdbg.h>
int main(int argc, char * argv[]) {
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

	// initialize GLUT
	glutInit(&argc, argv);
	// initialize the application
	init("hello!");
	g_game.init();
	// keeps track of when the update method was called last
	int lastTime, now = glutGet(GLUT_ELAPSED_TIME), passed;
	int msDelay = 1000 / 50;		// the ideal game loop rate (50 FPS)
	// the glut game loop
	while (g_game.gameLoopRunning) {
		lastTime = now;				// remember the last time that update ran?
		now = glutGet(GLUT_ELAPSED_TIME);	// what time is it now?
		passed = now - lastTime;	// how much time has passed since update was called?

		glutMainLoopEvent();		// update the glut system
		update(passed);				// update the game

		glutPostRedisplay();		// request to refresh the display

		passed = glutGet(GLUT_ELAPSED_TIME) - now;	// how long did update take to process?
		int wait = msDelay-passed-1;// factor that into the next wait
		Sleep((wait > 0)?wait:0);	// wait (if it makes sense to), to throttle the game engine

		// if something caused the game to take a LONG time, ignore that time duration.
		if (wait < -1000) now = glutGet(GLUT_ELAPSED_TIME);
	}
	return 0;
}
