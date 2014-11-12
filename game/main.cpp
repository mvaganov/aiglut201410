// properties->"VC++ Directories"->"Include Directories"->add "$(ProjectDir)\freeglut\include"
// properties->"VC++ Directories"->"Library Directories"->add "$(ProjectDir)\freeglut\lib"
// properties->"Debugging"->"Environment"->set to "PATH=%PATH%;$(ProjectDir)\freeglut\bin"

#include <stdio.h>
#include <GL/freeglut.h>
#include "codegiraffe/glutrenderingcontext.h"
#include "game.h"
#include "bullet.h"

/** the screen pixels, cartesian plane min, cartesian plane max */
GLUTRenderingContext g_screen(V2f(600, 600), V2f(-3.0, -3.0), V2f(8.0, 8.0));
/** whether or not to keep running the game loop */
bool g_gameLoopRunning = true;
bool g_paused = false;

Game g_game;

/**
 * @param key what keyboard key was pressed (ASCII codes)
 * @param x/y where the mouse was at the time
 * @note: the following may be helpful:<code>
 int state=glutGetModifiers();
 if (state & GLUT_ACTIVE_SHIFT)	// shift is pressed</code>
 */
void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27:	g_gameLoopRunning = false;  break;
	case ' ':
		{
			Agent * a = g_game.selected;
			if(a != NULL) {
				Agent * bullet = new Bullet(CircF(a->body.center, 0.2f), &g_game,
					a->direction, 2);
				bullet->parent = a;
				g_game.agents.add(bullet);
			}
		}
		break;
	case 'p':
		g_paused = !g_paused;
		break;
	}
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
	V2f click = g_screen.convertPixelsToCartesian(V2f((float)x, (float)y));
	switch (button) {
	case 3: g_screen.zoom(1/1.03125f, click); break;
	case 4: g_screen.zoom(1.03125f, click); break;
	case GLUT_MIDDLE_BUTTON:
		switch (state) {
		case GLUT_DOWN: g_screen.scrollStart(click); break;
		case GLUT_UP:   g_screen.scrollStop(click);  break;
		}
		break;
	case GLUT_LEFT_BUTTON:
		if(state == GLUT_DOWN) {
			Agent * a = g_game.getAgentAt(click);
			if(a != NULL) {
				if(g_game.selected != NULL) {
					g_game.selected->playerControlled = false;
					g_game.selected->behavior = Agent::BEHAVIOR_AGGRO;
				}
				g_game.selected = a;
				a->behavior = Agent::BEHAVIOR_SEEK;
				a->playerControlled = true;
			} else {
				g_game.mouseClick = click;
				if(g_game.selected != NULL) {
					a = g_game.selected;
					V2f delta = click - a->body.center;
					delta.normalize();
					a->direction = delta; // normalized move vector
					a->targetPosition = click;
				}
			}
		}
		break;
	case GLUT_RIGHT_BUTTON:
		if (state == GLUT_DOWN) {
			Agent * a = g_game.getAgentAt(click);
			if (a != NULL) {
				a->showDebugLines = !a->showDebugLines;
			}
		}
		break;
	}
	glutPostRedisplay();
}

/** @param x/y the coordinate of where the mouse is */
void passiveMotion(int x, int y) {
	g_game.mousePosition = g_screen.convertPixelsToCartesian(V2f((float)x, (float)y));
}

/** @param x/y the coordinate of where the mouse is being dragged */
void draggedMotion(int x, int y) {
	V2f drag = g_screen.convertPixelsToCartesian(V2f((float)x, (float)y));
	g_screen.scrollDrag(drag);
	g_game.mouseDragged = drag;
	glutPostRedisplay();
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

	g_game.display(g_screen);

	glFlush();						// print everything to the (off?) screen buffer
	glutSwapBuffers();				// swap the draw buffers
}

void init(const char * windowName) {
	// setup the screen using the graphics context
	g_screen.glSetupScreen();
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
	// keeps track of when the update method was called last
	int lastTime, now = glutGet(GLUT_ELAPSED_TIME), passed;
	int msDelay = 1000 / 50;		// the ideal game loop rate (50 FPS)
	// the glut game loop
	while (g_gameLoopRunning) {
		lastTime = now;				// remember the last time that update ran?
		now = glutGet(GLUT_ELAPSED_TIME);	// what time is it now?
		passed = now - lastTime;	// how much time has passed since update was called?

		glutMainLoopEvent();		// update the glut system
		if(!g_paused)
			update(passed);				// update the game

		glutPostRedisplay();		// request to refresh the display

		passed = glutGet(GLUT_ELAPSED_TIME) - now;	// how long did update take to process?
		int wait = msDelay-passed-1;// factor that into the next wait
		Sleep((wait > 0)?wait:0);	// wait (if it makes sense to), to throttle the game engine
	}
	return 0;
}
