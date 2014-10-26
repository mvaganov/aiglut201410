aiglut201410
============

a C / C++ game using freeGLUT. Intended for student use, using Microsoft Visual Studio (Express) 2012 or higher.

lesson plan TODO:
//---wk1day1 --- basic display/input/update loop, display/update list
// draw colored lines and circle
// mouse and key input
// draw player
// draw ai
// select with mouse
// set velocity with key presses
	// turn and strafe with key presses
// game loop

//---wk1day2 --- basic steering behavior
// player moves (seek)
/*
vec2	seek(vec2	target,	Agent	agent)	{
				vec2	desired	=	target	-	agent.position;
				desired	*=	agent.maxSpeed	/	desired.mag();
				vec2	force	=	desired	-	agent.velocity;
				return	force	*	(agent.maxForce	/	agent.maxSpeed);
}
*/
// player moves (arrive)
// ai moves (arrive)
// ai random-walks
// multiple ai
// ai seeks/flees
	// pursue
	// evade
	// offset pursuit

//---wk2day3 --- fuzzy state machines for multiple steering impulses
// player attacks with dangerous projectiles
	// context map?
// player deals damage to ai
// ai gets killed
// ai respawns
// ai flees projectiles
// ai attacks when in range
// player gets killed
// player respawns

//---wk2day4 --- environment aware steering behaviors (more fuzzy logic)
// obstacles
// obstacle collision detection
// ai avoid obstacles
	// ai avoid other ai agents
	// ai try to keep some sort of comfort zone
// ai flee behind cover when low on HP

//---wk3day5 --- finite state machines
// ai follows you when you get in range
// ai flees when low on hp
// ai regenerates HP
// ai heal each other when close
	// loneliness meter causes ai to seek out friendlies, which reduces loneliness
	// aggressive meter causes ai to attack if player (or other AI) in range
	// hungriness meter causes ai to seek out food resources

//---wk3day6 --- graph navigation
// create a graph
// draw the nav graph
// ai move using the nav graph (A*)

//---wk4day7 --- group behavior
// ai separation
/*
vec2	Separation(Agent	agent,	List<Agent>	neighbors)	{
		if	(neighbors.empty())	return	vec2.zero;
		vec2	totalForce	=	vec2.zero;
		for	(Agent	neighbor	:	neighbors)	{
					vec2	pushForce	=	(agent.pos	-	neighbor.pos)
					totalForce	+=	1-	(pushForce	/	agent.neighborRadius);
		}
		
		totalForce	/=	neighbors.count();
		totalForce	*=	agent.maxForce;
}
*/
// ai cohesion
/*
vec2	Cohesion(Agent	agent,	List<Agent>	neighbors)	{
		if	(neighbors.empty())	return	vec2.zero;
		vec2	centerOfMass	=	agent.position;
		for	(Agent	neighbor	:	neighbors)
					centerOfMass	+=	neighbor.position;
		centerOfMass	/=	neighbors.count();
		vec2	desired	=	centerOfMass	-	agent.position;
		desired	*=	agent.maxSpeed	/	desired.mag();
		vec2	force	=	desired	-	agent.velocity;
		return	force	*	(agent.maxForce	/	agent.maxSpeed);
}
*/
// ai alignment
/*
vec2	Alignment(Agent	agent,	List<Agent>	neighbors)	{
		if	(neighbors.empty())	return	vec2.zero;
		vec2	avgHeading	=	norm(	agent.velocity	);
		for	(Agent	neighbor	:	neighbors)
				avgHeading	+=	norm(	neighbor.velocity	);
		avgHeading	/=	neighbors.count();
		vec2	desired	=	avgHeading	*	agent.maxSpeed;
		vec2	force	=	desired	-	agent.velocity;
		return	force	*	(agent.maxForce	/	agent.maxSpeed);
}
*/
// ai flocking

//---wk4day8 --- midterm review

//---wk5day9 --- midterm

//---wk6day10 --- group project
	//group project
	//markov chain

//---wk6day11 --- behavior trees
//---wk7day12 --- genetic algorithm, gini impurity, neural network
//---wk7day13 --- project presentations
