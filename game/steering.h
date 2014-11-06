#pragma once

#include "codegiraffe/v2.h"

class Agent; // class protoype, a "forward declaration"

V2f seek(V2f target, Agent * agent);

V2f stop(Agent * agent, int a_ms);

V2f flee(V2f target, Agent * agent);

