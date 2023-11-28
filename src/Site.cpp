#include "Site.hpp"
#include "Edge.hpp"


Site::Site(bool inside_the_system, int x, int y, unsigned int state, vector<Edge*> edges, vector<Site*> neighbors) : 
inside_the_system(inside_the_system), x(x), y(y), state(state), edges(edges), neighbors(neighbors) {}
