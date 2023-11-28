#include "Edge.hpp"
#include "Site.hpp"


Edge::Edge(bool vertical, vector<Site*> sites, int active_index) : 
is_vertical(is_vertical), sites(sites), active_index(active_index) {}
