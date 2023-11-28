#ifndef EDGE_HPP
#define EDGE_HPP

#include <vector>

using namespace std;


class Site;

class Edge
{
public:
  Edge() {}
  Edge(bool is_vertical, vector<Site*> sites, int active_index);
  
  bool is_vertical;
  vector<Site*> sites;
  int active_index; // Index in the active_edges array (equal to -1 if the edge is not active)
};

#endif
