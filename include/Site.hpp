#ifndef SITE_HPP
#define SITE_HPP

#include <vector>
#include <iostream>

using namespace std;


class Edge;

class Site
{
public:
  Site() {}
  Site(bool inside_the_system, int x, int y, unsigned int state, vector<Edge*> edges, vector<Site*> neighbors);
  
  bool inside_the_system;
  int x;
  int y;
  unsigned int state;
  vector<Edge*> edges;
  vector<Site*> neighbors;
};

#endif

