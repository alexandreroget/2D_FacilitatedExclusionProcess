#include "ParticleSystem.hpp"

using namespace std;


ParticleSystem::ParticleSystem(unsigned int horizontal_direction, unsigned int vertical_direction, unsigned int dimension, bool is_periodic, vector<unsigned int> boundaries) : 
L1(horizontal_direction), L2(vertical_direction), dimension(dimension), boundaries(boundaries), n_particles(0), n_active_edges(0), n_reservoirs(0), n(2*dimension), last_jump(nullptr)
{
  empty_site = new Site(0,L1,L2,EMPTY,vector<Edge*>({nullptr}),vector<Site*>(n,nullptr));
  empty_site->neighbors = vector<Site*>(n,empty_site);
  occupied_site = new Site(0,L1,L2,OCCUPIED,vector<Edge*>({nullptr}),vector<Site*>(n,empty_site));
  closed_site = new Site(0,L1,L2,CLOSED,vector<Edge*>({nullptr}),vector<Site*>(n,empty_site));
  
  system_size = L1*L2;
  sites = vector<Site>(system_size);
  
  n_occupied_sites = vector<unsigned int>(L1,0);
  n_active_sites = vector<unsigned int>(L1,0);
  activity = vector<unsigned int>(L1,0);
  
  setNumberOfEdges(is_periodic);
  edges = vector<Edge>(n_edges,Edge(1,vector<Site*>({nullptr,nullptr}),-1));

  random_device rd;
  rng = new mt19937(rd());
}


void ParticleSystem::initSystem(const vector<bool> particles_array)
{
  initSites(particles_array);
  initEdges();
  initReservoirs();
  
  initNumberOfOccupiedSites();
  initNumberOfActiveSites();
  initActivity();
}


void ParticleSystem::initSystem(const double particle_density)
{
  vector<bool> particles_array = randomConfiguration(particle_density);
  initSystem(particles_array);
}


void ParticleSystem::initSystem(const double density_a, const double density_b)
{
  vector<bool> particles_array = linearInterpolation(density_a,density_b);
  initSystem(particles_array);
}


void ParticleSystem::initSites(const vector<bool> particles_array)
{
  n_particles = 0;
  
  unsigned int i = 0;
  for(unsigned int y = 0 ; y < L2 ; y++) {
    for(unsigned int x = 0 ; x < L1 ; x++) {
      if(particles_array[i]) {
        n_particles++;
      }
      
      sites[i] = Site(1, x, y, particles_array[i], vector<Edge*>(n,nullptr), vector<Site*>(n,nullptr));
      setEdges(sites[i]);
      setNeighbors(sites[i]);
      
      i++;
    }
  }
}


void ParticleSystem::setNumberOfEdges(bool is_periodic)
{
  if(is_periodic) {
    if(dimension == 1) {
      n_vertical_edges = L1;
      n_horizontal_edges = 0;
    }
    else { // dimension == 2
      n_vertical_edges = (boundaries[LEFT] == OPENED) ? L1*L2 : (L1+1)*L2;
      n_horizontal_edges = (boundaries[UP] == OPENED) ? L1*L2 : L1*(L2+1);
    }
  }
  else {
    if(dimension == 1) {
      n_vertical_edges = L1+1;
      n_horizontal_edges = 0;
    }
    else {
      n_vertical_edges = (L1+1)*L2;
      n_horizontal_edges = L1*(L2+1);
    }
  }
  
  n_edges = n_vertical_edges + n_horizontal_edges;
}


void PeriodicSystem::setVerticalEdges(Site& site)
{
  unsigned int x0 = (boundaries[LEFT] == OPENED) ? L1*site.y : (L1+1)*site.y;
  unsigned int i = x0 + site.x;
  
  site.edges[LEFT] = &edges[i];
  
  unsigned int j;
  if(site.x < L1-1) {
    j = i+1;
  }
  else {
    j = (boundaries[LEFT] == OPENED) ? x0 : i+1;
  }
  
  site.edges[RIGHT] = &edges[j];
  
  /*  Sites numerotation         Edges numerotation
  
      ┌───────┬───────┐          ┌────4──┬────6──┐
      │       │       │          │       │       │
      │   0   │   1   │          0       1       0
      │       │       │          │       │       │
      ├───────┼───────┤          ├────5──┼────7──┤
      │       │       │          │       │       │
      │   2   │   3   │          2       3       2
      │       │       │          │       │       │
      └───────┴───────┘          └────4──┴────6──┘
      
      Example:  L1 = 2 ; L2 = 2
                Site 2 => x = 0 ; y = 1
                Left edge index = x + L1*y = 2
                Right edge index = [Left edge index] + 1 = 3   */
}


void NonPeriodicSystem::setVerticalEdges(Site& site)
{
  unsigned int i = site.x + (L1+1)*site.y;
  site.edges[LEFT] = &edges[i];
  site.edges[RIGHT] = &edges[i+1];
  
  /*  Sites numerotation         Edges numerotation
  
      ┌───────┬───────┐          ┌────6──┬────9──┐
      │       │       │          │       │       │
      │   0   │   1   │          0       1       2
      │       │       │          │       │       │
      ├───────┼───────┤          ├────7──┼───10──┤
      │       │       │          │       │       │
      │   2   │   3   │          3       4       5
      │       │       │          │       │       │
      └───────┴───────┘          └────8──┴───11──┘
      
      Example:  L1 = 2 ; L2 = 2
                Site 2 => x = 0 ; y = 1
                Left edge index = x + (L1+1)*y = 3
                Right edge index = [Left edge index] + 1 = 4   */
}


void ParticleSystem_1D::initEdges()
{
  initVerticalEdges();
  initActiveEdges();
}


void ParticleSystem_2D::initEdges()
{
  initVerticalEdges();
  initHorizontalEdges();
  initActiveEdges();
}


void PeriodicSystem::initVerticalEdges()
{
  unsigned int i = 0;
  for(unsigned int y = 0 ; y < L2 ; y++) { 
    for(unsigned int x = 0 ; x < L1 ; x++) {
      edges[i].is_vertical = 1;
      Site* site = getSite(x,y);
      edges[i].sites[0] = site->neighbors[LEFT];
      edges[i].sites[1] = site;
      i++;
    }
    
    if((boundaries[LEFT] == CLOSED) || (boundaries[LEFT] == RESERVOIR)) {
      edges[i].is_vertical = 1;
      Site* site = getSite(L1-1,y);
      edges[i].sites[0] = site;
      edges[i].sites[1] = site->neighbors[RIGHT];
      i++;
    }
  }
}


void NonPeriodicSystem::initVerticalEdges()
{
  unsigned int i = 0;
  for(unsigned int y = 0 ; y < L2 ; y++) {    
    for(unsigned int x = 0 ; x < L1 ; x++) {
      edges[i].is_vertical = 1;
      Site* site = getSite(x,y);
      edges[i].sites[0] = site->neighbors[LEFT];
      edges[i].sites[1] = site;
      i++;
    }
    
    {
      edges[i].is_vertical = 1;
      Site* site = getSite(L1-1,y);
      edges[i].sites[0] = site;
      edges[i].sites[1] = site->neighbors[RIGHT];
    } i++;
  }
}


void PeriodicSystem_2D::initHorizontalEdges()
{
  unsigned int i = n_vertical_edges;
  for(unsigned int x = 0 ; x < L1 ; x++) {
    for(unsigned int y = 0 ; y < L2 ; y++) {
      edges[i].is_vertical = 0;
      Site* site = getSite(x,y);
      edges[i].sites[0] = site->neighbors[UP];
      edges[i].sites[1] = site;
      i++;
    }
    
    if((boundaries[UP] == CLOSED) || (boundaries[UP] == RESERVOIR)) {
      edges[i].is_vertical = 0;
      Site* site = getSite(x,L2-1);
      edges[i].sites[0] = site;
      edges[i].sites[1] = site->neighbors[DOWN];
      i++;
    }
  }
}


void NonPeriodicSystem_2D::initHorizontalEdges()
{
  unsigned int i = n_vertical_edges;
  for(unsigned int x = 0 ; x < L1 ; x++) {
    for(unsigned int y = 0 ; y < L2 ; y++) {
      edges[i].is_vertical = 0;
      Site* site = getSite(x,y);
      edges[i].sites[0] = site->neighbors[UP];
      edges[i].sites[1] = site;
      i++;
    }
    
    {
      edges[i].is_vertical = 0;
      Site* site = getSite(x,L2-1);
      edges[i].sites[0] = site;
      edges[i].sites[1] = site->neighbors[DOWN];
    } i++;
  }
}


void ParticleSystem::initActiveEdges()
{
  active_edges.clear();
  n_active_edges = 0;

  for(unsigned int i = 0 ; i < n_edges ; i++) {
    if(isActive(&edges[i])) {
      addActiveEdge(&edges[i]);
    }
  }
}


void ParticleSystem::initReservoirs()
{
  // Vertical reservoirs
  for(unsigned int i = 0 ; i < 2 ; i++) {
    if(boundaries[i] == RESERVOIR) {
      reservoirs.insert(pair<unsigned int, Reservoir>(i,{L2,bernoulli_distribution(0.)}));
      n_reservoirs += L2;
    }
  }
  // Horizontal reservoirs
  for(unsigned int i = 2 ; i < 4 ; i++) {
    if(boundaries[i] == RESERVOIR) {
      reservoirs.insert(pair<unsigned int, Reservoir>(i,{L1,bernoulli_distribution(0.)}));
      n_reservoirs += L1;
    }
  }
}


void ParticleSystem::setReservoirDensity(unsigned int reservoir_position, double diffusion_rate)
{
  map<unsigned int, Reservoir>::iterator it = reservoirs.find(reservoir_position);
  it->second.bernoulli = bernoulli_distribution(diffusion_rate);
}


void ParticleSystem::addParticle(Site* site)
{ 
  site->state = OCCUPIED;
  n_particles++;
  
  increaseNumberOfOccupiedSites(site);
  if(getNumberOfOccupiedNeighbors(site) > 0) {
    increaseNumberOfActiveSites(site);
  }
  updateNewSiteNeighborsEdges(site,-1);
}


void ParticleSystem::removeParticle(Site* site)
{
  site->state = EMPTY;
  n_particles--;
  
  decreaseNumberOfOccupiedSites(site);
  if(getNumberOfOccupiedNeighbors(site) > 0) {
    decreaseNumberOfActiveSites(site);
  }
  updateOldSiteNeighborsEdges(site,-1);
}


void ParticleSystem::addActiveEdge(Edge* edge)
{
  active_edges.push_back(edge);
  n_active_edges++;
  increaseActivity(edge);
  
  unsigned int i = n_active_edges - 1;  // last index of active_particles
  edge->active_index = i;
}


void ParticleSystem::removeActiveEdge(Edge* edge)
{
  unsigned int i = edge->active_index;
  edge->active_index = -1;
  unsigned int j = n_active_edges - 1;  // last index of active_particles
  
  active_edges[i] = active_edges[j];
  active_edges[i]->active_index = i;
  
  active_edges.pop_back();
  n_active_edges--;
  decreaseActivity(edge);
}


double ParticleSystem::singleStep()
{
  unsigned int n_actions = n_active_edges + n_reservoirs;
  double dt = n_actions > 0 ? 1./n_actions : 0.;
  
  unsigned int i = getRandomInt(0,n_actions-1);
  
  if(i < n_active_edges) {
    moveParticle(active_edges[i]);
  }
  else {
    i -= n_active_edges;
    activateReservoir(i);
  }
  
  return dt;
}


void ParticleSystem::moveParticle(Edge* active_edge)
{
  Site* old_site;
  Site* new_site;
  if(isAnEmptySite(active_edge->sites[1])) {
    old_site = active_edge->sites[0];
    new_site = active_edge->sites[1];
  }
  else {
    old_site = active_edge->sites[1];
    new_site = active_edge->sites[0];
  }

  unsigned int new_site_direction = getNewSiteDirection(active_edge);
  
  if(isInsideTheSystem(old_site)) {
    old_site->state = EMPTY;
    
    decreaseNumberOfOccupiedSites(old_site);
    decreaseNumberOfActiveSites(old_site);
  }
  if(isInsideTheSystem(new_site)) {
    new_site->state = OCCUPIED;
    
    increaseNumberOfOccupiedSites(new_site);
    if(getNumberOfOccupiedNeighbors(new_site) > 0) {
      increaseNumberOfActiveSites(new_site);
    }
  }
  last_jump = active_edge;
  
  if(isInsideTheSystem(old_site)) {
    updateOldSiteNeighborsEdges(old_site,new_site_direction);
  }
  
  if(isInsideTheSystem(new_site)) {
    // Check if the edge is still active
    if(isActive(active_edge) == 0) {
      removeActiveEdge(active_edge);
    }

    unsigned int old_site_direction = oppositeDirection(new_site_direction);
    updateNewSiteNeighborsEdges(new_site,old_site_direction);
  }
  else {
    n_particles--;
    removeActiveEdge(active_edge);
  }
}


void ParticleSystem::activateReservoir(unsigned int i)
{
  // Select reservoir
  Reservoir* reservoir = nullptr;
  unsigned int reservoir_position;
  
  for(map<unsigned int,Reservoir>::iterator it = reservoirs.begin() ; it != reservoirs.end() ; ++it) {
    reservoir = &it->second;
    
    if(i < reservoir->size) {
      reservoir_position = it->first;
      break;
    }
    
    i -= reservoir->size;
  }

  // Select site
  unsigned int x, y;
  switch(reservoir_position) {
    case LEFT:
      x = 0 ; y = i;
      break;
    case RIGHT:
      x = L1-1 ; y = i;
      break;
    case UP:
      x = i ; y = 0;
      break;
    case DOWN:
      x = i ; y = L2-1;
      break;        
  }
  
  Site* site = getSite(x,y);
  last_jump = site->edges[reservoir_position];
  
  // Change site state
  bool old_state = isAnOccupiedSite(site);
  bool new_state = reservoir->bernoulli(*rng);
  
  if(old_state ^ new_state) {
    old_state == 0 ? addParticle(site) : removeParticle(site);
  }
}


bool ParticleSystem::isAnEmptySite(const Site* site) const
{
  return (site->state == EMPTY);
}


bool ParticleSystem::isAnOccupiedSite(const Site* site) const
{
  return (site->state == OCCUPIED);
}


bool ParticleSystem::isNotAClosedSite(const Site* site) const
{
  return (site->state != CLOSED);
}


unsigned int ParticleSystem::getNumberOfOccupiedNeighbors(const Site* site) const
{
  unsigned int n_occupied = 0;

  for(unsigned int i = 0 ; i < n ; i++) {
    if(isAnOccupiedSite(site->neighbors[i])) {
      n_occupied++;
    }
  }
  
  return n_occupied;
}


unsigned int ParticleSystem::getNewSiteDirection(Edge* active_edge) const
{
  unsigned int direction;

  if(active_edge->is_vertical) {
    direction = isAnEmptySite(active_edge->sites[0]) ? LEFT : RIGHT;
  }
  else {
    direction = isAnEmptySite(active_edge->sites[0]) ? UP : DOWN;
  }
  
  return direction;
}


unsigned int ParticleSystem::getLength() const
{
  return L1;
}


unsigned int ParticleSystem::getWidth() const
{
  return L2;
}


unsigned int ParticleSystem::getSize() const
{
  return system_size;
}


unsigned int ParticleSystem::getDimension() const
{
  return dimension;
}


unsigned int ParticleSystem::getTotalNumberOfOccupiedSites() const
{
  return n_particles;
}


unsigned int ParticleSystem::getTotalNumberOfActiveSites() const
{
  unsigned int sum_active_sites = 0;
  
  for(unsigned int x = 0 ; x < L1 ; x++) {
    sum_active_sites += n_active_sites[x];
  }
  
  return sum_active_sites;
}


unsigned int ParticleSystem::getTotalNumberOfActiveEdges() const
{
  return n_active_edges;
}


Edge* ParticleSystem::getActiveEdge(const unsigned int i) const
{
  return active_edges[i];
}


Site* ParticleSystem::getSite(const unsigned int x, const unsigned int y)
{
  unsigned int i = x + L1*y;
  return &sites[i];
}


Edge* ParticleSystem::getLastJump() const
{
  return last_jump;
}


vector<Site*> ParticleSystem::getSites(const Edge* edge) const
{
  return edge->sites;
}


vector<bool> ParticleSystem::getParticlesArray() const
{
  vector<bool> particles_array;
  
  for(unsigned int i = 0 ; i < system_size ; i++) {
    if(isAnOccupiedSite(&sites[i])) {
      particles_array.push_back(i);
    }
  }
  
  return particles_array;
}


void ParticleSystem::saveParticlesArray(const string filename) const
{
  ofstream file;
  file.open(filename,ios::out|ios::trunc);
  
  for(unsigned int i = 0 ; i < system_size ; i++) {
    if(isAnOccupiedSite(&sites[i])) {
      file<<i<<endl;
    }
  }
    
  file.close();
}


unsigned int ParticleSystem::getRandomInt(unsigned int from, unsigned int to) const
{
  uniform_int_distribution<> distribution(from,to);
  return distribution(*rng);
}


bool ParticleSystem::isInsideTheSystem(const Site* site) const
{
  return site->inside_the_system;
}


vector<bool> ParticleSystem::randomConfiguration(const double particle_density) const
{
  vector<bool> particles_array(system_size,0);
  
  unsigned int n_particles = floor(particle_density*system_size);
  for(unsigned int i = 0 ; i < n_particles ; i++) {
    particles_array[i] = 1;
  }
  
  randomPermutation(particles_array,system_size,2*system_size);
  
  return particles_array;
}


void ParticleSystem::randomPermutation(vector<bool>& tab, const unsigned int tab_size, const unsigned int n_permutations) const
{
  uniform_int_distribution<> distribution(0,tab_size-1);

  for(unsigned int i = 0 ; i < n_permutations ; i++) {
    const int a = distribution(*rng);
    const int b = distribution(*rng);
    bool p = tab[b];
    tab[b] = tab[a];
    tab[a] = p;
  }
}


unsigned int ParticleSystem::oppositeDirection(unsigned int direction) const
{
  unsigned int opposite_direction;
  switch(direction) {
    case LEFT:
      opposite_direction = RIGHT;
      break;
    case RIGHT:
      opposite_direction = LEFT;
      break;
    case UP:
      opposite_direction = DOWN;
      break;
    case DOWN:
      opposite_direction = UP;
      break;
  }
  return opposite_direction;
}

// --------------------------------------------------------------------------------------------- //

void PeriodicSystem::setNeighborOfSiteAtBorder(Site& site, const unsigned int direction)
{
  switch(boundaries[direction]) {
    case OPENED:
      switch(direction) {
        case LEFT:
          site.neighbors[direction] = getSite(L1-1,site.y);
          break;
        case RIGHT:
          site.neighbors[direction] = getSite(0,site.y);
          break;
        case UP:
          site.neighbors[direction] = getSite(site.x,L2-1);
          break;
        case DOWN:
          site.neighbors[direction] = getSite(site.x,0);
          break;
      }
      break;
    case CLOSED:
      site.neighbors[direction] = closed_site;
      break;
    case RESERVOIR:
      site.neighbors[direction] = occupied_site;
      break;
  }
}

// --------------------------------------------------------------------------------------------- //

void NonPeriodicSystem::setNeighborOfSiteAtBorder(Site& site, const unsigned int direction)
{
  switch(boundaries[direction]) {
    case OPENED:
      site.neighbors[direction] = empty_site; 
      break;
    case CLOSED:
      site.neighbors[direction] = closed_site;
      break;
    case RESERVOIR:
      site.neighbors[direction] = occupied_site;
      break;
  }
}

// --------------------------------------------------------------------------------------------- //

vector<bool> ParticleSystem_1D::linearInterpolation(const double min_density, const double max_density) const
{
  double particle_density = (min_density+max_density)/2;
  return randomConfiguration(particle_density);
}


void ParticleSystem_1D::setNeighbors(Site& site)
{
  setLeftNeighbor(site);
  setRightNeighbor(site);
}


void ParticleSystem_1D::setEdges(Site& site)
{
  setVerticalEdges(site);
}

// --------------------------------------------------------------------------------------------- //

void PeriodicSystem_1D::setLeftNeighbor(Site& site)
{
  site.neighbors[LEFT] = (site.x > 0) ?  &sites[site.x-1] : &sites[L1-1];
}


void PeriodicSystem_1D::setRightNeighbor(Site& site)
{
  site.neighbors[RIGHT] = (site.x < L1-1) ? &sites[site.x+1] : &sites[0];
}


void PeriodicSystem_1D::initActivity()
{
  for(unsigned int x = 0 ; x < L1 ; x++) {
    Site* site = getSite(x,0);
      
    Edge* right_edge = site->edges[RIGHT];
    isActive(right_edge) ? activity[x] = 1 : activity[x] = 0;
  }
}

// --------------------------------------------------------------------------------------------- //

void NonPeriodicSystem_1D::setLeftNeighbor(Site& site)
{
  if(site.x > 0) {
    site.neighbors[LEFT] = &sites[site.x-1];
  }
  else {
    setNeighborOfSiteAtBorder(site,LEFT);
  }
}


void NonPeriodicSystem_1D::setRightNeighbor(Site& site)
{
  if(site.x < L1-1) {
    site.neighbors[RIGHT] = &sites[site.x+1];
  }
  else {
    setNeighborOfSiteAtBorder(site,RIGHT);
  }
}


void NonPeriodicSystem_1D::initActivity()
{
  Site* site_0 = getSite(0,0);
      
  if(boundaries[LEFT] == OPENED) {
    Edge* left_edge_0 = site_0->edges[LEFT];
    isActive(left_edge_0) ? activity[0] = 1 : activity[0] = 0;
  }
      
  Edge* right_edge_0 = site_0->edges[RIGHT];
  if(isActive(right_edge_0)) {
    activity[0]++;
  }

  for(unsigned int x = 1 ; x < L1 ; x++) {
    Site* site = getSite(x,0);

    Edge* right_edge = site->edges[RIGHT];
    isActive(right_edge) ? activity[x] = 1 : activity[x] = 0;
  }
}

// --------------------------------------------------------------------------------------------- //

vector<bool> ParticleSystem_2D::linearInterpolation(const double min_density, const double max_density) const
{ 
  vector<bool> particles_array(L1*L2,0);

  for(unsigned int i = 0 ; i < L1 ; i++) {
    double particle_density = min_density + ((max_density - min_density) * i/(L2-1));

    unsigned int n = floor(particle_density*L2);
    
    vector<bool> column(L2,0);
    for(unsigned int j = 0 ; j < n ; j++) {
      column[j] = 1;
    }
    randomPermutation(column,L2,2*L2);
    
    for(unsigned int j = 0 ; j < L1 ; j++) {
      unsigned int k = i + j*L1;
      particles_array[k] = column[j];
    }
  }
  
  return particles_array;
}


void ParticleSystem_2D::setNeighbors(Site& site)
{
  setLeftNeighbor(site);
  setRightNeighbor(site);
  setTopNeighbor(site);
  setBottomNeighbor(site);
}


void ParticleSystem_2D::setLeftNeighbor(Site& site)
{
  if(site.x > 0) {
    unsigned int i = (site.x-1) + L1*site.y;
    site.neighbors[LEFT] = &sites[i];
  }
  else {
    setNeighborOfSiteAtBorder(site,LEFT);
  }
}


void ParticleSystem_2D::setRightNeighbor(Site& site)
{
  if(site.x < L1-1) {
    unsigned int i = (site.x+1) + L1*site.y;
    site.neighbors[RIGHT] = &sites[i];
  }
  else {
    setNeighborOfSiteAtBorder(site,RIGHT);
  }
}


void ParticleSystem_2D::setTopNeighbor(Site& site)
{
  if(site.y > 0) {
    unsigned int i = site.x + L1*(site.y-1);
    site.neighbors[UP] = &sites[i];
  }
  else {
    setNeighborOfSiteAtBorder(site,UP);
  }
}


void ParticleSystem_2D::setBottomNeighbor(Site& site)
{
  if(site.y < L2-1) {
    unsigned int i = site.x + L1*(site.y+1);
    site.neighbors[DOWN] = &sites[i];
  }
  else {
    setNeighborOfSiteAtBorder(site,DOWN);
  }
}


void ParticleSystem_2D::setEdges(Site& site)
{
  setVerticalEdges(site);
  setHorizontalEdges(site);
}


void ParticleSystem::initNumberOfOccupiedSites()
{
  for(unsigned int x = 0 ; x < L1 ; x++) {
    for(unsigned int y = 0 ; y < L2 ; y++) {
      Site* site = getSite(x,y);
      if(isAnOccupiedSite(site)) {
        n_occupied_sites[x]++;
      }
    }
  }
}


void ParticleSystem::increaseNumberOfOccupiedSites(Site* new_occupied_site)
{
  unsigned int x = new_occupied_site->x;
  n_occupied_sites[x]++;
}


void ParticleSystem::decreaseNumberOfOccupiedSites(Site* old_occupied_site)
{
  unsigned int x = old_occupied_site->x;
  n_occupied_sites[x]--;
}


unsigned int ParticleSystem::getNumberOfOccupiedSites(const unsigned int x) const
{
  return n_occupied_sites[x];
}


void ParticleSystem::initNumberOfActiveSites()
{
  for(unsigned int x = 0 ; x < L1 ; x++) {
    for(unsigned int y = 0 ; y < L2 ; y++) {
      Site* site = getSite(x,y);
      if(isAnOccupiedSite(site)) {
        if(getNumberOfOccupiedNeighbors(site) > 0) {
          n_active_sites[x]++;
        }
      }
    }
  }
}


void ParticleSystem::increaseNumberOfActiveSites(Site* new_active_site)
{
  unsigned int x = new_active_site->x;
  n_active_sites[x]++;
}


void ParticleSystem::decreaseNumberOfActiveSites(Site* old_active_site)
{
  unsigned int x = old_active_site->x;
  n_active_sites[x]--;
}


unsigned int ParticleSystem::getNumberOfActiveSites(const unsigned int x) const
{
  return n_active_sites[x];
}


void ParticleSystem::increaseActivity(Edge* new_active_edge)
{
  Site* site;
  // Vertical edge -> get Left site ; Horizontal edge -> get Bottom site
  new_active_edge->is_vertical ? site = new_active_edge->sites[0] : site = new_active_edge->sites[1];
  
  unsigned int x = site->x;
  activity[x]++;
}


void ParticleSystem::decreaseActivity(Edge* old_active_edge)
{
  Site* site;
  // Vertical edge -> get Left site ; Horizontal edge -> get Bottom site
  old_active_edge->is_vertical ? site = old_active_edge->sites[0] : site = old_active_edge->sites[1];
  
  unsigned int x = site->x;
  activity[x]--;
}


unsigned int ParticleSystem::getActivity(const unsigned int x) const
{
  return activity[x];
}

// --------------------------------------------------------------------------------------------- //

void PeriodicSystem_2D::setHorizontalEdges(Site& site)
{
  unsigned int y0 = (boundaries[UP] == OPENED) ? n_vertical_edges + L2*site.x : n_vertical_edges + (L2+1)*site.x;
  unsigned int i = y0 + site.y;
  
  site.edges[UP] = &edges[i];
  
  unsigned int j;
  if(site.y < L2-1) {
    j = i+1;
  }
  else {
    j = (boundaries[UP] == OPENED) ? y0 : i+1;
  }
  
  site.edges[DOWN] = &edges[j];

  /*  Sites numerotation         Edges numerotation
  
      ┌───────┬───────┐          ┌────4──┬────6──┐
      │       │       │          │       │       │
      │   0   │   1   │          0       1       0
      │       │       │          │       │       │
      ├───────┼───────┤          ├────5──┼────7──┤
      │       │       │          │       │       │
      │   2   │   3   │          2       3       2
      │       │       │          │       │       │
      └───────┴───────┘          └────4──┴────6──┘
      
      Example:  L1 = 2 ; L2 = 2
                Site 1 => x = 1 ; y = 0
                Top edge index = (L1+x)*L2 + y = 6
                Bottom edge index = [Top Edge index] + 1 = 7   */
}


void PeriodicSystem_2D::initActivity()
{
  for(unsigned int x = 0 ; x < L1 ; x++) {
    activity[x] = 0;
    for(unsigned int y = 0 ; y < L2 ; y++) {
      Site* site = getSite(x,y);
      
      Edge* right_edge = site->edges[RIGHT];
      if(isActive(right_edge)) {
        activity[x]++;
      }
      
      Edge* top_edge = site->edges[UP];
      if(isActive(top_edge)) {
        activity[x]++;
      }
    }
  }
}

// --------------------------------------------------------------------------------------------- //

void NonPeriodicSystem_2D::setHorizontalEdges(Site& site)
{
  unsigned int i = n_vertical_edges + (L2+1)*site.x + site.y;

  site.edges[UP] = &edges[i];
  site.edges[DOWN] = &edges[i+1];
  
  /*  Sites numerotation         Edges numerotation
  
      ┌───────┬───────┐          ┌────6──┬────9──┐
      │       │       │          │       │       │
      │   0   │   1   │          0       1       2
      │       │       │          │       │       │
      ├───────┼───────┤          ├────7──┼───10──┤
      │       │       │          │       │       │
      │   2   │   3   │          3       4       5
      │       │       │          │       │       │
      └───────┴───────┘          └────8──┴───11──┘
      
      Example:  L1 = 2 ; L2 = 2
                Site 2 => x = 0 ; y = 1
                Top edge index = (L1+1)*L2 + (L2+1)*x + y = 7
                Bottom edge index = [Top Edge index] + 1 = 8   */
}


void NonPeriodicSystem_2D::initActivity()
{
  activity[0] = 0;
  for(unsigned int y = 0 ; y < L2 ; y++) {
    Site* site = getSite(0,y);
      
    if(boundaries[LEFT] == OPENED) {
      Edge* left_edge = site->edges[LEFT];
      if(isActive(left_edge)) {
        activity[0]++;
      }
    }
      
    Edge* right_edge = site->edges[RIGHT];
    if(isActive(right_edge)) {
      activity[0]++;
    }
      
    Edge* top_edge = site->edges[UP];
    if(isActive(top_edge)) {
      activity[0]++;
    }
  }

  for(unsigned int x = 0 ; x < L1 ; x++) {
    activity[x] = 0;
    for(unsigned int y = 0 ; y < L2-1 ; y++) {
      Site* site = getSite(x,y);

      Edge* right_edge = site->edges[RIGHT];
      if(isActive(right_edge)) {
        activity[x]++;
      }
      
      Edge* top_edge = site->edges[UP];
      if(isActive(top_edge)) {
        activity[x]++;
      }
    }
    
    Site* site = getSite(x,L2-1);
    
    if(boundaries[DOWN] == OPENED) {
      Edge* bottom_edge = site->edges[DOWN];
      if(isActive(bottom_edge)) {
        activity[x]++;
      }
    }
      
    Edge* right_edge = site->edges[RIGHT];
    if(isActive(right_edge)) {
      activity[x]++;
    }
      
    Edge* top_edge = site->edges[UP];
    if(isActive(top_edge)) {
      activity[x]++;
    }
  }
}

// --------------------------------------------------------------------------------------------- //

bool FacilitatedExclusionProcess::isActive(const Edge* edge) const
{ 
  if(isNotAClosedSite(edge->sites[0]) && isNotAClosedSite(edge->sites[1])) {
    bool a, b, c, d;

    a = isAnOccupiedSite(edge->sites[0]);
    b = isAnOccupiedSite(edge->sites[1]);
    if(edge->is_vertical) {
      c = isAnOccupiedSite(edge->sites[0]->neighbors[LEFT]);
      d = isAnOccupiedSite(edge->sites[1]->neighbors[RIGHT]);
    }
    else {
      c = isAnOccupiedSite(edge->sites[0]->neighbors[UP]);
      d = isAnOccupiedSite(edge->sites[1]->neighbors[DOWN]);
    }
  
    /*   ┌───┬───┬───┬───┐
         │ c │ a │ b │ d │
         └───┴───┴───┴───┘
                         ↑
            Active edge?
       
         ┌───┬───┬───┬───┐    ┌───┬───┬───┬───┐
         │ 1 │ 1 │ 0 │ d │    │ c │ 0 │ 1 │ 1 │
         └───┴───┴───┴───┘    └───┴───┴───┴───┘
                         ↑                             ↑
               Active               Active       */
      
    return ((a && !b && c) || (!a && b && d));
  }
  else {
    return 0;
  }
}


void FacilitatedExclusionProcess::updateOldSiteNeighborsEdges(Site* old_site, int new_site_direction)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    if(i != new_site_direction) {
      Site* neighbor = old_site->neighbors[i];

      if(isAnOccupiedSite(neighbor)) {
        if(isInsideTheSystem(neighbor)) {
          isAnEmptySite(neighbor->neighbors[i]) ? removeActiveEdge(neighbor->edges[i]) : addActiveEdge(old_site->edges[i]);
      
          /*   Example: Check Right Neighbor
      
               - Case 1: isAnEmptySite(neighbor->neighbors[RIGHT]) == 1
      
                 Before jump         After jump
                ┌───┬───┬───┐      ┌───┬───┬───┐
                │ 0 │   │   │      │ 1 │   │   │
                ├───┼───┼───┤      ├───┼───┼───┤
                │ 1 │ 1 │ 0 │      │ 0 │ 1 │ 0 │
                └───┴───┴───┘      └───┴───┴───┘
                                                                ↑
                                   Edge no more active
          
               - Case 2: isAnEmptySite(neighbor->neighbors[RIGHT]) == 0
       
                 Before jump         After jump   
                ┌───┬───┬───┐      ┌───┬───┬───┐
                │ 0 │   │   │      │ 1 │   │   │
                ├───┼───┼───┤      ├───┼───┼───┤
                │ 1 │ 1 │ 1 │      │ 0 │ 1 │ 1 │
                └───┴───┴───┘      └───┴───┴───┘
                                                          ↑
                                New active edge            */
                                
                                
          if(getNumberOfOccupiedNeighbors(neighbor) == 0) {
            decreaseNumberOfActiveSites(neighbor);
          }
        }
      }
      else {
        unsigned j = oppositeDirection(i);
        if(isAnOccupiedSite(old_site->neighbors[j])) {
          removeActiveEdge(old_site->edges[i]);
          
          /*   Example: Check Right Neighbor
      
                 Before jump         After jump
                ┌───┬───┬───┐      ┌───┬───┬───┐
                │   │ 0 │   │      │   │ 1 │   │
                ├───┼───┼───┤      ├───┼───┼───┤
                │ 1 │ 1 │ 0 │      │ 1 │ 0 │ 0 │
                └───┴───┴───┘      └───┴───┴───┘
                                                                ↑
                                   Edge no more active     */
        }
      }
    }
  }
}


void FacilitatedExclusionProcess::updateNewSiteNeighborsEdges(Site* new_site, int old_site_direction)
{
  for(unsigned int i = 0 ; i < n ; i++) {
    if(i != old_site_direction) {
      Site* neighbor = new_site->neighbors[i];

      if(isAnOccupiedSite(neighbor)) {
        if(isInsideTheSystem(neighbor)) {
          isAnEmptySite(neighbor->neighbors[i]) ? addActiveEdge(neighbor->edges[i]) : removeActiveEdge(new_site->edges[i]);

          /*   Example: Check Right Neighbor
      
               - Case 1: isAnEmptySite(neighbor->neighbors[RIGHT]) == 1
      
                 Before jump         After jump
                ┌───┬───┬───┐      ┌───┬───┬───┐
                │ 1 │   │   │      │ 0 │   │   │
                ├───┼───┼───┤      ├───┼───┼───┤
                │ 0 │ 1 │ 0 │      │ 1 │ 1 │ 0 │
                └───┴───┴───┘      └───┴───┴───┘
                                                                ↑
                                     New active edge
          
               - Case 2: isAnEmptySite(neighbor->neighbors[RIGHT]) == 0
       
                 Before jump         After jump
                ┌───┬───┬───┐      ┌───┬───┬───┐
                │ 1 │   │   │      │ 0 │   │   │
                ├───┼───┼───┤      ├───┼───┼───┤
                │ 0 │ 1 │ 1 │      │ 1 │ 1 │ 1 │
                └───┴───┴───┘      └───┴───┴───┘
                                                          ↑
                                Edge no more active           */
                             
          if(getNumberOfOccupiedNeighbors(neighbor) == 1) {
            increaseNumberOfActiveSites(neighbor);
          }
        }
      }
      else {
        unsigned j = oppositeDirection(i);
        if(isAnOccupiedSite(new_site->neighbors[j])) {
          addActiveEdge(new_site->edges[i]);
          
          /*   Example: Check Right Neighbor
      
                 Before jump         After jump
                ┌───┬───┬───┐      ┌───┬───┬───┐
                │   │ 1 │   │      │   │ 0 │   │
                ├───┼───┼───┤      ├───┼───┼───┤
                │ 1 │ 0 │ 0 │      │ 1 │ 1 │ 0 │
                └───┴───┴───┘      └───┴───┴───┘
                                                                ↑
                                     New active edge         */
        }
      }
    }
  }
}

// --------------------------------------------------------------------------------------------- //

bool ConstrainedLatticeGases::isActive(const Edge* edge) const
{
  if(isNotAClosedSite(edge->sites[0]) && isNotAClosedSite(edge->sites[1])) {
    bool a = isAnOccupiedSite(edge->sites[0]);
    bool b = isAnOccupiedSite(edge->sites[1]);
    
    if(a ^ b) {
      if(a == 1) {
        if(getNumberOfOccupiedNeighbors(edge->sites[0]) > 0) {
          return 1;
        }
        else {
          return 0;
        }
      }
      else {  // (b == 1)
        if(getNumberOfOccupiedNeighbors(edge->sites[1]) > 0) {
          return 1;
        }
        else {
          return 0;
        }
      }
    }
    else {
      return 0;
    }
  }
  else {
    return 0;
  }
}


void ConstrainedLatticeGases::updateOldSiteNeighborsEdges(Site* old_site, int new_site_direction)
{
  unsigned int n_occupied_neighbors = getNumberOfOccupiedNeighbors(old_site);

  for(unsigned int i = 0 ; i < n ; i++) {
    if(i != new_site_direction) {
      Site* neighbor = old_site->neighbors[i];

      if(isAnOccupiedSite(neighbor)) {
        if(isInsideTheSystem(neighbor)) {
          if(getNumberOfOccupiedNeighbors(neighbor) > 0) {
            addActiveEdge(old_site->edges[i]);

            /*    Before jump         After jump
                 ┌───┬───┬───┐      ┌───┬───┬───┐
                 │   │ a │ 0 │      │   │ a │ 1 │
                 ├───┼───┼───┤      ├───┼───┼───┤
                 │ b │ 1 │ 1 │      │ b │ 1 │ 0 │     (a + b + c = 1)
                 ├───┼───┼───┤      ├───┼───┼───┤
                 │   │ c │   │      │   │ c │   │
                 └───┴───┴───┘      └───┴───┴───┘
                                 Add one active edge         */
          }   
          else {
            decreaseNumberOfActiveSites(neighbor);
            
            unsigned int old_site_direction = oppositeDirection(i);
            for(unsigned int j = 0 ; j < n ; j++) {
              if(isAnEmptySite(neighbor->neighbors[j]) && (j != old_site_direction)) {
                removeActiveEdge(neighbor->edges[j]);
              }
            }

            /*    Before jump         After jump
                 ┌───┬───┬───┐      ┌───┬───┬───┐
                 │   │ 0 │ 0 │      │   │ 0 │ 1 │
                 ├───┼───┼───┤      ├───┼───┼───┤
                 │ 0 │ 1 │ 1 │      │ 0 │ 1 │ 0 │
                 ├───┼───┼───┤      ├───┼───┼───┤
                 │   │ 0 │   │      │   │ 0 │   │
                 └───┴───┴───┘      └───┴───┴───┘
                               Remove three active edges     */
          }
        }
      }
      else {
        if(isNotAClosedSite(neighbor) && (n_occupied_neighbors > 0)) {
          removeActiveEdge(old_site->edges[i]);
          
          /*      Before jump         After jump
                 ┌───┬───┬───┐      ┌───┬───┬───┐
                 │   │ a │   │      │   │ a │   │
                 ├───┼───┼───┤      ├───┼───┼───┤
                 │ 0 │ 1 │ b │      │ 0 │ 0 │ b │     (a + b = 1)
                 ├───┼───┼───┤      ├───┼───┼───┤
                 │ 1 │ 0 │   │      │ 0 │ 1 │   │
                 └───┴───┴───┘      └───┴───┴───┘
                                Remove one active edge       */
        }
      }
    }
  }
}


void ConstrainedLatticeGases::updateNewSiteNeighborsEdges(Site* new_site, int old_site_direction)
{
  unsigned int n_occupied_neighbors = getNumberOfOccupiedNeighbors(new_site);

  for(unsigned int i = 0 ; i < n ; i++) {
    if(i != old_site_direction) {
      Site* neighbor = new_site->neighbors[i];

      if(isAnOccupiedSite(neighbor)) {
        if(isInsideTheSystem(neighbor)) {
          if(getNumberOfOccupiedNeighbors(neighbor) > 1) {
            removeActiveEdge(new_site->edges[i]);
        
            /*    Before jump         After jump
                 ┌───┬───┬───┐      ┌───┬───┬───┐
                 │   │ a │   │      │   │ a │   │
                 ├───┼───┼───┤      ├───┼───┼───┤
                 │ 0 │ 1 │ b │      │ 1 │ 1 │ b │     (a + b + c = 1)
                 ├───┼───┼───┤      ├───┼───┼───┤
                 │ 1 │ c │   │      │ 0 │ c │   │
                 └───┴───┴───┘      └───┴───┴───┘
                                Remove one active edge       */
          }
          else {
            increaseNumberOfActiveSites(neighbor);
          
            unsigned int new_site_direction = oppositeDirection(i);
            for(unsigned int j = 0 ; j < n ; j++) {
              if(isAnEmptySite(neighbor->neighbors[j]) && (j != new_site_direction)) {
                addActiveEdge(neighbor->edges[j]);
              }
            }
        
            /*    Before jump         After jump
                 ┌───┬───┬───┐      ┌───┬───┬───┐
                 │   │ 0 │   │      │   │ 0 │   │
                 ├───┼───┼───┤      ├───┼───┼───┤
                 │ 0 │ 1 │ 0 │      │ 1 │ 1 │ 0 │
                 ├───┼───┼───┤      ├───┼───┼───┤
                 │ 1 │ 0 │   │      │ 0 │ 0 │   │
                 └───┴───┴───┘      └───┴───┴───┘
                                Add three active edges       */
          }
        }
      }
      else {
        if(isNotAClosedSite(neighbor) && (n_occupied_neighbors > 0)) {
          addActiveEdge(new_site->edges[i]);
          
          /*    Before jump         After jump
               ┌───┬───┬───┐      ┌───┬───┬───┐
               │   │ a │   │      │   │ a │   │
               ├───┼───┼───┤      ├───┼───┼───┤
               │ 0 │ 0 │ b │      │ 0 │ 1 │ b │     (a + b = 1)
               ├───┼───┼───┤      ├───┼───┼───┤
               │ 0 │ 1 │   │      │ 1 │ 0 │   │
               └───┴───┴───┘      └───┴───┴───┘
                               Add one active edge       */
        }
      }
    }
  }
}

// --------------------------------------------------------------------------------------------- //

ParticleSystem* createSystem(const unsigned int horizontal_direction, const unsigned int vertical_direction, const SystemConfiguration system_config)
{
  ParticleSystem* system = NULL;
  
  unsigned int config = 2*system_config.directional_fep + system_config.is_periodic;
  
  switch(config) {
    case 0:
      system = new CLG_NonPeriodicSystem(horizontal_direction,vertical_direction,
                   system_config.vertical_boundaries,system_config.vertical_boundaries,
                   system_config.horizontal_boundaries,system_config.horizontal_boundaries);
      break;
    case 1:
      system = new CLG_PeriodicSystem(horizontal_direction,vertical_direction,
                   system_config.vertical_boundaries,system_config.horizontal_boundaries);
      break;
    case 2:
      system = new FEP_NonPeriodicSystem(horizontal_direction,vertical_direction,
                   system_config.vertical_boundaries,system_config.vertical_boundaries,
                   system_config.horizontal_boundaries,system_config.horizontal_boundaries);
      break;
    case 3:
      system = new FEP_PeriodicSystem(horizontal_direction,vertical_direction,
                   system_config.vertical_boundaries,system_config.horizontal_boundaries);
      break;
  }
  
  return system;
}

// --------------------------------------------------------------------------------------------- //

vector<bool> loadSystemArray(const string filename, const unsigned int system_size, const unsigned int system_dimension)
{
  unsigned int N;

  if(system_dimension == 1) {
    N = system_size;
  }
  else { // 2D
    N = system_size * system_size;
  }

  vector<bool> system_array(N,0);

  ifstream file(filename);  
  if(file.is_open()) {
    string line;  
    while(getline(file,line)) {
      unsigned int i = stoi(line);
      system_array[i] = 1;
    }
    file.close();
  }
  else {
    cerr<<"Unable to open file."<<endl;
  }
  
  return system_array;
}

// --------------------------------------------------------------------------------------------- //
