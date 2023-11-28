#ifndef PARTICLE_SYSTEM_HPP
#define PARTICLE_SYSTEM_HPP

#include <omp.h>

#include <math.h>
#include <random>

#include <iostream>
#include <fstream>

#include <map>

#include "Site.hpp"
#include "Edge.hpp"

#define EMPTY 0
#define OCCUPIED 1

#define OPENED 1
#define CLOSED 2
#define RESERVOIR 3

#define LEFT 0
#define RIGHT 1
#define UP 2
#define DOWN 3

using namespace std;


typedef struct {
  bool directional_fep = 1; // 0 = CLG ; 1 = FEP
  bool is_periodic = 1;
  unsigned int horizontal_boundaries = 1; // Right and Left boundaries
  unsigned int vertical_boundaries = 1; // Top and Bottom boundaries
} SystemConfiguration;

typedef struct {
  unsigned int size;
  bernoulli_distribution bernoulli;
} Reservoir;


class ParticleSystem
{
public:
  ParticleSystem() {}
  ParticleSystem(unsigned int horizontal_direction, unsigned int vertical_direction, unsigned int dimension, bool is_periodic, vector<unsigned int> boundaries);
  virtual ~ParticleSystem() {}
  
  void initSystem(const vector<bool> particles_array);
  void initSystem(const double particle_density);
  void initSystem(const double density_a, const double density_b);
  
  void setReservoirDensity(unsigned int reservoir_position, double diffusion_rate);
  
  void addParticle(Site* site);
  void removeParticle(Site* site);
  
  double singleStep();
  void moveParticle(Edge* active_edge);
  
  unsigned int getLength() const;
  unsigned int getWidth() const;
  unsigned int getSize() const;
  unsigned int getDimension() const;
  unsigned int getTotalNumberOfOccupiedSites() const;
  unsigned int getTotalNumberOfActiveSites() const;
  unsigned int getTotalNumberOfActiveEdges() const;
  
  Edge* getActiveEdge(const unsigned int i) const;
  Site* getSite(const unsigned int x, const unsigned int y);
  Edge* getLastJump() const;
  vector<Site*> getSites(const Edge* edge) const;
  
  unsigned int getNumberOfOccupiedSites(const unsigned int x) const;
  unsigned int getNumberOfActiveSites(const unsigned int x) const;
  unsigned int getActivity(const unsigned int x) const;
  
  bool isAnEmptySite(const Site* site) const;
  bool isAnOccupiedSite(const Site* site) const;
  bool isNotAClosedSite(const Site* site) const;
  bool isInsideTheSystem(const Site* site) const;
  
  vector<bool> getParticlesArray() const;
  void saveParticlesArray(const string filename) const;
    
  unsigned int getRandomInt(const unsigned int from, const unsigned int to) const;
  
protected:
  void initSites(const vector<bool> particles_array);
  void setNumberOfEdges(bool is_periodic);
  virtual void setEdges(Site& site) = 0;
  virtual void setVerticalEdges(Site& site) = 0;
  virtual void setNeighbors(Site& site) = 0;
  virtual void setNeighborOfSiteAtBorder(Site& site, const unsigned int direction) = 0;

  virtual void initEdges() = 0;
  virtual void initVerticalEdges() = 0;
  void initActiveEdges();
  
  void initReservoirs();

  void initNumberOfOccupiedSites();
  void increaseNumberOfOccupiedSites(Site* new_occupied_site);
  void decreaseNumberOfOccupiedSites(Site* old_occupied_site);
    
  void initNumberOfActiveSites();
  void increaseNumberOfActiveSites(Site* new_active_site);
  void decreaseNumberOfActiveSites(Site* old_active_site);

  virtual void initActivity() = 0;
  void increaseActivity(Edge* new_active_edge);
  void decreaseActivity(Edge* old_active_edge);

  vector<bool> randomConfiguration(const double particle_density) const;
  void randomPermutation(vector<bool>& tab, const unsigned int tab_size, const unsigned int n_permutations) const;
  virtual vector<bool> linearInterpolation(const double min_density, const double max_density) const = 0;

  void addActiveEdge(Edge* edge);
  void removeActiveEdge(Edge* edge);
  
  void activateReservoir(unsigned int i);
  
  unsigned int getNumberOfOccupiedNeighbors(const Site* site) const;
  unsigned int getNewSiteDirection(Edge* active_edge) const;
      
  virtual bool isActive(const Edge* edge) const = 0;
  
  virtual void updateOldSiteNeighborsEdges(Site* old_site, int new_site_direction) = 0;
  virtual void updateNewSiteNeighborsEdges(Site* new_site, int old_site_direction) = 0;
  
  unsigned int oppositeDirection(unsigned int direction) const;
  
protected:
  unsigned int L1; 		// horizontal direction
  unsigned int L2; 		// vertical direction
  unsigned int system_size; 	// = L1*L2
  
  unsigned int dimension;	// = 1 or 2
  vector<unsigned int> boundaries;
  
  vector<Site> sites;
  vector<Edge> edges;
  vector<Edge*> active_edges;
  
  map<unsigned int, Reservoir> reservoirs;
  
  unsigned int n_particles;
  unsigned int n_edges;
  unsigned int n_vertical_edges;
  unsigned int n_horizontal_edges;
  unsigned int n_active_edges;
  unsigned int n_reservoirs;
  unsigned int n;
  
  vector<unsigned int> n_occupied_sites;
  vector<unsigned int> n_active_sites;
  vector<unsigned int> activity;
  
  Edge* last_jump;
  
  Site* empty_site;
  Site* occupied_site;
  Site* closed_site;
  
  mt19937* rng;
};


class PeriodicSystem : public virtual ParticleSystem
{
public:
  PeriodicSystem() {}
  virtual ~PeriodicSystem() {}
  
protected:
  void setVerticalEdges(Site& site);
  void initVerticalEdges();
  
  void setNeighborOfSiteAtBorder(Site& site, const unsigned int direction);
};


class NonPeriodicSystem : public virtual ParticleSystem
{
public:
  NonPeriodicSystem() {}
  virtual ~NonPeriodicSystem() {}
  
protected:
  void setVerticalEdges(Site& site);
  void initVerticalEdges();
  
  void setNeighborOfSiteAtBorder(Site& site, const unsigned int direction);
};


class ParticleSystem_1D : public virtual ParticleSystem
{
public:
  ParticleSystem_1D() {}
  virtual ~ParticleSystem_1D() {}

protected:
  vector<bool> linearInterpolation(const double min_density, const double max_density) const;

  void setNeighbors(Site& site);
  void setEdges(Site& site);
  
  void initEdges();

  virtual void setLeftNeighbor(Site& site) = 0;
  virtual void setRightNeighbor(Site& site) = 0;
};



class PeriodicSystem_1D : public ParticleSystem_1D, public PeriodicSystem
{
public:
  PeriodicSystem_1D() {}
  virtual ~PeriodicSystem_1D() {}
  
protected:
  void setLeftNeighbor(Site& site);
  void setRightNeighbor(Site& site);
  
  void initActivity();
};


class NonPeriodicSystem_1D : public ParticleSystem_1D, public NonPeriodicSystem
{
public:
  NonPeriodicSystem_1D() {}
  virtual ~NonPeriodicSystem_1D() {}

protected:
  void setLeftNeighbor(Site& site);
  void setRightNeighbor(Site& site);
  
  void initActivity();
};


class ParticleSystem_2D : public virtual ParticleSystem
{
public:
  ParticleSystem_2D() {}
  virtual ~ParticleSystem_2D() {}
  
protected:
  vector<bool> linearInterpolation(const double min_density, const double max_density) const;

  void setNeighbors(Site& site);
  void setLeftNeighbor(Site& site);
  void setRightNeighbor(Site& site);
  void setTopNeighbor(Site& site);
  void setBottomNeighbor(Site& site);
  
  void setEdges(Site& site);
  virtual void setHorizontalEdges(Site& site) = 0;
  
  void initEdges();
  virtual void initHorizontalEdges() = 0;
};


class PeriodicSystem_2D : public ParticleSystem_2D, public PeriodicSystem
{
public:
  PeriodicSystem_2D() {}
  virtual ~PeriodicSystem_2D() {}
  
protected:
  void setHorizontalEdges(Site& site);
  void initHorizontalEdges();
  
  void initActivity();
};


class NonPeriodicSystem_2D : public ParticleSystem_2D, public NonPeriodicSystem
{
public:
  NonPeriodicSystem_2D() {}
  virtual ~NonPeriodicSystem_2D() {}
  
protected:
  void setHorizontalEdges(Site& site);
  void initHorizontalEdges();
  
  void initActivity();
};


class FacilitatedExclusionProcess : public virtual ParticleSystem
{
public:
  FacilitatedExclusionProcess() {}
  virtual ~FacilitatedExclusionProcess() {}
  
private:
  bool isActive(const Edge* edge) const;
  
  void updateOldSiteNeighborsEdges(Site* old_site, int new_site_direction);
  void updateNewSiteNeighborsEdges(Site* new_site, int old_site_direction);
};


class ConstrainedLatticeGases : public virtual ParticleSystem
{
public:
  ConstrainedLatticeGases() {}
  virtual ~ConstrainedLatticeGases() {}
  
private:
  bool isActive(const Edge* edge) const;

  void updateOldSiteNeighborsEdges(Site* old_site, int new_site_direction);
  void updateNewSiteNeighborsEdges(Site* new_site, int old_site_direction);
};


class SEP_PeriodicSystem : public PeriodicSystem_1D, public FacilitatedExclusionProcess
{
public:
  SEP_PeriodicSystem(unsigned int horizontal_direction) : ParticleSystem(horizontal_direction,1,1,1,vector<unsigned int>({OPENED,OPENED})), PeriodicSystem_1D(), FacilitatedExclusionProcess() {}
  virtual ~SEP_PeriodicSystem() {}
};


class SEP_NonPeriodicSystem : public NonPeriodicSystem_1D, public FacilitatedExclusionProcess
{
public:
  SEP_NonPeriodicSystem(unsigned int horizontal_direction, unsigned int left_boundary, unsigned int right_boundary) : ParticleSystem(horizontal_direction,1,1,0,vector<unsigned int>({left_boundary, right_boundary})), NonPeriodicSystem_1D(), FacilitatedExclusionProcess() {}
  virtual ~SEP_NonPeriodicSystem() {}
};


class FEP_PeriodicSystem : public PeriodicSystem_2D, public FacilitatedExclusionProcess
{
public:
  FEP_PeriodicSystem(unsigned int horizontal_direction, unsigned int vertical_direction, unsigned int vertical_boundaries, unsigned int horizontal_boundaries) : ParticleSystem(horizontal_direction,vertical_direction,2,1,vector<unsigned int>({vertical_boundaries,vertical_boundaries,horizontal_boundaries,horizontal_boundaries})), PeriodicSystem_2D(), FacilitatedExclusionProcess() {}
  virtual ~FEP_PeriodicSystem() {}
};


class FEP_NonPeriodicSystem : public NonPeriodicSystem_2D, public FacilitatedExclusionProcess
{
public:
  FEP_NonPeriodicSystem(unsigned int horizontal_direction, unsigned int vertical_direction, unsigned int left_boundary, unsigned int right_boundary, unsigned int top_boundary, unsigned int bottom_boundary) : ParticleSystem(horizontal_direction,vertical_direction,2,0,vector<unsigned int>({left_boundary,right_boundary,top_boundary,bottom_boundary})), NonPeriodicSystem_2D(), FacilitatedExclusionProcess() {}
  virtual ~FEP_NonPeriodicSystem() {}
};


class CLG_PeriodicSystem : public PeriodicSystem_2D, public ConstrainedLatticeGases
{
public:
  CLG_PeriodicSystem(unsigned int horizontal_direction, unsigned int vertical_direction, unsigned int  vertical_boundaries, unsigned int horizontal_boundaries) : ParticleSystem(horizontal_direction,vertical_direction,2,1,vector<unsigned int>({vertical_boundaries,vertical_boundaries,horizontal_boundaries,horizontal_boundaries})), PeriodicSystem_2D(), ConstrainedLatticeGases() {}
  virtual ~CLG_PeriodicSystem() {}
};


class CLG_NonPeriodicSystem : public NonPeriodicSystem_2D, public ConstrainedLatticeGases
{
public:
  CLG_NonPeriodicSystem(unsigned int horizontal_direction, unsigned int vertical_direction, unsigned int left_boundary, unsigned int right_boundary, unsigned int top_boundary, unsigned int bottom_boundary) : ParticleSystem(horizontal_direction,vertical_direction,2,0,vector<unsigned int>({left_boundary,right_boundary,top_boundary,bottom_boundary})), NonPeriodicSystem_2D(), ConstrainedLatticeGases() {}
  virtual ~CLG_NonPeriodicSystem() {}
};


ParticleSystem* createSystem(const unsigned int horizontal_direction, const unsigned int vertical_direction, const SystemConfiguration system_config);
                             
vector<bool> loadSystemArray(const string filename, const unsigned int system_size, const unsigned int system_dimension);

#endif
