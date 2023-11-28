#include "ParticleSystem.hpp"
#include <filesystem>

using namespace std;

void computeActiveSitesDensity(const unsigned int L, const double rho_c, const string save_folder, const unsigned int n_jumps);
/*----------------------------------------------------------------------------------------------------------------------*
 *                                                                                                                      *
 *  Compute rho, rho_active and the activity as a function of x.                                                        *
 *                                                                                                                      *
 *  1 ─ Create a closed particle system of size L with two reservoirs.                                                  *
 *  2 ─ Set the particle density to the critical value rho_c.                                                           *
 *  3 ─ Set the left reservoir density at 0. and the right reservoir density at 1.                                      *
 *                                                                                                                      *
 *  4 ─ Set t = 0                                                                                                       *
 *  5 ─ While t < L²:                                                                                                   *
 *      ├─ Select a random active edge and move the particle to the empty site OR Activate a reservoir.                 *
 *      ├─ Set dt = 1 / number_of_actions                                                                               *
 *      └─ t = t + dt                                                                                                   *
 *                                                                                                                      *
 *  6 ─ For i = 0 to n_jumps:                                                                                           *
 *      └─ For x = 0 to L-1:                                                                                            *
 *         └─ Compute the number of occupied sites, active sites, and active edges at section x.                        *
 *                                                                                                                      * 
 *  7 ─ For x = 0 to L-1:                                                                                               *
 *      ├─ Compute rho(x) = [number of occupied sites at section x] / (L * n_jumps)                                     * 
 *      ├─ Compute rho_active(x) = [number of active sites at section x] / (L * n_jumps)                                *
 *      ├─ Compute activity(x) = [number of active edges at section x] / (L * n_jumps)                                  *
 *      └─ Write x, rho(x), rho_active(x) and the activity(x) in the results.txt file.                                  *
 *                                                                                                                      *
 *----------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
  string save_folder = "data/";
  filesystem::create_directory(save_folder);
    
  unsigned int L = 300; 
  double rho_c = 0.329356;
  unsigned int n_jumps = 20000000;

  computeActiveSitesDensity(L, rho_c, save_folder, n_jumps);

  return 0;
}


void computeActiveSitesDensity(const unsigned int L, const double rho_c, const string save_folder, const unsigned int n_jumps)
{
  SystemConfiguration system_config;
  system_config.directional_fep = 0;
  system_config.is_periodic = 0;
  system_config.vertical_boundaries = RESERVOIR;
  system_config.horizontal_boundaries = CLOSED;

  unsigned int L2 = L*L;
  
  ParticleSystem* system = NULL;
  double t;

  do {
    system = createSystem(L, L, system_config);
    system->initSystem(rho_c);
    system->setReservoirDensity(LEFT, 0.);
    system->setReservoirDensity(RIGHT, 1.);
  
    t = 0.;
    do {
      double dt = system->singleStep();
      t += dt;
    }while((t < L2) && (system->getTotalNumberOfActiveEdges() > 0));
    
  }while(t < L2);
  
  string filename = save_folder + "system.txt";
  system->saveParticlesArray(filename);
  
  vector<unsigned int> sum_occupied_sites(L,0.);
  vector<unsigned int> sum_active_sites(L,0.);
  vector<unsigned int> sum_activity(L,0.);
  
  for(unsigned int i = 0 ; i < n_jumps ; i++) {
    system->singleStep();
    
    for(unsigned x = 0 ; x < L ; x++) {
      unsigned int n_occupied_sites = system->getNumberOfOccupiedSites(x);
      sum_occupied_sites[x] += n_occupied_sites;
      
      unsigned int n_active_sites = system->getNumberOfActiveSites(x);
      sum_active_sites[x] += n_active_sites;
      
      unsigned int activity = system->getActivity(x);
      sum_activity[x] += activity;
    }
  }
  
  filename = save_folder + "results.txt";
  ofstream file;
  file.open(filename,ios::out|ios::trunc);
   
  for(unsigned int x = 0 ; x < L ; x++) {
    double rho = (double) (sum_occupied_sites[x] / L) / n_jumps;
    double rho_active = (double) (sum_active_sites[x] / L) / n_jumps;
    double A = (double) (sum_activity[x] / L) / n_jumps;
    file<<rho<<" ; "<<rho_active<<" ; "<<A<<endl;
  }
  
  file.close();
  
  delete system;
}
