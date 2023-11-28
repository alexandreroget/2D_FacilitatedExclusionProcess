#include "ParticleSystem.hpp"
#include <filesystem>

using namespace std;

double computeAbsorptionTime(const unsigned int L, const double rho, const SystemConfiguration system_config);
/*----------------------------------------------------------------------------------------------------------------------*
 *                                                                                                                      *
 *  Compute the absorption time of a particle system.                                                                   *
 *                                                                                                                      *
 *  1 ─ Create a particle system of size L with an uniform particle density rho.                                        *
 *  2 ─ Set t = 0                                                                                                       *
 *  3 ─ While there are active edges in the system:                                                                     *
 *      ├─ Select a random active edge and move the particle to the empty site.                                         *
 *      ├─ Set dt = 1 / number_of_active_edges                                                                          *
 *      └─ t = t + dt                                                                                                   *
 *                                                                                                                      *
 *  Return t.                                                                                                           *
 *                                                                                                                      *
 *----------------------------------------------------------------------------------------------------------------------*/

double getCriticalDensity(const unsigned int L, const SystemConfiguration system_config, const unsigned int n_jumps_max);


int main(int argc, char** argv)
{
  string save_folder = "data/";
  filesystem::create_directory(save_folder);
  
  unsigned int L[7] = {5, 10, 20, 40, 60, 80, 100};
  double rho_c = 0.329356;
  unsigned int N = 10000;
  
  SystemConfiguration system_config;
  system_config.directional_fep = 0;
  system_config.is_periodic = 1;
  system_config.vertical_boundaries = OPENED;
  system_config.horizontal_boundaries = OPENED;
  
  for(unsigned int i = 0 ; i < size(L) ; i++) {
    for(unsigned int j = 0 ; j < 2 ; j++) {
      double rho;
      string filename = save_folder + "L" + to_string(L[i]);
    
      if(j == 0) {
        rho = rho_c - 0.05;
        filename += "_subcritical.txt";
      }
      else {
        rho = rho_c + 0.005;
        filename += "_supercritical.txt";
      }

      ofstream file;
      file.open(filename,ios::out|ios::trunc);

      file<<rho<<endl;
      for(unsigned int k = 0 ; k < N ; k++) {
        double abs_time = computeAbsorptionTime(L[i], rho, system_config);
        file<<abs_time<<endl;
      }
      
      file.close();
    }
  }
  
  return 0;
}


double computeAbsorptionTime(const unsigned int L, const double rho, const SystemConfiguration system_config)
{
  ParticleSystem* system = createSystem(L, L, system_config);
  system->initSystem(rho);

  double t = 0;
  while(system->getTotalNumberOfActiveEdges() > 0) {
    double dt = system->singleStep();
    t += dt;
  }
  
  delete system;
  
  return t;
}


double getCriticalDensity(const unsigned int L, const SystemConfiguration system_config, const unsigned int n_jumps_max)
{
  ParticleSystem* system = createSystem(L, L, system_config);
  system->initSystem(0.);

  do {
    Site* site;
    do {
      unsigned int x = system->getRandomInt(0,L-1);
      unsigned int y = system->getRandomInt(0,L-1);
      site = system->getSite(x,y);
    }while(system->isAnOccupiedSite(site));
            
    system->addParticle(site);

    unsigned int i = 0;
    while((system->getTotalNumberOfActiveEdges() != 0) && (i < n_jumps_max)) {
      system->singleStep();
      i++;
    }

  }while(system->getTotalNumberOfActiveEdges() == 0);
  
  unsigned int n_particles = system->getTotalNumberOfOccupiedSites();
  double rho_c = (double) n_particles/(L*L);
  
  delete system;
  
  return rho_c;
}
