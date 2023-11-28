#include "ParticleSystem.hpp"
#include <filesystem>

using namespace std;

double getParticleDensity(const unsigned int L, const double lambda, const double t_max);
/*----------------------------------------------------------------------------------------------------------------------*
 *                                                                                                                      *
 *  Compute the particle density rho(lambda).                                                                           *
 *                                                                                                                      *
 *  1 ─ Create a closed particle system of size L with two reservoirs.                                                  *
 *  2 ─ Set the particle density at 0.                                                                                  *
 *  3 ─ Set the two reservoirs densities at lambda.                                                                     *
 *                                                                                                                      *
 *  4 ─ Set t = 0                                                                                                       *
 *  5 ─ While t < t_max:                                                                                                *
 *      ├─ Select a random active edge and move the particle to the empty site OR Activate a reservoir.                 *
 *      ├─ Set dt = 1 / number_of_actions                                                                               *
 *      └─ t = t + dt                                                                                                   *
 *                                                                                                                      *
 *  6 ─ Return the particle density.                                                                                    *
 *                                                                                                                      *
 *----------------------------------------------------------------------------------------------------------------------*/


void computeDiffusionCurrent(const unsigned int L, double lambda, double epsilon, const double t_max, const string filename, const unsigned int n_records);
/*----------------------------------------------------------------------------------------------------------------------*
 *                                                                                                                      *
 *  Compute the diffusion current.                                                                                      *
 *                                                                                                                      *
 *  1 ─ Create a periodic particle system of size L with two reservoirs.                                                *
 *  2 ─ Set the particle density at rho(lambda).                                                                        *
 *  3 ─ Set the left reservoir density at lambda.                                                                       *
 *  4 ─ Set the right reservoir density at lambda + epsilon.                                                            *
 *                                                                                                                      *
 *  5 ─ Set t = 0                                                                                                       *
 *  6 ─ While t < t_max:                                                                                                *
 *      ├─ Select a random active edge and move the particle to the empty site OR Activate a reservoir.                 *
 *      ├─ Set dt = 1 / number_of_actions                                                                               *
 *      ├─ t = t + dt                                                                                                   *
 *      └─ Update the diffusion current value.                                                                          *
 *                                                                                                                      * 
 *  7 ─ Write the diffusion current values in the filename.txt file.                                                    *
 *                                                                                                                      *
 *----------------------------------------------------------------------------------------------------------------------*/


int main(int argc, char** argv)
{
  const string save_folder = "data/";
  filesystem::create_directory(save_folder);
  
  unsigned int L = 100;
  double lambda[4] = {0.1, 0.25, 0.5, 0.75};
  double epsilon = 0.01;
  double t_max = 100*(L*L);  
    
  #pragma omp parallel for shared(lambda, save_folder)
  for(unsigned int i = 0 ; i < 4 ; i++) {
    string filename = save_folder + "diffusion_current_lambda0" + to_string((int) (lambda[i]*100)) + ".txt";    
    computeDiffusionCurrent(L, lambda[i], epsilon, t_max, filename, 400000);
  }
  
  return 0;
}


double getParticleDensity(const unsigned int L, const double lambda, const double t_max)
{
  SystemConfiguration system_config;
  system_config.directional_fep = 0;
  system_config.is_periodic = 0;
  system_config.vertical_boundaries = CLOSED;
  system_config.horizontal_boundaries = RESERVOIR;
  
  ParticleSystem* system = createSystem(L, L, system_config);
  system->initSystem(0.);
  system->setReservoirDensity(LEFT, lambda);
  system->setReservoirDensity(RIGHT, lambda);
  
  double t = 0.;
  while(t < t_max) {
    t += system->singleStep();
  }
  
  double rho = system->getTotalNumberOfOccupiedSites() / (L*L);

  delete system;
  
  return rho;
}


void computeDiffusionCurrent(const unsigned int L, double lambda, double epsilon, const double t_max, const string filename, const unsigned int n_records)
{
  SystemConfiguration system_config = {0, 1, RESERVOIR, OPENED};
  ParticleSystem* system = createSystem(L, L, system_config);
  
  double rho = getParticleDensity(L, lambda, L*L);
  system->initSystem(rho, rho+epsilon);
  
  system->setReservoirDensity(LEFT, lambda);
  system->setReservoirDensity(RIGHT, lambda+epsilon);

  ofstream file;
  file.open(filename,ios::out|ios::trunc);
  
  unsigned int x_diffusion = system->getLength()/2;
  
  double t = 0.;
  int N = 0;
  file<<t<<" ; "<<N<<endl;

  double dt = t_max/n_records;
  double t_next_record = dt;
  
  do {
    t += system->singleStep();

    Edge* jump = system->getLastJump();
    vector<Site*> sites = system->getSites(jump);
    if(jump->is_vertical && (sites[0]->x == x_diffusion)) {
      system->isAnOccupiedSite(sites[0]) ? N++ : N--;
    }
    
    if(t >= t_next_record) {
      file<<t<<" ; "<<N<<endl;
      t_next_record += dt;
    }
  }while(t < t_max);

  file.close();
  
  delete system;
}
