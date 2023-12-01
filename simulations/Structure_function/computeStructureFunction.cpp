#include "ParticleSystem.hpp"
#include <filesystem>

#include <complex>
#include <fftw3.h>

// #include <omp.h>

using namespace std;

void computeStructureFunction(const unsigned int L, const double rho, const unsigned int n_samples, const string filename);
/*----------------------------------------------------------------------------------------------------------------------*
 *                                                                                                                      *
 *  Compute the structure function S_rho(k).                                                                            *
 *                                                                                                                      *
 *  1 ─ For i = 0 to n_samples:                                                                                         *
 *      ├─ Create a periodic particle system of size L with an uniform particle density rho.                            *
 *      ├─ Set t = 0                                                                                                    *
 *      ├─ While t < L²:                                                                                                *
 *      │  ├─ Select a random active edge and move the particle to the empty site.                                      *
 *      │  ├─ Set dt = 1 / number_of_active_edges                                                                       *
 *      │  └─ t = t + dt                                                                                                *
 *      └─ Compute S^i_rho(k) = |FFT(eta(x) - rho)|²                                                                    *
 *  2 ─ Compute S_rho(k) = sum_{i = 0 to n_samples-1} S^i_rho(k) / n_samples                                            *
 *  3 ─ Write S_rho(k) values in the [filename] file.                                                                   *
 *                                                                                                                      *
 *----------------------------------------------------------------------------------------------------------------------*/

int main()
{
  string save_folder = "data/";
  filesystem::create_directory(save_folder);

  unsigned int L = 100;
  double rho_c = 0.336;
  unsigned int N = 20;
  unsigned int n_samples = 100;
    
  // fftw_init_threads();

  // #pragma omp parallel for shared(L,rho_c,save_folder)
  for(unsigned int i = 0 ; i < N ; i++) {
    double rho = rho_c + i*0.0005;
    string filename = save_folder + "S" + to_string(i) + ".txt";
    
    computeStructureFunction(L,rho,n_samples,filename);
  }
  
  return 0;
}


void computeStructureFunction(const unsigned int L, const double rho, const unsigned int n_samples, const string filename)
{
  SystemConfiguration system_config;
  system_config.directional_fep = 0;
  system_config.is_periodic = 1;
  system_config.vertical_boundaries = OPENED;
  system_config.horizontal_boundaries = OPENED;

  complex<double> i(0.,1.);
  unsigned int L2 = L * L;
  double rho2 = rho * rho;
  
  vector<double> S_k(L2,0.);
  
  fftw_complex* eta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*L2);
  fftw_complex* eta_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*L2);

  fftw_plan plan = fftw_plan_dft_1d(L2,eta,eta_hat,FFTW_FORWARD,FFTW_ESTIMATE);

  bool simulation_completed = false;
  do {
    ParticleSystem* system = createSystem(L,L,system_config);
    system->initSystem(rho);
  
    S_k = vector<double>(L2,0.);
  
    unsigned int m = 0;
    while((m < n_samples) && (system->getTotalNumberOfActiveEdges() > 0)) {
      double t = 0.;
      do {
        double dt = system->singleStep();
        t += dt;
      }while((t < L2) && (system->getTotalNumberOfActiveEdges() > 0));
    
      if(system->getTotalNumberOfActiveEdges() > 0) {
        for(unsigned int x = 0 ; x < L ; x++) {
          for(unsigned int y = 0 ; y < L ; y++) {
            Site* site = system->getSite(x,y);
            unsigned int i = L*x + y;
            system->isAnOccupiedSite(site) ? eta[i][0] = (1. - rho) : eta[i][0] = -rho;
            eta[i][1] = 0.;
          }
        }
    
        fftw_execute_dft(plan,eta,eta_hat);
    
        for(unsigned int i = 0 ; i < L2 ; i++) {
          S_k[i] += eta_hat[i][0]*eta_hat[i][0] + eta_hat[i][1]*eta_hat[i][1];
        }
      
        m++;
      }
    }
    
    delete system;
    
    simulation_completed = (m == n_samples);
  }while(simulation_completed == false);
  
  fftw_destroy_plan(plan);
  fftw_free(eta);
  fftw_free(eta_hat);
    
  ofstream file;
  file.open(filename,ios::out|ios::trunc);
  
  file<<rho<<endl;
  
  for(unsigned int x = 0 ; x < L ; x++) {
    for(unsigned int y = 0 ; y < L ; y++) {
      unsigned int i = L*x + y;
      S_k[i] /= (n_samples * L2);
      file<<x<<","<<y<<","<<S_k[i]<<endl;
    }
  }
  
  file.close();
}
