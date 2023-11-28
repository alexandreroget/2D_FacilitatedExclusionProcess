# OVERVIEW

**2D_FacilitatedExclusionProcess** is a C++ software that simulates an exclusion model on a 2-dimensional lattice, where active particles jump randomly at rate 1 to each empty nearest neighbor. A particle is considered active if at least one of its neighboring sites is also occupied.

Four simulation programs are currently available in the **simulations** folder. Each subfolder contains a C++ source code, a short description of the code and a Python post process script.

If you wish to conduct your own simulations, please consult the documentation that describes functions, modules, and objects included in  **2D_FacilitatedExclusionProcess**.

# INSTALLATION GUIDE

To get started with the code, follow these steps for a successful installation:

1. Clone the repository to your local machine using the following command:
    ```bash
    git clone https://gitlab.inria.fr/aroget/2d_facilitatedexclusionprocess.git
    ```

2. Navigate to the project directory:
    ```bash
    cd 2D_FacilitatedExclusionProcess
    ```

3. Create a 'build' directory for the build process:
    ```bash
    mkdir build
    ```

4. Change into the 'build' directory:
    ```bash
    cd build
    ```

5. Configure and Build the environment:
   Ensure that you have **CMake** installed for configuring the build environment. If not, you can install it using:
    ```bash
    sudo apt-get install cmake  # For Debian/Ubuntu
    ```
   The makefile created by CMake compiles all source files in the **simulations** folder, generating the corresponding executable files.

   Please note that the **FFTW** library is required specifically for the code `computeStructureFunction.cpp`. Ensure that you have it installed before proceeding. You can install FFTW using:
    ```bash
    sudo apt-get install libfftw3-dev  # For Debian/Ubuntu
    ```

   You can then configure and build the project:
    ```bash
    cmake ..
    make
    ```

   After the compilation is complete, you can proceed to run a simulation.

# RUN A SIMULATION

To run a simulation, follow these steps:

1. Navigate to the 'simulations' folder within the project directory:
    ```bash
    cd ../simulations
    ```

2. Choose the simulation you want to run and locate its corresponding executable file. For example, if you want to run the simulation that computes the Absorption Time, use:
    ```bash
    cd Absorption_time
    ./computeAbsorptionTime
    ```

3. Once the simulation is complete, run the Python post-process script to analyze the results.
    ```bash
    python3 post_process.py
    ```

Feel free to explore different simulations and adjust parameters as needed.

# CONTACTS

**2D_FacilitatedExclusionProcess** has been developed in the Paradyse team of Inria by Alexandre Roget. You can contact him at alexandre.roget@inria.fr.
