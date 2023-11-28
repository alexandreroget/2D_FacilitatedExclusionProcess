cmake_minimum_required(VERSION 3.12)
project(2D_FEP)

include_directories(include)

set(SOURCES
    src/ParticleSystem.cpp
    src/Site.cpp
    src/Edge.cpp
)

add_executable(computeAbsorptionTime simulations/Absorption_time/computeAbsorptionTime.cpp ${SOURCES})
set_target_properties(computeAbsorptionTime PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simulations/Absorption_time)

add_executable(computeActiveSitesDensity simulations/Active_sites_density/computeActiveSitesDensity.cpp ${SOURCES})
set_target_properties(computeActiveSitesDensity PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simulations/Active_sites_density)

add_executable(computeDiffusionCurrent simulations/Diffusion_coefficient/computeDiffusionCurrent.cpp ${SOURCES})
set_target_properties(computeDiffusionCurrent PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simulations/Diffusion_coefficient)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/findFFTW")
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})

add_executable(computeStructureFunction simulations/Structure_function/computeStructureFunction.cpp ${SOURCES})
set_target_properties(computeStructureFunction PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/simulations/Structure_function)

target_link_libraries(computeStructureFunction PRIVATE ${FFTW_LIBRARIES})
