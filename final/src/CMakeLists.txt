cmake_minimum_required(VERSION 3.12)
project(AcousticSimulation VERSION 1.0 LANGUAGES C)

find_package(MPI REQUIRED)
find_package(PETSc REQUIRED)

add_library(driver_functions driver_functions.c driver_functions.h)

add_executable(acoustic_simulation driver.c)

include_directories(${MPI_INCLUDE_PATH} ${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})

add_subdirectory(lib)

target_link_libraries(acoustic_simulation PUBLIC
	${MPI_C_LIBRARIES}
	${PETSC_LIBRARIES}
	driver_functions
	acoustic_problem)
