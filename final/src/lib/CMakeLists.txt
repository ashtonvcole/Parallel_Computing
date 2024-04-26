add_library(acoustic_problem "")

target_sources(acoustic_problem
	PRIVATE
	${CMAKE_CURRENT_LIST_DIR}/acoustic_problem.c
	PUBLIC
	${CMAKE_CURRENT_LIST_DIR}/acoustic_problem.h)

target_include_directories(acoustic_problem PUBLIC
	${CMAKE_CURRENT_LIST_DIR}
	${MPI_C_INCLUDE_DIRS})