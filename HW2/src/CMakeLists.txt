cmake_minimum_required(VERSION 3.12)
project(HW2 VERSION 1.0 LANGUAGES C)

find_package(MPI)

add_executable(HW2 HW2.c)
add_library(tools STATIC
	lib/tools.c
	lib/tools.h)
target_include_directories(tools PUBLIC
	lib
	${MPI_C_INCLUDE_DIRS})
target_link_libraries(HW2 PUBLIC
	tools
	${MPI_C_LIBRARIES})
install(TARGETS HW2 DESTINATION .)
