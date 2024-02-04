cmake_minimum_required(VERSION 3.12)
project(HW3 VERSION 1.0 LANGUAGES C)

find_package(MPI)

add_executable(HW3 HW3.c)

add_subdirectory(tools)

target_link_libraries(HW3 PUBLIC
	tools
	${MPI_C_LIBRARIES})
