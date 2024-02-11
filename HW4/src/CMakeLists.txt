cmake_minimum_required(VERSION 3.12)
project(HW4 VERSION 1.0 LANGUAGES C)

find_package(MPI)

add_executable(HW4 HW4.c)

include_directories(${MPI_INCLUDE_PATH}) # Addition specific to MacOS

target_link_libraries(HW4 PUBLIC
	${MPI_C_LIBRARIES})
