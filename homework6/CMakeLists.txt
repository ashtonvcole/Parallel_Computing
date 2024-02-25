cmake_minimum_required(VERSION 3.12)
project(HW6 VERSION 1.0 LANGUAGES C)

find_package(MPI)

add_executable(HW6 HW6.c)

include_directories(${MPI_INCLUDE_PATH}) # Addition specific to MacOS

target_link_libraries(HW6 PUBLIC
	${MPI_C_LIBRARIES})
