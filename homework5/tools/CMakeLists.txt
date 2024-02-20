add_library(tools "")

target_sources(tools
	PRIVATE
	${CMAKE_CURRENT_LIST_DIR}/tools.c
	PUBLIC
	${CMAKE_CURRENT_LIST_DIR}/tools.h)

target_include_directories(tools PUBLIC
	${CMAKE_CURRENT_LIST_DIR}
	${MPI_C_INCLUDE_DIRS})
