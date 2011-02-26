cmake_minimum_required(VERSION 2.8)

add_subdirectory(jsoncpp-0.5.0)
#
SET(SOURCES
	cj_main.cpp
	drive_simulation.cpp
	messages.cpp 
	mesh_manager.cpp
	stk_mesh.cpp
  )

SET(HEADERS
  )

set(CJ_TRILINOS_LIBS ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

include_directories (${Trilinos_INCLUDE_DIRS}
                     .
                     ${Trilinos_TPL_INCLUDE_DIRS})

link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

add_library(cjLib ${SOURCES} ${HEADERS})

find_package(Boost 1.46.0 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Now build executables from this library, tailored main()s, and Trilinos

add_executable(Cj cj_main.cpp)

target_link_libraries(Cj cjLib) 
target_link_libraries(Cj stk_io) 
target_link_libraries(Cj stk_mesh_base) 
target_link_libraries(Cj stk_mesh_fem) 
target_link_libraries(Cj stk_util_diag) 
target_link_libraries(Cj stk_util_env) 
target_link_libraries(Cj stk_util_parallel) 
target_link_libraries(Cj stk_util_util) 
target_link_libraries(Cj shards) 
target_link_libraries(Cj curl hdf5 hdf5_hl exodus netcdf z)
target_link_libraries(Cj ${Boost_LIBRARIES})
