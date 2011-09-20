cmake_minimum_required(VERSION 2.8)


set(JsonCpp_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/jsoncpp-0.5.0/)
set(JsonCpp_LIBRARY_DIRS ${CMAKE_CURRENT_BINARY_DIR}/jsoncpp-0.5.0/)
set(JsonCpp_LIBRARY jsoncpp)
add_subdirectory(jsoncpp-0.5.0)

#
SET(SOURCES
	cj_main.cpp
	drive_simulation.cpp
	log.cpp
	messages.cpp 
	mesh_manager.cpp
	stk_mesh.cpp
  )

SET(HEADERS
  )

set(CJ_TRILINOS_LIBS ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

include_directories (
    ${Trilinos_INCLUDE_DIRS}
    .
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${JsonCpp_INCLUDE_DIRS}
)

link_directories (
    ${Trilinos_LIBRARY_DIRS} 
    ${Trilinos_TPL_LIBRARY_DIRS} 
    ${JsonCpp_LIBRARY_DIRS}
)

add_library(cjLib ${SOURCES} ${HEADERS})

find_package(Boost 1.46.0 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Now build executables from this library, tailored main()s, and Trilinos

add_executable(Cj cj_main.cpp)

SET(Cj_LIBRARIES
    cjLib
#    stk_io
    stk_mesh_base
    stk_mesh_fem
    stk_util_diag
    stk_util_env
    stk_util_parallel
    stk_util_util
    shards
    curl hdf5 hdf5_hl exodus netcdf z
    ${JsonCpp_LIBRARY}
)

target_link_libraries(Cj ${Cj_LIBRARIES} ${Boost_LIBRARIES}) 

add_subdirectory(gui)
