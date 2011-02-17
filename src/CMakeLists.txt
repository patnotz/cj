
#
SET(SOURCES
	cj_main.cpp
	messages.cpp 
	mesh_manager.cpp
	stk_mesh.cpp
  )

SET(HEADERS
  )

set(CJ_TRILINOS_LIBS ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

include_directories (${Trilinos_INCLUDE_DIRS}
                     ../include 
                     ${Trilinos_TPL_INCLUDE_DIRS})

link_directories (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

add_library(cjLib ${SOURCES} ${HEADERS})
# Now build executables from this library, tailored main()s, and Trilinos

add_executable(Cj cj_main.cpp)

target_link_libraries(Cj          cjLib ${ALB_TRILINOS_LIBS} ${ALB_CUBIT_LIBS}) 
