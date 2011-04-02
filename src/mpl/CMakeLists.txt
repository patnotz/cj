cmake_minimum_required(VERSION 2.8)


SET(SOURCES
	mpl_main.cpp
  )

SET(HEADERS
  )

include_directories (
    .
)

#find_package(Boost 1.46.0 COMPONENTS program_options REQUIRED)
find_package(Boost 1.46.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Now build executables from this library, tailored main()s, and Trilinos

#add_library(mplLib ${SOURCES} ${HEADERS})
add_executable(mpl mpl_main.cpp)

#SET(mpl_LIBRARIES
#    mplLib
#)

#target_link_libraries(mpl ${mpl_LIBRARIES} ${Boost_LIBRARIES}) 
target_link_libraries(mpl ${Boost_LIBRARIES}) 
