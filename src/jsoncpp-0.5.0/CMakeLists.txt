project(jsoncpp-0.5.0)
cmake_minimum_required(VERSION 2.8.0)
#
SET(SOURCES
    json_reader.cpp
    json_value.cpp
    json_writer.cpp
  )

SET(HEADERS
    json_batchallocator.h
    json_internalarray.inl
    json_internalmap.inl
    json_valueiterator.inl
    json/autolink.h
    json/config.h
    json/features.h
    json/forwards.h
    json/json.h
    json/reader.h
    json/value.h
    json/writer.h
  )

include_directories (.)

link_directories ()

add_library(jsoncpp ${SOURCES} ${HEADERS})

if(GTEST_FOUND)
   include_directories(${GTEST_INCLUDE_DIRS})
   add_executable(json_unit_tests json_unit_tests.cpp)
   target_link_libraries(json_unit_tests jsoncpp)
   target_link_libraries(json_unit_tests ${GTEST_BOTH_LIBRARIES})
   add_test(AllJsonTests json_unit_tests)
endif()
