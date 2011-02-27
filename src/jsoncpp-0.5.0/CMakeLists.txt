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

include_directories (${JsonCpp_INCLUDE_DIRS})
link_directories (${JsonCpp_LIBRARY_DIRS})

add_library(${JsonCpp_LIBRARY} ${SOURCES} ${HEADERS})

if(GTEST_FOUND)
   enable_testing(true)
   include_directories(${GTEST_INCLUDE_DIRS} ${JsonCpp_INCLUDE_DIRS})
   add_executable(json_unit_tests json_unit_tests.cpp)
   target_link_libraries(json_unit_tests ${JsonCpp_LIBRARY})
   target_link_libraries(json_unit_tests ${GTEST_BOTH_LIBRARIES})
   add_test(JsonUnitTests json_unit_tests)
endif()
