project(qcj)

cmake_minimum_required(VERSION 2.8.0)
find_package(Qt4)

IF(NOT QT_FOUND)
    MESSAGE("Qt Not Found - skipping GUI build")
    RETURN()
ENDIF()

include(${QT_USE_FILE})

include_directories(
    .
    ${JsonCpp_INCLUDE_DIRS}
    )

SET(qcj_SRCS
    gui_main.cpp
    gui_mainWindow.cpp
)

set(qcj_MOC_SRCS
    gui_mainWindow.h
)

set(qcj_UIS
    #input.ui
    #main_window.ui
)

set(qcj_RCCS
#    input.qrc
)

qt4_wrap_cpp(qcj_MOCS ${qcj_MOC_SRCS})
qt4_wrap_ui(qcj_UIS_H ${qcj_UIS})
qt4_wrap_cpp(qcj_MOC_UI ${qcj_UIS_H})
qt4_add_resources(qcj_RCC_SRCS ${qcj_RCCS})
add_definitions(-DQT_NO_DEBUG)
add_executable(qcj
    ${qcj_SRCS}
    ${qcj_MOCS}
    ${qcj_MOC_UI}
    ${qcj_RCC_SRCS}
)

link_directories (
    ${JsonCpp_LIBRARY_DIRS}
)

target_link_libraries(qcj ${QT_LIBRARIES} ${Cj_LIBRARIES} pthread)
