project(qcj)

cmake_minimum_required(VERSION 2.8.0)
find_package(Qt4 REQUIRED)
include(${QT_USE_FILE})

include_directories(.)

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
target_link_libraries(qcj ${QT_LIBRARIES} pthread)
