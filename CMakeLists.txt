cmake_minimum_required(VERSION 3.10)
project(Fluid3D)

##### Setup ######
set(CMAKE_CXX_STANDARD 98)

##### External Libraries ######
find_package(OpenGL)
find_package(GLUT)

#find_package(PkgConfig)
#
## GTest Setup
#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS} include)
#
## GMock Setup
#pkg_check_modules(GMOCK REQUIRED)
#include_directories(${GMOCK_INCLUDE_DIRS})
#
## General Testing Setup
#set(TESTING_LIBS ${GTEST_BOTH_LIBRARIES} pthread)
#enable_testing()
#
###### External Libraries ######
#add_subdirectory(lib)

###### GTK ######
#find_package(PkgConfig REQUIRED)
#pkg_check_modules("GTKMM" REQUIRED "gtkmm-3.0")
#include_directories("${GTKMM_INCLUDE_DIRS}")
#link_directories("${GTKMM_LIBRARY_DIRS}")

###### Header Files ######
include_directories(
    ./
    ${GLU_INCLUDES}
    )

##### Executables #####

# To avoid naming conflicts (can't have a "viewpls3D" directory AND file in
# the same directory) output all the executables to their own folder
set(EXECUTABLE_OUTPUT_PATH "bin")

## Solver ##
# Find all the source (.cpp, .c) files in the current directory
aux_source_directory(. SOLVER_SRC)

# Build the executable
add_executable(solver
        ${SOLVER_SRC}
        )

## viewpls3D ##
# Find all the source (.cpp, .c) files in the current directory
aux_source_directory(./viewpls3D VIEW_PLS_3D_SRC)

# Build the executable
add_executable(viewpls3D
        ${VIEW_PLS_3D_SRC}
        )

# Link the executable to external libraries
target_link_libraries(viewpls3D
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    )

