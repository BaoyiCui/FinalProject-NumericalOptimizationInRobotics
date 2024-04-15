cmake_minimum_required(VERSION 3.0.2)
project(TOPP)

## Compile as C++11, supported in ROS Kinetic and newer
add_compile_options(-std=c++14)

set(MOSEK_DIR "/home/star/mosek/10.1")
set(MOSEK_LIBRARIES ${MOSEK_DIR}/tools/platform/linux64x86/bin)
set(MOSEK_INCLUDE_DIR ${MOSEK_DIR}/tools/platform/linux64x86/h)

## Find catkin macros and libraries
find_package(catkin REQUIRED COMPONENTS
  geometry_msgs
  roscpp
  std_msgs
  visualization_msgs
)

## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)
find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(PythonLibs REQUIRED)

include_directories(${MOSEK_DIR}/tools/platform/linux64x86/h)

find_library(MOSEK_BASE
    NAMES mosek64  # 库的名称，不带前缀"lib"和扩展名".so"
    HINTS ${MOSEK_LIBRARIES}  # 提供库文件所在的目录
)
find_library(MOSEK_FUSION
    NAMES fusion64  # 库的名称，不带前缀"lib"和扩展名".so"
    HINTS ${MOSEK_LIBRARIES}  # 提供库文件所在的目录
)

# 检查是否找到库文件
if(NOT MOSEK_BASE)
    message(FATAL_ERROR "MOSEK not found")
else()
    message(STATUS "Found MOSEK library: ${MOSEK_BASE}")
endif()

if(NOT MOSEK_FUSION)
    message(FATAL_ERROR "MOSEK fusion not found")
else()
    message(STATUS "Found MOSEK library: ${MOSEK_FUSION}")
endif()

# link_directories(${MOSEK_LIBRARIES_DIR})

catkin_package(
 INCLUDE_DIRS include
#  LIBRARIES TOPP
#  CATKIN_DEPENDS geometry_msgs roscpp std_msgs
#  DEPENDS system_lib
)


include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${MOSEK_INCLUDE_DIR}
  ${PYTHON_INCLUDE_DIRS}
)


add_executable(global_planner_node 
src/TOPP_node.cpp
src/global_planner.cpp
src/map_gen.cpp
)


target_link_libraries(global_planner_node
  ${catkin_LIBRARIES}
  CGAL::CGAL
  ${PYTHON_LIBRARIES}
  ${MOSEK_BASE}
  ${MOSEK_FUSION}
)

