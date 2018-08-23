cmake_minimum_required( VERSION 3.5 )
project(SCC)

set( CMAKE_BUILD_TYPE "Release" )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_CXX_STANDARD 11)

# thirdparty libs 
find_package(Eigen3 REQUIRED)
find_package(OpenCV 2.4.0 REQUIRED)
include_directories(
#${EIGEN3_INCLUDE_DIR}
${OpenCV_INCLUDE_DIRS}
${EIGEN3_INCLUDE_DIRS}
)
add_subdirectory(SCC)

add_executable(example example.cpp)
target_link_libraries(example ${OpenCV_LIBS} scc_motion_seg )

