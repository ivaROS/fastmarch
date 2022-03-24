cmake_minimum_required(VERSION 2.8.3)
project(fastmarch)

add_compile_options(-Wall)


## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
)

find_package(OpenCV REQUIRED)

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if your package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME} ${PROJECT_NAME}_dt
  CATKIN_DEPENDS
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
  include
  ${OpenCV_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
)


add_library(${PROJECT_NAME} 
  external_projects/fastmarch/src/fastmarch.cc
  external_projects/fastmarch/src/heap.cc
  src/fastmarch_helper.cpp
  )
add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
#set_property(TARGET ego_circle_lib PROPERTY CXX_STANDARD 14)
#target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_include_directories(${PROJECT_NAME} PRIVATE 
  external_projects/fastmarch/include
  external_projects/include
  )
target_compile_definitions(${PROJECT_NAME} PRIVATE _MATRIX_ _DEBUG_ _IMFLOAT_ _INTERP_RI_  
  #_RESTRICT_  _USENB_ #These do not appear to do anything
  )
target_link_libraries(${PROJECT_NAME}  
  ${OpenCV_LIBRARIES}
  ${catkin_LIBRARIES}
  )
  
  
  

add_library(${PROJECT_NAME}_dt 
  src/dt_helper.cpp
  )
add_dependencies(${PROJECT_NAME}_dt ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
#set_property(TARGET ego_circle_lib PROPERTY CXX_STANDARD 14)
#target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_include_directories(${PROJECT_NAME}_dt PRIVATE 
  external_projects/dt
  )
target_link_libraries(${PROJECT_NAME}_dt 
  ${OpenCV_LIBRARIES}
  ${catkin_LIBRARIES}
  )
  
  
  

add_executable(${PROJECT_NAME}_test
    src/test.cpp
  )
target_link_libraries(${PROJECT_NAME}_test
  ${PROJECT_NAME}
  ${PROJECT_NAME}_dt
  ${OpenCV_LIBRARIES}
  ${catkin_LIBRARIES}
  )



add_executable(${PROJECT_NAME}_stest
    src/stest.cpp
  )
# target_link_libraries(${PROJECT_NAME}_stest
#   ${PROJECT_NAME}
#   ${PROJECT_NAME}_dt
#   ${OpenCV_LIBRARIES}
#   ${catkin_LIBRARIES}
#   )
#   
target_include_directories(${PROJECT_NAME}_stest PRIVATE 
  external_projects/dt/dt
  )
target_compile_options(${PROJECT_NAME}_stest PRIVATE
 -march=native -O3 -ftree-vectorize -fopt-info-vec-all=vectorized.txt #-fassociative-math -fno-signed-zeros -fno-trapping-math  -fno-signaling-nans -ffast-math 
 ) 
target_link_libraries(${PROJECT_NAME}_stest
  ${PROJECT_NAME}_dt
  ${OpenCV_LIBRARIES}
  ${catkin_LIBRARIES}
  )
  
  
