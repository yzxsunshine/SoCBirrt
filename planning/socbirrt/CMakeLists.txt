if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing ")
endif()

if (DEFINED ENV{ROS_ROOT})

  cmake_minimum_required(VERSION 2.4.6)
  include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)
  rosbuild_init()
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
  rosbuild_add_library(socbirrt socbirrtmain.cpp TaskSpaceRegion.cpp socbirrt.cpp socbirrtproblem.cpp)
  target_link_libraries(socbirrt qhull newmat)

else (DEFINED ENV{ROS_ROOT})

  cmake_minimum_required (VERSION 2.4)
  if(POLICY CMP0015)
    cmake_policy(SET CMP0015 NEW)
  endif()
  project (socbirrt)

  if( NOT DEFINED ENV{PLUGIN_INSTALL_DIR} )
    set(PLUGIN_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../plugins")
  else( NOT DEFINED ENV{PLUGIN_INSTALL_DIR} )
    set(PLUGIN_INSTALL_DIR "$ENV{PLUGIN_INSTALL_DIR}")
  endif( NOT DEFINED ENV{PLUGIN_INSTALL_DIR} )

  include(FindOpenRAVE.cmake)

  if( NOT OPENRAVE_FOUND )
    message(FATAL_ERROR "could not find openrave, check for openrave-config in path")
  endif( NOT OPENRAVE_FOUND )

  link_directories(${OPENRAVE_LINK_DIRS})
  add_library(socbirrt SHARED socbirrtmain.cpp TaskSpaceRegion.cpp socbirrt.cpp socbirrtproblem.cpp)
  set_target_properties(socbirrt PROPERTIES
                        LINK_FLAGS ${OPENRAVE_LINK_FLAGS}
                        COMPILE_FLAGS ${OPENRAVE_CXXFLAGS})
  include_directories(${OPENRAVE_INCLUDE_DIRS} ${THIRDPARTY_DIR}/inc/newmat)
  target_link_libraries(socbirrt ${OPENRAVE_LIBRARY} qhull newmat)

  install(TARGETS socbirrt DESTINATION ${PLUGIN_INSTALL_DIR} )

endif (DEFINED ENV{ROS_ROOT})
