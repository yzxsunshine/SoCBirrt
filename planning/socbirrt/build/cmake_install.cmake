# Install script for directory: /home/jason/RBE595/CBiRRT/planning/socbirrt

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}/home/jason/RBE595/CBiRRT/planning/socbirrt/../../plugins/libsocbirrt.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/jason/RBE595/CBiRRT/planning/socbirrt/../../plugins/libsocbirrt.so")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/jason/RBE595/CBiRRT/planning/socbirrt/../../plugins/libsocbirrt.so"
         RPATH "")
  ENDIF()
  list(APPEND CPACK_ABSOLUTE_DESTINATION_FILES
   "/home/jason/RBE595/CBiRRT/planning/socbirrt/../../plugins/libsocbirrt.so")
FILE(INSTALL DESTINATION "/home/jason/RBE595/CBiRRT/planning/socbirrt/../../plugins" TYPE SHARED_LIBRARY FILES "/home/jason/RBE595/CBiRRT/planning/socbirrt/build/libsocbirrt.so")
  IF(EXISTS "$ENV{DESTDIR}/home/jason/RBE595/CBiRRT/planning/socbirrt/../../plugins/libsocbirrt.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/jason/RBE595/CBiRRT/planning/socbirrt/../../plugins/libsocbirrt.so")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/jason/RBE595/CBiRRT/planning/socbirrt/../../plugins/libsocbirrt.so")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/home/jason/RBE595/CBiRRT/planning/socbirrt/build/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/home/jason/RBE595/CBiRRT/planning/socbirrt/build/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
