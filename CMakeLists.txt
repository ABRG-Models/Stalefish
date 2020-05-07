# 3.1 required for the CMAKE_CXX_STANDARD variable
cmake_minimum_required(VERSION 3.1)
project(Stalefish)
set(VERSION 0.1)
set(SOVERSION 1.0)

set(CMAKE_CXX_STANDARD 17)

message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  (This can be changed with `cmake -DCMAKE_INSTALL_PREFIX=/some/place`")

# From CMAKE_SYSTEM work out which of __OSX__, __GLN__, __NIX__ or
# __WIN__ are required in case code needs to know.
message(STATUS "Operating system: " ${CMAKE_SYSTEM})
if(CMAKE_SYSTEM MATCHES Linux.*)
  set(STALEFISH_HOST_DEFINITION "-D__GLN__")
elseif(CMAKE_SYSTEM MATCHES BSD.*)
  set(STALEFISH_HOST_DEFINITION "-D__NIX__")
elseif(APPLE)
  set(STALEFISH_HOST_DEFINITION "-D__OSX__")
elseif(CMAKE_SYSTEM MATCHES Win.*)
  set(STALEFISH_HOST_DEFINITION "-D__WIN__")
else()
  message(ERROR "Operating system not supported: " ${CMAKE_SYSTEM})
endif()
# Add the host definition to CFLAGS
set(CMAKE_CXX_FLAGS "${STALEFISH_HOST_DEFINITION}")
set(CMAKE_C_FLAGS "${STALEFISH_HOST_DEFINITION}")

if(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -Wfatal-errors")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -Wfatal-errors -g")
endif()

# Following `cmake --help-policy CMP0072`
set(OpenGL_GL_PREFERENCE "GLVND")

# Additional GL compiler flags.
# maybe -DGL_GLEXT_PROTOTYPES
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")

# Lib finding
find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(X11 REQUIRED)
find_package(LAPACK REQUIRED)
find_package(PkgConfig)

# glfw3 provides a cmake find_package script. Because it doesn't set
# lib paths, we go on and use pkg-config to do this.
find_package(glfw3 QUIET)
if(${glfw3_FOUND})
  message(INFO ": cmake's find_package() found glfw3, using pkg-config to get details...")
  if(${PKG_CONFIG_FOUND})
    # Use pkg-config to check GLFW3 is present:
    pkg_search_module(GLFW REQUIRED glfw3)
    if(${GLFW_FOUND})
      message(INFO ": pkg-config found glfw(3). Excellent.")
    else(${GLFW_FOUND})
      message(ERROR ": glfw(3) was NOT found by pkg-config. Uh oh.")
    endif(${GLFW_FOUND})
  else(${PKG_CONFIG_FOUND})
    message(WARNING "There's no pkg-config on this system to find exact details for glfw3. You may need to try `sudo apt-get install pkg-config` or similar")
  endif(${PKG_CONFIG_FOUND})
endif(${glfw3_FOUND})

if(APPLE)
  find_package (Armadillo 9.20 HINTS /opt/local/lib REQUIRED)
else()
  find_package(Armadillo REQUIRED)
endif(APPLE)

include_directories(${OpenCV_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR} ${X11_INC_SEARCH_PATH})

# Find libmorphologica using pkgconfig
if(PKG_CONFIG_FOUND)
  message(STATUS "Have pkg_config, searching for libmorphologica...")

  if(APPLE)
    # Manually set PKG_CONFIG_PATH on Mac, as mac ports will by default search only in /opt/local
    set(ENV{PKG_CONFIG_PATH}  "ENV{PKG_CONFIG_PATH}:${CMAKE_INSTALL_PREFIX}/lib/pkgconfig")
  endif(APPLE)

  pkg_check_modules(MORPH REQUIRED libmorphologica)
  if (MORPH_FOUND)
    message(STATUS "pkg_config found libmorphologica")
    include_directories(${MORPH_INCLUDEDIR})
    include_directories(${MORPH_INCLUDEDIR}/morph)
    message(STATUS "Included morph includedir:  ${MORPH_INCLUDEDIR} and ${MORPH_INCLUDEDIR}/morph")
  else()
    message(FATAL_ERROR "You need libmorphologica from github.com/ABRG_Models/morphologica")
  endif(MORPH_FOUND)
else()
  message(WARNING "There's no pkg-config on this system to check for libmorphologica. "
    "You can set MORPH_LIB with `cmake -DMORPH_LIB=/pathto/lib` ...")
endif()

# Library code is compiled up as a shared library in lib/ (could also
# be compiled static if needed)
add_subdirectory(src)
add_subdirectory(allen)

# Unit testing using the ctest framework
# enable_testing()
# add_subdirectory(tests)

# For debugging of variables:
set(DEBUG_VARIABLES ON)
if(DEBUG_VARIABLES)
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endif(DEBUG_VARIABLES)
