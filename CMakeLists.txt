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

# Tell clang to be quiet about brace initialisers
if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-missing-braces")
endif()

# Following `cmake --help-policy CMP0072`
set(OpenGL_GL_PREFERENCE "GLVND")

# Additional GL compiler flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")
if(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL_SILENCE_DEPRECATION")
endif()

# Lib finding
find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
find_package(HDF5 REQUIRED)
find_package(jsoncpp REQUIRED)
find_package(LAPACK REQUIRED)
find_package(glfw3 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(morphologica REQUIRED)

# Library code is compiled up as a shared library in lib/ (could also
# be compiled static if needed)
add_subdirectory(src)
add_subdirectory(allen)

# Unit testing using the ctest framework
# enable_testing()
# add_subdirectory(tests)

# For debugging of variables:
option(DEBUG_VARIABLES ON)
if(DEBUG_VARIABLES)
  get_cmake_property(_variableNames VARIABLES)
  foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
  endforeach()
endif(DEBUG_VARIABLES)
