# AUTOMOC requires 2.8.6. Per-target include directories would require 2.8.8
cmake_minimum_required(VERSION 2.8.8)
project(Stalefish)
set(VERSION 0.1)
set(SOVERSION 1.0)

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
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c++11")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -g -std=c++11")
endif()

# Lib finding
find_package(OpenCV REQUIRED)
find_package(X11 REQUIRED)

if(APPLE)
  find_package (HDF5 REQUIRED)
  set (HDF5LIBLINK ${HDF5_C_LIBRARY_hdf5})
  set (JSONLIBLINK /usr/local/lib/libjsoncpp.dylib)
  # If omp is to be used:
  # include_directories(/opt/local/include/libomp)
  # link_libraries(-L/opt/local/lib/libomp omp)
else()
  set (LIB_TYPE SHARED) # STATIC or SHARED
  string(TOLOWER ${LIB_TYPE} SEARCH_TYPE)
  find_package (HDF5 NAMES hdf5 COMPONENTS C ${SEARCH_TYPE} REQUIRED)
  set (HDF5LIBLINK ${HDF5_C_${LIB_TYPE}_LIBRARY})
  set (JSONLIBLINK jsoncpp)
endif(APPLE)

include_directories(${OpenCV_INCLUDE_DIRS} ${X11_INC_SEARCH_PATH} ${HDF5_INCLUDE_DIR})

# Library code is compiled up as a shared library in lib/ (could also
# be compiled static if needed)
add_subdirectory(src)

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
