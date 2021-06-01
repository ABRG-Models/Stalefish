# Define collections of includes that have to be made for morphologica
set(MORPH_INC_CORE ${ARMADILLO_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIRS} ${HDF5_INCLUDE_DIR})
set(MORPH_INC_GL ${OpenCV_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR} ${GLFW3_INCLUDE_DIR} ${FREETYPE_INCLUDE_DIRS})
include_directories(${MORPH_INC_CORE} ${MORPH_INC_GL})

# Define collections of links that have to be made for morphologica
set(MORPH_LIBS_CORE ${ARMADILLO_LIBRARY} ${ARMADILLO_LIBRARIES} ${LAPACK_LIBRARIES} ${HDF5_C_LIBRARIES} jsoncpp_lib)
set(MORPH_LIBS_GL ${OpenCV_LIBS} OpenGL::GL glfw Freetype::Freetype)

add_executable(stalefish stalefish.cpp)
target_link_libraries(stalefish ${MORPH_LIBS_CORE} ${OpenCV_LIBS})

# On Mac, add -headerpad_max_install_names so that the executable can have its paths reconfigured for packaging
if(APPLE)
  set(MAC_HEADERPAD -headerpad_max_install_names)
endif()

add_executable(sfresave sfresave.cpp)
target_link_libraries(sfresave ${MORPH_LIBS_CORE} ${OpenCV_LIBS} ${MAC_HEADERPAD})

add_executable(sfgetjson sfgetjson.cpp)
target_link_libraries(sfgetjson ${MORPH_LIBS_CORE} ${OpenCV_LIBS} ${MAC_HEADERPAD})

add_executable(sfview sfview.cpp)
# Assume popt is in /usr/local/lib. Can't use target_link_directories without specifying minimum cmake version 3.13
if(APPLE AND OpenMP_CXX_FOUND)
  target_link_libraries(sfview OpenMP::OpenMP_CXX ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL} -L/usr/local/lib popt ${MAC_HEADERPAD})
else()
  target_link_libraries(sfview ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL} -L/usr/local/lib popt)
endif()

add_executable(sfview_unwrap_fig sfview_unwrap_fig.cpp)
if(APPLE AND OpenMP_CXX_FOUND)
  target_link_libraries(sfview_unwrap_fig OpenMP::OpenMP_CXX ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL} -L/usr/local/lib popt ${MAC_HEADERPAD})
else()
  target_link_libraries(sfview_unwrap_fig ${MORPH_LIBS_CORE} ${MORPH_LIBS_GL} -L/usr/local/lib popt)
endif()

install(TARGETS stalefish sfview sfresave sfgetjson DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
