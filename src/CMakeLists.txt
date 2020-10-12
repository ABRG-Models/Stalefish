include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${ARMADILLO_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIRS})
include_directories(${GLFW3_INCLUDE_DIR})

add_executable(stalefish stalefish.cpp)
target_link_libraries(stalefish ${OpenCV_LIBS} ${HDF5_C_LIBRARIES} jsoncpp_lib)
target_link_libraries(stalefish ${ARMADILLO_LIBRARY} ${ARMADILLO_LIBRARIES} ${LAPACK_LIBRARIES})

# To be removed; superceded by sfview
add_executable(viewer viewer.cpp)

add_executable(sfview sfview.cpp)
target_link_libraries(sfview ${OpenCV_LIBS} ${HDF5_C_LIBRARIES} jsoncpp_lib OpenGL::GL)
target_link_libraries(sfview ${ARMADILLO_LIBRARY} ${ARMADILLO_LIBRARIES} ${LAPACK_LIBRARIES})
target_link_libraries(sfview glfw)
# Assume popt is in /usr/local/lib. Can't use target_link_directories without specifying minimum cmake version 3.13
target_link_libraries(sfview -L/usr/local/lib popt)

install(TARGETS stalefish sfview viewer DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
