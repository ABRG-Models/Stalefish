
add_executable(stalefish stalefish.cpp)
target_link_libraries(stalefish ${OpenCV_LIBS} ${MORPH_LDFLAGS} ${ARMADILLO_LIBRARIES})

add_executable(viewer viewer.cpp)
target_link_libraries(viewer ${OpenCV_LIBS} ${MORPH_LDFLAGS} ${ARMADILLO_LIBRARIES} ${GLFW_LIBRARIES} ${OPENGL_LIBRARIES})
