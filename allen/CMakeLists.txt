add_executable(compare compare.cpp)
target_link_libraries(compare ${OpenCV_LIBS} ${MORPH_LDFLAGS} ${ARMADILLO_LIBRARIES})
