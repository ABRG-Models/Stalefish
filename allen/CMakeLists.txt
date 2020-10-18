include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIRS})
add_executable(compare compare.cpp)
target_link_libraries(compare ${OpenCV_LIBS} ${ARMADILLO_LIBRARIES})
