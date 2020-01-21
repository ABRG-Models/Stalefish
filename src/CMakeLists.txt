
add_executable(stalefish stalefish.cpp)
#target_compile_definitions(stalefish PUBLIC)
target_link_libraries(stalefish ${OpenCV_LIBS} ${X11_X11_LIB})
# To add ${HDF5LIBLINK} ${JSONLIBLINK}
