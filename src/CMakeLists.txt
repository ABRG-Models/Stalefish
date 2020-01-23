
add_executable(stalefish stalefish.cpp)
#target_compile_definitions(stalefish PUBLIC)
target_link_libraries(stalefish ${OpenCV_LIBS} ${MORPH_LDFLAGS})
# To add ${HDF5LIBLINK} ${JSONLIBLINK}  ${X11_X11_LIB}
