
add_executable(stalefish stalefish.cpp)
target_link_libraries(stalefish morphologica)

add_executable(viewer viewer.cpp)
target_link_libraries(viewer morphologica)

add_executable(viewer2 viewer2.cpp)
target_link_libraries(viewer2 popt morphologica)

install(TARGETS viewer stalefish DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
