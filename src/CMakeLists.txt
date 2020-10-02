
add_executable(stalefish stalefish.cpp)
target_link_libraries(stalefish morphologica)

# To be removed; superceded by sfview
add_executable(viewer viewer.cpp)
target_link_libraries(viewer morphologica)

add_executable(sfview sfview.cpp)
target_link_libraries(sfview popt morphologica)

install(TARGETS stalefish sfview viewer DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
