
add_executable(stalefish stalefish.cpp)
target_link_libraries(stalefish morphologica)

# To be removed; superceded by sfview
add_executable(viewer viewer.cpp)
target_link_libraries(viewer morphologica)

add_executable(sfview sfview.cpp)
# Assume popt is in /usr/local/lib. Can't use target_link_directories without specifying minimum cmake version 3.13
link_directories(sfview PUBLIC /usr/local/lib)
target_link_libraries(sfview popt morphologica)

install(TARGETS stalefish sfview viewer DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
