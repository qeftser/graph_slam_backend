cmake_minimum_required(VERSION 3.6)
project(graph_slam_backend)

add_library(graph_slam_backend SHARED src/graph_slam.c src/hash.c src/homogeneous.c src/solver.c src/set.c
                                      src/graph_slam.h src/hash.h src/homogeneous.h src/solver.h src/set.h)
add_library(GraphSlamBackend::graph_slam_backend ALIAS graph_slam_backend)

include(GNUInstallDirs)

install(
   TARGETS graph_slam_backend
   EXPORT GraphSlamBackendTargets
   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(
   EXPORT GraphSlamBackendTargets
   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GraphSlamBackend
   NAMESPACE GraphSlamBackend::
   FILE GraphSlamBackendTargets.cmake
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
   "Config.cmake.in"
   "GraphSlamBackendConfig.cmake"
   INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GraphSlamBackend
   PATH_VARS
      CMAKE_INSTALL_LIBDIR
)

write_basic_package_version_file(
   ${CMAKE_CURRENT_BINARY_DIR}/GraphSlamBackendConfigVersion.cmake
   VERSION 1.0.0
   COMPATIBILITY SameMajorVersion
)

install(
   FILES "${CMAKE_CURRENT_BINARY_DIR}/GraphSlamBackendConfig.cmake"
         "${CMAKE_CURRENT_BINARY_DIR}/GraphSlamBackendConfigVersion.cmake"
   DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GraphSlamBackend"
)

