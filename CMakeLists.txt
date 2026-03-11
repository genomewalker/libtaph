cmake_minimum_required(VERSION 3.10)
project(libdart-damage VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Library target ────────────────────────────────────────────────────────────
add_library(dart-damage
    src/damage_estimation.cpp
    src/damage_c_api.cpp
)

target_include_directories(dart-damage PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_options(dart-damage PRIVATE -O2 -Wall -Wextra)

# ── Install ───────────────────────────────────────────────────────────────────
include(GNUInstallDirs)

install(TARGETS dart-damage
    EXPORT dart-damage-targets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# CMake package config (for find_package support)
install(EXPORT dart-damage-targets
    FILE dart-damage-targets.cmake
    NAMESPACE dart::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dart-damage
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/dart-damage-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/dart-damage-config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/dart-damage"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/dart-damage-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/dart-damage-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/dart-damage-config-version.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/dart-damage"
)
