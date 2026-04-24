cmake_minimum_required(VERSION 3.10)
project(libtaph VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Library target ────────────────────────────────────────────────────────────
add_library(taph
    src/damage_estimation.cpp
    src/damage_c_api.cpp
    src/length_stratified_profile.cpp
    src/log_length_gmm.cpp
    src/library_interpretation.cpp
)

target_include_directories(taph PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_options(taph PRIVATE -O2 -Wall -Wextra)

# ── Install ───────────────────────────────────────────────────────────────────
include(GNUInstallDirs)

install(TARGETS taph
    EXPORT taph-targets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# CMake package config (for find_package support)
install(EXPORT taph-targets
    FILE taph-targets.cmake
    NAMESPACE taph::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/taph
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/taph-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/taph-config.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/taph"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/taph-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/taph-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/taph-config-version.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/taph"
)

# ── Tests ─────────────────────────────────────────────────────────────────────
option(TAPH_BUILD_TESTS "Build libtaph unit tests" OFF)
if(TAPH_BUILD_TESTS)
    enable_testing()
    add_executable(test_damage_context_profile tests/test_damage_context_profile.cpp)
    target_link_libraries(test_damage_context_profile PRIVATE taph)
    target_compile_options(test_damage_context_profile PRIVATE -O2 -Wall -Wextra)
    add_test(NAME damage_context_profile COMMAND test_damage_context_profile)
endif()
