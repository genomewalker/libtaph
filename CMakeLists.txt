cmake_minimum_required(VERSION 3.10)
project(libtaph VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Dev setup: install git hooks ──────────────────────────────────────────────
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_SOURCE_DIR}/scripts/hooks/pre-commit"
            "${CMAKE_SOURCE_DIR}/.git/hooks/pre-commit"
    )
    execute_process(
        COMMAND chmod +x "${CMAKE_SOURCE_DIR}/.git/hooks/pre-commit"
    )
endif()

# ── Library target ────────────────────────────────────────────────────────────
add_library(taph
    src/damage_estimation_update.cpp
    src/damage_estimation_finalize.cpp
    src/damage_estimation_finalize_init.cpp
    src/damage_estimation_finalize_decay.cpp
    src/damage_estimation_finalize_oxidation.cpp
    src/damage_estimation_finalize_context.cpp
    src/damage_estimation_finalize_libtype.cpp
    src/damage_estimation_finalize_dmax.cpp
    src/damage_estimation_finalize_preservation.cpp
    src/damage_estimation_channels.cpp
    src/bulk_damage_model.cpp
    src/read_ancient_llr.cpp
    src/finalize_tau.cpp
    src/finalize_scission.cpp
    src/finalize_epsilon.cpp
    src/finalize_gc_depletion.cpp
    src/finalize_oxidation_comovement.cpp
    src/damage_c_api.cpp
    src/length_stratified_profile.cpp
    src/log_length_gmm.cpp
    src/library_interpretation_hex.cpp
    src/library_interpretation_scores.cpp
    src/library_interpretation_summary.cpp
    src/damage_profiler.cpp
    src/profile_json.cpp
    src/fraction_fitting.cpp
    src/llr_table.cpp
    src/ancient_fraction.cpp
    src/oxog_score.cpp
)

add_library(taph::taph ALIAS taph)

target_include_directories(taph PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

set(TAPH_MARCH "" CACHE STRING "Value passed to -march=; empty = portable. Set to 'native' for local optimisation.")
target_compile_options(taph PRIVATE -O2 -Wall -Wextra)
if(NOT "${TAPH_MARCH}" STREQUAL "")
    target_compile_options(taph PRIVATE -march=${TAPH_MARCH})
endif()

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

    # Taphonomy A-bar: the bulk damage kernel must differ for ds vs ss at the 5' terminus (Wave-3).
    add_executable(test_ss_overhang_kernel tests/test_ss_overhang_kernel.cpp)
    target_link_libraries(test_ss_overhang_kernel PRIVATE taph)
    target_compile_options(test_ss_overhang_kernel PRIVATE -O2 -Wall -Wextra)
    add_test(NAME ss_overhang_kernel COMMAND test_ss_overhang_kernel)

    add_executable(test_damage_regressions tests/test_damage_regressions.cpp)
    target_link_libraries(test_damage_regressions PRIVATE taph)
    target_compile_options(test_damage_regressions PRIVATE -O2 -Wall -Wextra)
    add_test(NAME damage_regressions COMMAND test_damage_regressions)

    # Layer-0 count-table golden gate: re-profiles fixtured samples and diffs the pre-clamp count
    # tables (catches cap-masked z regressions the emitted-JSON gate cannot see). Self-skips when the
    # external panel inputs / fqdup binary are absent, so it is green in a data-less CI checkout.
    add_test(NAME count_table_golden
             COMMAND ${CMAKE_COMMAND} -E env bash
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/golden/run_count_table_golden.sh)
    # Architecture A-bar: every ChannelSpec field must have a reader (wire-iff-consumer; no data needed).
    add_test(NAME registry_field_readers
             COMMAND ${CMAKE_COMMAND} -E env bash
                     ${CMAKE_CURRENT_SOURCE_DIR}/tests/registry_field_readers.sh)
endif()
