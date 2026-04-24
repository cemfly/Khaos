# =============================================================================
# Interactive Semiconductor Bandgap & Doping Simulator
# -----------------------------------------------------------------------------
# Phase 2 build:
#   * Main SFML application      (always built)
#   * GoogleTest-based unit test suite  (opt-in with -DBUILD_TESTS=ON)
#
# Linux (GCC / Clang):
#     sudo apt install libsfml-dev        # Debian / Ubuntu
#     cmake -S . -B build -DBUILD_TESTS=ON
#     cmake --build build -j
#     ./build/SemiconductorSim
#     ctest --test-dir build --output-on-failure      # run all tests
#
# Windows (MSVC / MinGW):
#     cmake -S . -B build -DSFML_DIR=<path>/SFML/lib/cmake/SFML
#  or  cmake -S . -B build -DFETCH_SFML=ON            # download SFML
#     cmake --build build --config Release
# =============================================================================

cmake_minimum_required(VERSION 3.16)

project(SemiconductorSim
        VERSION     2.0.0
        LANGUAGES   CXX
        DESCRIPTION "Interactive Semiconductor Bandgap & Doping Simulator")

# -----------------------------------------------------------------------------
# Language / compiler settings
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()
if(MSVC)
    add_compile_options(/W4 /permissive-)
endif()

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------
option(FETCH_SFML  "Download and build SFML with FetchContent if not found" OFF)
option(BUILD_TESTS "Build the GoogleTest unit-test suite"                   OFF)

# =============================================================================
# SFML dependency
#
#   1) find_package() first (system install, vcpkg, Conan, ...).
#   2) If FETCH_SFML=ON, fall back to FetchContent (pulls SFML 2.6.1).
# =============================================================================
find_package(SFML 2.5 COMPONENTS graphics window system QUIET)

if(NOT SFML_FOUND)
    if(FETCH_SFML)
        message(STATUS "SFML not found -- fetching with FetchContent")
        include(FetchContent)
        FetchContent_Declare(
            SFML
            GIT_REPOSITORY https://github.com/SFML/SFML.git
            GIT_TAG        2.6.1
        )
        FetchContent_MakeAvailable(SFML)
    else()
        message(FATAL_ERROR
            "SFML 2.5+ not found.\n"
            "  * Install it with your package manager "
            "(apt/dnf/pacman/brew), or\n"
            "  * pass -DSFML_DIR=<path>/SFML/lib/cmake/SFML, or\n"
            "  * re-run cmake with -DFETCH_SFML=ON")
    endif()
endif()

# =============================================================================
# Shared physics object library
#
# The simulator executable and the unit-test executable both need
# PhysicsEngine.cpp, so we compile it once as an OBJECT library and reuse the
# object files in both targets.  This avoids double compilation and guarantees
# that the tests exercise *exactly* the same object code as the application.
# =============================================================================
add_library(semicon_physics OBJECT
    src/PhysicsEngine.cpp
)
target_include_directories(semicon_physics PUBLIC src)
set_target_properties(semicon_physics PROPERTIES POSITION_INDEPENDENT_CODE ON)

# =============================================================================
# Main application
# =============================================================================
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/CrystalLattice.cpp
    src/BandDiagram.cpp
    src/UIPanel.cpp
    $<TARGET_OBJECTS:semicon_physics>
)

target_include_directories(${PROJECT_NAME} PRIVATE src)

target_link_libraries(${PROJECT_NAME} PRIVATE
    sfml-graphics
    sfml-window
    sfml-system
)

# Copy `assets/` next to the binary so the font loader can find it when the
# executable is run from the build directory.
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/assets
        $<TARGET_FILE_DIR:${PROJECT_NAME}>/assets
)

# On Windows with shared SFML, ship the DLLs next to the binary.
if(WIN32 AND TARGET sfml-graphics)
    add_custom_command(
        TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:sfml-graphics>
            $<TARGET_FILE:sfml-window>
            $<TARGET_FILE:sfml-system>
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
        COMMAND_EXPAND_LISTS
    )
endif()


# =============================================================================
# Unit tests  (GoogleTest via FetchContent)
# -----------------------------------------------------------------------------
# Enabled with -DBUILD_TESTS=ON.  The FetchContent block is only evaluated in
# that case so developers without a network connection can still build the
# application.
# =============================================================================
if(BUILD_TESTS)
    message(STATUS "Unit tests enabled -- fetching GoogleTest")

    include(FetchContent)

    # gtest_force_shared_crt is required on Windows/MSVC so the gtest runtime
    # matches the rest of the project.
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.14.0
    )
    FetchContent_MakeAvailable(googletest)

    enable_testing()

    add_executable(physics_tests
        tests/test_physics.cpp
        $<TARGET_OBJECTS:semicon_physics>
    )
    target_include_directories(physics_tests PRIVATE src)
    target_link_libraries(physics_tests PRIVATE
        GTest::gtest
        GTest::gtest_main
    )

    # Register every TEST() declaration with CTest.
    include(GoogleTest)
    gtest_discover_tests(physics_tests)
endif()
