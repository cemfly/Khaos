# =============================================================================
# Semiconductor Analysis & Simulation Platform
#
#   Author : dex / cemfly-april2026
#   License: MIT
# -----------------------------------------------------------------------------
# Build dependencies:
#   * SFML 3              -- windowing + 2D rendering
#   * Dear ImGui (docking)-- dockable panels
#   * ImGui-SFML v3       -- ImGui<->SFML 3 binding
#   * ImPlot v0.16        -- live oscilloscope, parameter sweeps, heatmaps
#   * GoogleTest (opt-in) -- physics regression suite
#
# Linux (GCC / Clang):
#     sudo apt install libsfml-dev cmake ninja-build
#     cmake -S . -B build -G Ninja -DBUILD_TESTS=ON
#     cmake --build build -j
#
# Windows (MSYS2 UCRT64 -- the workflow we use):
#     pacman -S mingw-w64-ucrt-x86_64-{gcc,cmake,ninja,sfml}
#     cmake -S . -B build -G Ninja -DBUILD_TESTS=ON
#     cmake --build build -j
# =============================================================================

cmake_minimum_required(VERSION 3.20)

project(SemiconductorSim
        VERSION     4.0.0
        LANGUAGES   CXX
        DESCRIPTION "Semiconductor Analysis & Simulation Platform")

# -----------------------------------------------------------------------------
# Compiler / language
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
    add_compile_options(/W4 /permissive- /utf-8)
endif()

option(FETCH_SFML  "Download SFML if not found"             OFF)
option(BUILD_TESTS "Build the GoogleTest unit-test suite"   OFF)
option(USE_SANITIZER
       "Enable AddressSanitizer + UndefinedBehaviorSanitizer (GCC/Clang only)"
       OFF)

# -----------------------------------------------------------------------------
# Sanitizers (opt-in -- AddressSanitizer + UndefinedBehaviorSanitizer)
# -----------------------------------------------------------------------------
# Toggle with `-DUSE_SANITIZER=ON`.  Adds `-fsanitize=address,undefined` to
# both compile *and* link lines (sanitizers must see every translation unit
# *and* the final link to get usable diagnostics).
#
# Notes
#   * GCC / Clang only.  MSVC ships its own AddressSanitizer with different
#     compile/link flags and a separate runtime; it is intentionally skipped
#     here so this option remains a one-toggle switch.
#   * Pair with `-DCMAKE_BUILD_TYPE=Debug` for line-accurate stack traces.
#   * Set `ASAN_OPTIONS` / `UBSAN_OPTIONS` env vars at run time, e.g.:
#         ASAN_OPTIONS=detect_leaks=1:abort_on_error=1
#         UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1
#   * Sample full configure line:
#         cmake -S . -B build -G Ninja \
#               -DCMAKE_BUILD_TYPE=Debug \
#               -DBUILD_TESTS=ON \
#               -DUSE_SANITIZER=ON
# -----------------------------------------------------------------------------
if(USE_SANITIZER)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "Sanitizers ENABLED: address,undefined "
                       "(compiler=${CMAKE_CXX_COMPILER_ID})")
        add_compile_options(
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
            -fno-sanitize-recover=all
            -g)
        add_link_options(-fsanitize=address,undefined)
    else()
        message(WARNING
            "USE_SANITIZER=ON is GCC/Clang only; ignored for "
            "${CMAKE_CXX_COMPILER_ID}.")
    endif()
endif()

include(FetchContent)

# =============================================================================
# 1) SFML 3
# =============================================================================
find_package(SFML 3 COMPONENTS Graphics Window System QUIET)

if(NOT SFML_FOUND)
    if(FETCH_SFML)
        message(STATUS "SFML not found -- fetching")
        FetchContent_Declare(
            SFML
            GIT_REPOSITORY https://github.com/SFML/SFML.git
            GIT_TAG        3.0.0
        )
        FetchContent_MakeAvailable(SFML)
    else()
        message(FATAL_ERROR
            "SFML 3 not found.  Install via your package manager or "
            "re-run cmake with -DFETCH_SFML=ON")
    endif()
endif()


# =============================================================================
# 2) Dear ImGui (docking branch -- DockSpace + MultiViewport API)
# -----------------------------------------------------------------------------
# We pin to a tagged docking release so behaviour is reproducible across
# machines. The library does not ship a CMakeLists.txt of its own; we let
# imgui-sfml build it (see step 3) by exposing IMGUI_DIR.
# =============================================================================
FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.91.6-docking
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
)
FetchContent_MakeAvailable(imgui)


# =============================================================================
# 3) ImGui-SFML v3 (binding for SFML 3)
# -----------------------------------------------------------------------------
# - IMGUI_DIR        : tells imgui-sfml where the imgui sources live so it
#                      compiles them into the same translation unit set.
# - IMGUI_SFML_FIND_SFML=OFF : we already located SFML above (or fetched it).
# - IMGUI_SFML_IMGUI_DEMO=ON : ImGui::ShowDemoWindow() available -- handy
#                              while the UI is being built up.
# =============================================================================
FetchContent_Declare(
    imgui-sfml
    GIT_REPOSITORY https://github.com/SFML/imgui-sfml.git
    GIT_TAG        v3.0
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
)
set(IMGUI_DIR             ${imgui_SOURCE_DIR} CACHE PATH "" FORCE)
set(IMGUI_SFML_FIND_SFML  OFF                 CACHE BOOL "" FORCE)
set(IMGUI_SFML_IMGUI_DEMO ON                  CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(imgui-sfml)


# =============================================================================
# 4) ImPlot v0.16 (no upstream CMakeLists -- we build a tiny static lib)
# =============================================================================
FetchContent_Declare(
    implot
    GIT_REPOSITORY https://github.com/epezent/implot.git
    GIT_TAG        v0.16
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
)
FetchContent_MakeAvailable(implot)

add_library(implot STATIC
    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp
    ${implot_SOURCE_DIR}/implot_demo.cpp
)
target_include_directories(implot PUBLIC ${implot_SOURCE_DIR})
# implot needs imgui.h + imgui_internal.h; both come transitively via
# the ImGui-SFML interface library.
target_link_libraries(implot PUBLIC ImGui-SFML::ImGui-SFML)


# =============================================================================
# 5) Physics core (object library reused by both app and tests)
# =============================================================================
add_library(semicon_physics OBJECT
    src/PhysicsEngine.cpp
    src/Material.cpp
)
target_include_directories(semicon_physics PUBLIC src)
set_target_properties(semicon_physics PROPERTIES POSITION_INDEPENDENT_CODE ON)


# =============================================================================
# 6) Main application
# -----------------------------------------------------------------------------
# CrystalLattice / BandDiagram / UIPanel still compile so the project keeps
# building during the SFML-UI -> ImGui transition.  They will be retired in
# subsequent phases as their responsibilities migrate to ImGui windows.
# =============================================================================
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/DriftDiffusion.cpp
    src/CrystalView.cpp
    src/BandView.cpp
    $<TARGET_OBJECTS:semicon_physics>
)
target_include_directories(${PROJECT_NAME} PRIVATE src)
target_link_libraries(${PROJECT_NAME} PRIVATE
    SFML::Graphics
    SFML::Window
    SFML::System
    ImGui-SFML::ImGui-SFML
    implot
)

# Copy the assets folder next to the binary (font, colormap, etc.).
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/assets
        $<TARGET_FILE_DIR:${PROJECT_NAME}>/assets
)

# Ship SFML DLLs alongside the binary on Windows when SFML is shared.
if(WIN32 AND TARGET SFML::Graphics)
    add_custom_command(
        TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:SFML::Graphics>
            $<TARGET_FILE:SFML::Window>
            $<TARGET_FILE:SFML::System>
            $<TARGET_FILE_DIR:${PROJECT_NAME}>
        COMMAND_EXPAND_LISTS
    )
endif()

# -----------------------------------------------------------------------------
# MinGW/UCRT64 runtime DLLs.
#
#   On MinGW (MSYS2 UCRT64) the executable depends on a chain of unsigned
#   DLLs (libstdc++-6.dll, libgcc_s_seh-1.dll, libwinpthread-1.dll, plus
#   anything pulled in by SFML such as freetype/png/harfbuzz/zlib/...).
#   When the user double-clicks the exe from Explorer, %PATH% does not
#   include C:\msys64\ucrt64\bin, so loading fails.
#
#   We resolve every required DLL once (via `ldd`) and copy whatever lives
#   under the toolchain's bin/ directory next to the binary. CMake calls
#   the small inline script after every successful link.
# -----------------------------------------------------------------------------
if(WIN32 AND MINGW)
    # objdump from the same toolchain is more reliable than ldd: it lists
    # the binary's declared imports without resolving them through %PATH%.
    get_filename_component(MINGW_BIN_DIR ${CMAKE_CXX_COMPILER} DIRECTORY)
    find_program(OBJDUMP_PROGRAM
        NAMES objdump.exe objdump
        PATHS ${MINGW_BIN_DIR}
        NO_DEFAULT_PATH)

    if(OBJDUMP_PROGRAM)
        add_custom_command(
            TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND}
                -DEXE=$<TARGET_FILE:${PROJECT_NAME}>
                -DOUT=$<TARGET_FILE_DIR:${PROJECT_NAME}>
                -DOBJDUMP=${OBJDUMP_PROGRAM}
                -DBINDIR=${MINGW_BIN_DIR}
                -P ${CMAKE_SOURCE_DIR}/cmake/copy_runtime_dlls.cmake
            COMMENT "Copying MinGW runtime DLLs next to ${PROJECT_NAME}"
        )
    else()
        message(WARNING "objdump not found; skipping MinGW DLL bundling")
    endif()
endif()


# =============================================================================
# 7) Unit tests (GoogleTest via FetchContent)
# =============================================================================
if(BUILD_TESTS)
    message(STATUS "Tests enabled -- fetching GoogleTest")
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

    include(GoogleTest)
    gtest_discover_tests(physics_tests)
endif()
