# =============================================================================
# copy_runtime_dlls.cmake
#
#   Author : dex / cemfly-april2026
#   License: MIT
# -----------------------------------------------------------------------------
# POST_BUILD action on MinGW/UCRT64. Walks the import table of EXE (via
# `objdump -p`), recursively follows DLLs that live under BINDIR, and
# copies every transitively required runtime DLL next to the executable.
#
# Why objdump and not ldd: once a DLL has been copied next to the
# executable, ldd starts resolving it from there instead of the toolchain
# directory, which breaks the next incremental build's post-step.
# objdump only reports the *declared* imports, so the resolution is
# stable across rebuilds.
#
# Inputs (passed via -D):
#   EXE     : path to the freshly-built binary  (Windows path)
#   OUT     : directory where DLLs should be deposited
#   OBJDUMP : path to objdump.exe from the same toolchain
#   BINDIR  : the toolchain bin/ directory
# =============================================================================

if(NOT EXE OR NOT OUT OR NOT OBJDUMP OR NOT BINDIR)
    message(FATAL_ERROR
        "copy_runtime_dlls.cmake requires EXE, OUT, OBJDUMP, BINDIR")
endif()

# -----------------------------------------------------------------------------
# scan_imports(<binary> <out var>)
#   Calls objdump -p, extracts every "DLL Name: foo.dll" line and stores
#   the lowercase names in <out var> as a CMake list.
# -----------------------------------------------------------------------------
function(scan_imports BINARY OUT_VAR)
    execute_process(
        COMMAND ${OBJDUMP} -p ${BINARY}
        OUTPUT_VARIABLE  RAW
        ERROR_QUIET
        RESULT_VARIABLE  RC
    )
    if(NOT RC EQUAL 0)
        set(${OUT_VAR} "" PARENT_SCOPE)
        return()
    endif()

    set(NAMES "")
    string(REPLACE "\n" ";" LINES "${RAW}")
    foreach(LINE IN LISTS LINES)
        if(LINE MATCHES "^[ \t]*DLL Name:[ \t]+([^ \r\n]+)")
            string(TOLOWER "${CMAKE_MATCH_1}" NAME_LOWER)
            list(APPEND NAMES "${NAME_LOWER}")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES NAMES)
    set(${OUT_VAR} "${NAMES}" PARENT_SCOPE)
endfunction()

# -----------------------------------------------------------------------------
# Iterative BFS over the dependency graph.
#   * Visited DLLs (basenames lowercased) live in VISITED.
#   * Frontier holds DLLs to scan next.
#   * For each DLL whose basename exists under BINDIR/, copy it to OUT and
#     enqueue its imports.
# -----------------------------------------------------------------------------
set(VISITED "")
scan_imports("${EXE}" FRONTIER)

set(COPIED 0)
while(FRONTIER)
    list(POP_FRONT FRONTIER NAME)
    list(FIND VISITED "${NAME}" ALREADY)
    if(NOT ALREADY EQUAL -1)
        continue()
    endif()
    list(APPEND VISITED "${NAME}")

    # Toolchain DLL? Test both lowercase and original-case (objdump preserves
    # case but Windows filesystems are case-insensitive).
    set(SRC "${BINDIR}/${NAME}")
    if(NOT EXISTS "${SRC}")
        continue()
    endif()

    file(COPY "${SRC}" DESTINATION "${OUT}")
    math(EXPR COPIED "${COPIED} + 1")

    scan_imports("${SRC}" SUB)
    foreach(D IN LISTS SUB)
        list(FIND VISITED "${D}" SEEN)
        if(SEEN EQUAL -1)
            list(APPEND FRONTIER "${D}")
        endif()
    endforeach()
endwhile()

message(STATUS "[runtime-deps] Copied ${COPIED} DLLs to ${OUT}")
