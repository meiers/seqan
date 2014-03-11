# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /core/tests/inplace_radix_sort/CMakeLists.txt
#
# CMakeLists.txt file for the inplace_radix_sort module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_core_tests_inplace_radix_sort)
message (STATUS "Configuring core/tests/inplace_radix_sort")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_inplace_radix_sort
                test_inplace_radix_sort.cpp
                test_inplace_radix_sort.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_inplace_radix_sort ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_inplace_radix_sort COMMAND $<TARGET_FILE:test_inplace_radix_sort>)
