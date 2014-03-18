# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/meiers/apps/benchmark_radixSort/CMakeLists.txt
#
# CMakeLists.txt file for benchmark_radixSort.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_meiers_apps_benchmark_radixSort)
message (STATUS "Configuring sandbox/meiers/apps/benchmark_radixSort")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (benchmark_radixSort benchmark_radixSort.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (benchmark_radixSort ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (benchmark_radixSort)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install benchmark_radixSort in ${PREFIX}/bin directory
install (TARGETS benchmark_radixSort
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/benchmark_radixSort for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (benchmark_radixSort)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:benchmark_radixSort")
  set (CPACK_PACKAGE_NAME "benchmark_radixSort")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "benchmark_radixSort")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Sascha Meiers <meiers@inf.fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Sascha Meiers <meiers@inf.fu-berlin.de>")

  seqan_configure_cpack_app (benchmark_radixSort "benchmark_radixSort")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:benchmark_radixSort")

