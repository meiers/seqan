# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/meiers/apps/benchmarkModCyclicShapeIterator/CMakeLists.txt
#
# CMakeLists.txt file for benchmarkModCyclicShapeIterator.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_meiers_apps_benchmarkModCyclicShapeIterator)
message (STATUS "Configuring sandbox/meiers/apps/benchmarkModCyclicShapeIterator")

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
add_executable (benchmarkModCyclicShapeIterator benchmarkModCyclicShapeIterator.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (benchmarkModCyclicShapeIterator ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (benchmarkModCyclicShapeIterator)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install benchmarkModCyclicShapeIterator in ${PREFIX}/bin directory
install (TARGETS benchmarkModCyclicShapeIterator
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/benchmarkModCyclicShapeIterator for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (benchmarkModCyclicShapeIterator)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:benchmarkModCyclicShapeIterator")
  set (CPACK_PACKAGE_NAME "benchmarkModCyclicShapeIterator")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "benchmarkModCyclicShapeIterator")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Sascha Meiers <meiers@inf.fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Sascha Meiers <meiers@inf.fu-berlin.de>")

  seqan_configure_cpack_app (benchmarkModCyclicShapeIterator "benchmarkModCyclicShapeIterator")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:benchmarkModCyclicShapeIterator")

