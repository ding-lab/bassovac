cmake_minimum_required(VERSION 2.6)

# .deb packaging
set(ARCH "i686")
if(${CMAKE_SIZEOF_VOID_P} MATCHES 8)
    set(ARCH "x86_64")
endif ()

# The format of the description field is a short summary line followed by a
# longer paragraph indented by a single space on each line
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A tool for detecting somatic mutations in normal/cancer .bam files
 Bassovac (BAyeSian SOmatic VAriants in Cancer) is a tool
 designed by Mike Wendl which analyzes intersections in
 pairs of sorted .bam or .sam files (tumor and normal).
 The output is a set of probability values related to
 whether or not each overlapping site is a somatic
 mutation or not.")

set(CPACK_PACKAGE_NAME "bassovac${EXE_VERSION_SUFFIX}")
set(CPACK_PACKAGE_VENDOR "The Genome Institute")
set(CPACK_PACKAGE_VERSION ${FULL_VERSION})
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Travis Abbott <tabbott@genome.wustl.edu>")
set(CPACK_SYSTEM_NAME "Linux-${ARCH}")
set(CPACK_TOPLEVEL_TAG "Linux-${ARCH}")
set(CPACK_DEBIAN_PACKAGE_SECTION science)
set(CPACK_DEBIAN_PACKAGE_PRIORITY optional)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libboost-program-options1.40.0 (>= 1.40.0-1), libc6 (>= 2.4), libgcc1 (>= 1:4.1.1), libstdc++6 (>= 4.4.0), zlib1g (>= 1:1.2.3.3.dfsg)")
if (CMAKE_BUILD_TYPE MATCHES package)
    set(CPACK_GENERATOR "DEB")
else(CMAKE_BUILD_TYPE MATCHES package)
    set(CPACK_GENERATOR "TGZ")
endif(CMAKE_BUILD_TYPE MATCHES package)

configure_file(debian/postinst.in debian/postinst @ONLY)
configure_file(debian/prerm.in debian/prerm @ONLY)
set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "debian/postinst;debian/prerm")

include(CPack)
