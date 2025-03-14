########################################################################
#
# ThirdParty configuration for Nektar++
#
# HDF5
#
########################################################################

OPTION(NEKTAR_USE_HDF5
    "Enable HDF5 I/O support." OFF)

IF (NEKTAR_USE_HDF5)
    IF (NOT NEKTAR_USE_MPI)
        MESSAGE(FATAL_ERROR "HDF5 requires Nektar++ to be configured with "
                "NEKTAR_USE_MPI for MPI support.")
    ENDIF()

    # Try to find parallel system HDF5 first.
    SET(HDF5_PREFER_PARALLEL ON)
    FIND_PACKAGE(HDF5 QUIET)

    IF (HDF5_FOUND AND NOT HDF5_IS_PARALLEL)
        MESSAGE(STATUS "Non-parallel system HDF5 detected: will build instead.")
        SET(BUILD_HDF5 ON)
    ELSEIF(HDF5_FOUND)
        SET(BUILD_HDF5 OFF)
    ELSE()
        SET(BUILD_HDF5 ON)
    ENDIF()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_HDF5
        "Build HDF5 from ThirdParty" ${BUILD_HDF5}
        "NEKTAR_USE_HDF5" OFF)

    IF(THIRDPARTY_BUILD_HDF5)
        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            hdf5-1.12.3
            PREFIX ${TPSRC}
            URL ${TPURL}/hdf5-1.12.3.tar.bz2
            URL_MD5 5d609bf2a74f980aa42dbe61de452185
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/hdf5-1.12.3
            BINARY_DIR ${TPBUILD}/hdf5-1.12.3 
            TMP_DIR ${TPBUILD}/hdf5-1.12.3-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
		-DDEFAULT_API_VERSION=v110
                -DHDF5_ENABLE_PARALLEL=ON
                -DHDF5_BUILD_CPP_LIB=OFF
                -DBUILD_TESTING=OFF
                -DHDF5_BUILD_TOOLS=OFF
                ${TPSRC}/hdf5-1.12.3
            )

        THIRDPARTY_LIBRARY(HDF5_LIBRARIES SHARED hdf5
            DESCRIPTION "HDF5 library")
        SET(HDF5_INCLUDE_DIRS ${TPDIST}/include CACHE FILEPATH
            "HDF5 include directory" FORCE)

        MESSAGE(STATUS "Build HDF5: ${HDF5_LIBRARIES}")

        SET(HDF5_CONFIG_INCLUDE_DIR ${TPINC})
    ELSE()
        MESSAGE(STATUS "Found HDF5: ${HDF5_LIBRARIES}")
        SET(HDF5_CONFIG_INCLUDE_DIR ${HDF5_INCLUDE_DIRS})
        ADD_CUSTOM_TARGET(hdf5-1.12.3 ALL)

        # Newer HDF5 versions have changed the API
        # We compile deprecated symbols using the old API if newer than 1.10.0
        IF(HDF5_VERSION VERSION_GREATER_EQUAL 1.10.0)
            ADD_DEFINITIONS(-DH5_USE_110_API)
        ENDIF()
    ENDIF()

    ADD_DEPENDENCIES(thirdparty hdf5-1.12.3)

    MARK_AS_ADVANCED(HDF5_LIBRARIES)
    MARK_AS_ADVANCED(HDF5_INCLUDE_DIRS)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
ENDIF()
