SET(MULTI_REGIONS_SOURCES
./AssemblyMap/AssemblyCommDG.cpp
./AssemblyMap/AssemblyMap.cpp
./AssemblyMap/AssemblyMapCG.cpp
./AssemblyMap/AssemblyMapDG.cpp
./AssemblyMap/InterfaceMapDG.cpp
./AssemblyMap/LocTraceToTraceMap.cpp
ContField.cpp
ContField3DHomogeneous1D.cpp
ContField3DHomogeneous2D.cpp
DisContField.cpp
DisContField3DHomogeneous1D.cpp
DisContField3DHomogeneous2D.cpp
ExpList.cpp
ExpListHomogeneous1D.cpp
ExpListHomogeneous2D.cpp
ExpList2DHomogeneous1D.cpp
ExpList2DHomogeneous2D.cpp
ExpList3DHomogeneous1D.cpp
ExpList3DHomogeneous2D.cpp
GJPStabilisation.cpp
GlobalLinSys.cpp
GlobalLinSysKey.cpp
GlobalLinSysDirect.cpp
GlobalLinSysDirectFull.cpp
GlobalLinSysDirectStaticCond.cpp
GlobalLinSysIterative.cpp
GlobalLinSysIterativeFull.cpp
GlobalLinSysIterativeStaticCond.cpp
GlobalLinSysStaticCond.cpp
GlobalMatrix.cpp
GlobalMatrixKey.cpp
Preconditioner.cpp
PreconditionerDiagonal.cpp
PreconditionerLowEnergy.cpp
PreconditionerBlock.cpp
SubStructuredGraph.cpp
)

SET(MULTI_REGIONS_HEADERS
ContField.h
ContField3DHomogeneous1D.h
ContField3DHomogeneous2D.h
DisContField.h
DisContField3DHomogeneous1D.h
DisContField3DHomogeneous2D.h
ExpList.h
ExpListHomogeneous1D.h
ExpListHomogeneous2D.h
ExpList2DHomogeneous1D.h
ExpList2DHomogeneous2D.h
ExpList3DHomogeneous1D.h
ExpList3DHomogeneous2D.h
GJPStabilisation.h
GlobalLinSys.h
GlobalLinSysKey.h
GlobalLinSysDirect.h
GlobalLinSysDirectFull.h
GlobalLinSysDirectStaticCond.h
GlobalLinSysIterative.h
GlobalLinSysIterativeFull.h
GlobalLinSysIterativeStaticCond.h
GlobalLinSysStaticCond.h
GlobalMatrix.h
GlobalMatrixKey.h
MultiRegions.hpp
MultiRegionsDeclspec.h
Preconditioner.h
PreconditionerDiagonal.h
PreconditionerLowEnergy.h
PreconditionerBlock.h
SubStructuredGraph.h
)

SET(ASSEMBLY_MAP_HEADERS
./AssemblyMap/AssemblyCommDG.h
./AssemblyMap/AssemblyMap.h
./AssemblyMap/AssemblyMapCG.h
./AssemblyMap/AssemblyMapDG.h
./AssemblyMap/InterfaceMapDG.h
./AssemblyMap/LocTraceToTraceMap.h
)

IF(NEKTAR_USE_MPI)
    SET(MULTI_REGIONS_HEADERS ${MULTI_REGIONS_HEADERS}
        GlobalLinSysXxt.h
        GlobalLinSysXxtFull.h
        GlobalLinSysXxtStaticCond.h
        PreconditionerLinear.h
        PreconditionerLinearWithLowEnergy.h
        PreconditionerLinearWithDiag.h
        PreconditionerLinearWithBlock.h
    )
    SET(MULTI_REGIONS_SOURCES ${MULTI_REGIONS_SOURCES}
        GlobalLinSysXxt.cpp
        GlobalLinSysXxtFull.cpp
        GlobalLinSysXxtStaticCond.cpp
        PreconditionerLinear.cpp
        PreconditionerLinearWithLowEnergy.cpp
        PreconditionerLinearWithDiag.cpp
        PreconditionerLinearWithBlock.cpp
    )
ENDIF(NEKTAR_USE_MPI)

IF(NEKTAR_USE_PETSC)
    SET(MULTI_REGIONS_HEADERS ${MULTI_REGIONS_HEADERS}
        GlobalLinSysPETSc.h
        GlobalLinSysPETScFull.h
        GlobalLinSysPETScStaticCond.h
    )
    SET(MULTI_REGIONS_SOURCES ${MULTI_REGIONS_SOURCES}
        GlobalLinSysPETSc.cpp
        GlobalLinSysPETScFull.cpp
        GlobalLinSysPETScStaticCond.cpp
    )
ENDIF(NEKTAR_USE_PETSC)


ADD_DEFINITIONS(-DMULTI_REGIONS_EXPORTS)

ADD_NEKTAR_LIBRARY(MultiRegions
    SOURCES ${MULTI_REGIONS_SOURCES}
    HEADERS ${MULTI_REGIONS_HEADERS}
    DEPENDS Collections
    SUMMARY "Nektar++ MultiRegions library"
    DESCRIPTION "This library provides global expansions on multi-element domains.")

# MultiRegions uses additional routines from GSMPI
IF( NEKTAR_USE_MPI )
    TARGET_LINK_LIBRARIES(MultiRegions LINK_PRIVATE ${GSMPI_LIBRARY} ${XXT_LIBRARY})
ENDIF( NEKTAR_USE_MPI )

IF( NEKTAR_USE_PETSC )
    TARGET_LINK_LIBRARIES(MultiRegions LINK_PRIVATE ${PETSC_LIBRARIES})
    ADD_DEPENDENCIES(MultiRegions petsc-3.19.3)
ENDIF( NEKTAR_USE_PETSC )

IF (NEKTAR_BUILD_PYTHON)
    SUBDIRS(Python)
ENDIF()

INSTALL(DIRECTORY ./ DESTINATION ${NEKTAR_INCLUDE_DIR}/MultiRegions COMPONENT dev FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")



