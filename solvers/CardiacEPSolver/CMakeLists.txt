CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_CARDIAC_EP
    "Build the Cardiac electrophysiology solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_CARDIAC_EP OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_CARDIAC_EP )
    SET(CardiacEPSolverSources
        EquationSystems/Monodomain.cpp
        EquationSystems/Bidomain.cpp
        EquationSystems/BidomainRoth.cpp
        CellModels/CellModel.cpp
        CellModels/FitzhughNagumo.cpp
        CellModels/AlievPanfilov.cpp
        CellModels/CourtemancheRamirezNattel98.cpp
        CellModels/FentonKarma.cpp
        CellModels/LuoRudy91.cpp
        CellModels/Fox02.cpp
        CellModels/PanditGilesDemir03.cpp
        CellModels/Winslow99.cpp
        CellModels/TenTusscher06.cpp
        Filters/FilterCheckpointCellModel.cpp
        Filters/FilterElectrogram.cpp
        Filters/FilterBenchmark.cpp
        Filters/FilterCellHistoryPoints.cpp
        Filters/FilterOffsetPhase.cpp
        Stimuli/Stimulus.cpp
        Stimuli/StimulusCircle.cpp
        Stimuli/StimulusRect.cpp
        Stimuli/StimulusPoint.cpp
        Stimuli/Protocol.cpp
        Stimuli/ProtocolSingle.cpp
        Stimuli/ProtocolS1.cpp
        Stimuli/ProtocolS1S2.cpp
        )

    IF( NEKTAR_USE_FFTW )
        SET(CardiacEPSolverSources ${CardiacEPSolverSources}
            Filters/FilterHilbertFFTPhase.cpp
        )
    ENDIF( NEKTAR_USE_FFTW )

    ADD_SOLVER_EXECUTABLE(CardiacEPSolver
        SOURCES CardiacEPSolver.cpp
        LIBRARY_SOURCES
        ${CardiacEPSolverSources})

    IF(NEKTAR_BUILD_SOLVER_LIBS)
        SET(CardiacEPSolverHeaders
            CellModels/AlievPanfilov.h
            CellModels/CellModel.h
            CellModels/CourtemancheRamirezNattel98.h
            CellModels/FentonKarma.h
            CellModels/FitzhughNagumo.h
            CellModels/Fox02.h
            CellModels/LuoRudy91.h
            CellModels/PanditGilesDemir03.h
            CellModels/TenTusscher06.h
            CellModels/Winslow99.h
            EquationSystems/Bidomain.h
            EquationSystems/BidomainRoth.h
            EquationSystems/Monodomain.h
            Filters/FilterBenchmark.h
            Filters/FilterCellHistoryPoints.h
            Filters/FilterCheckpointCellModel.h
            Filters/FilterElectrogram.h
            Filters/FilterHilbertFFTPhase.h
            Filters/FilterOffsetPhase.h
            Stimuli/Protocol.h
            Stimuli/ProtocolS1.h
            Stimuli/ProtocolS1S2.h
            Stimuli/ProtocolSingle.h
            Stimuli/StimulusCircle.h
            Stimuli/Stimulus.h
            Stimuli/StimulusPoint.h
            Stimuli/StimulusRect.h)
        # Create separate library. No need to pass source files, as the executable's
        # LIBRARY_SOURCES objects will be reused
        ADD_SOLVER_LIBRARY(CardiacEPSolver HEADERS ${CardiacEPSolverHeaders}
                            DEPENDS SolverUtils)
    ENDIF()

    ADD_SUBDIRECTORY(Utilities)
ENDIF( NEKTAR_SOLVER_CARDIAC_EP )
