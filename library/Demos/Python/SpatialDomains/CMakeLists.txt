ADD_NEKPY_EXECUTABLE(CreateMesh.py CreateMesh.py)

IF (NEKTAR_BUILD_TESTS)
    ADD_NEKPY_TEST(CreateMesh)
ENDIF()
