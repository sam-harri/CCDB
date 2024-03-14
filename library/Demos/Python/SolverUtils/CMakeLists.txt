IF (NEKTAR_BUILD_UNIT_TESTS)
    ADD_NEKPY_EXECUTABLE(Filter_UnitTest.py Filter_UnitTest.py)
    ADD_NEKPY_TEST(Filter_UnitTest)
ENDIF()
