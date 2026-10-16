cmake_minimum_required(VERSION 3.22)
project(GeneralExtreme320 LANGUAGES CXX)
enable_testing()
# An isolated entry: no existing test target or dataset is executed here.
set(BSPLINE2ANALYTIC_BUILD_TESTS OFF CACHE BOOL "Build the existing suites" FORCE)
add_subdirectory(../.. library)
add_executable(BSpline2AnalyticGeneralExtremeTests
    ../GeneralExtremeDataset.cpp ../GeneralExtremeCertificate.cpp ../GeneralExtremeTests.cpp
    ../GeneralDataset.cpp ../GeneralDatasetCertificate.cpp)
target_link_libraries(BSpline2AnalyticGeneralExtremeTests PRIVATE BSpline2Analytic)
target_include_directories(BSpline2AnalyticGeneralExtremeTests PRIVATE ../../src ..)
target_compile_options(BSpline2AnalyticGeneralExtremeTests PRIVATE /W4 /utf-8)
set_source_files_properties(../GeneralExtremeCertificate.cpp ../GeneralDatasetCertificate.cpp
    PROPERTIES COMPILE_OPTIONS "/fp:strict")
find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
add_test(NAME BSpline2Analytic.GeneralExtreme320 COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../RunGeneralExtreme.py"
    --exe "$<TARGET_FILE:BSpline2AnalyticGeneralExtremeTests>"
    --sdk "${SGGK_SDK_DIR}"
    --configuration "$<CONFIG>"
    --data "${CMAKE_CURRENT_SOURCE_DIR}/../data/general-extreme320"
    --output "${CMAKE_CURRENT_BINARY_DIR}/results.csv")
# Let the per-input watchdog finish the matrix and retain every result row.
set_tests_properties(BSpline2Analytic.GeneralExtreme320 PROPERTIES TIMEOUT 11000)

set_tests_properties(BSpline2Analytic.GeneralExtreme320 PROPERTIES DISABLED TRUE)
