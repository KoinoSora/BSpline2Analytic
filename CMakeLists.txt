cmake_minimum_required(VERSION 3.22)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(BSpline2AnalyticProject LANGUAGES CXX)
    include(CTest)
endif()
set(SGGK_SDK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../SGGK" CACHE PATH "SGK 1.7.11 SDK root (include and x64-win)")
get_filename_component(SGGK_SDK_DIR "${SGGK_SDK_DIR}" ABSOLUTE)
if(NOT EXISTS "${SGGK_SDK_DIR}/include/Geometry/3D/Surface/BSplineSurface.h")
    message(FATAL_ERROR "Set SGGK_SDK_DIR to your locally licensed SGK SDK root")
endif()
if(NOT WIN32 OR NOT MSVC)
    message(FATAL_ERROR "This SDK integration currently targets Windows x64 with MSVC")
endif()
add_library(BSpline2Analytic STATIC
    src/BSpline2Analytic.cpp src/SurfaceConversion.cpp src/BoundaryAnalysis.cpp src/AnalyticReconstruction.cpp
    src/PatchRecovery.cpp src/RecognitionValidation.cpp src/Numerics.cpp
    src/IntervalBounds.cpp src/ApproximateGeometry.cpp src/ApproximateNumerics.cpp src/IsoCurveAnalysis.cpp
    src/SurfaceSampling.cpp src/SurfaceSampleContext.cpp src/ApproximateEvaluation.cpp src/AnalyticJacobian.cpp src/PlaneOffsetRefinement.cpp src/CandidateRefinement.cpp src/ConversionValidation.cpp src/ApproximateConversion.cpp
    src/GeneralClassification.cpp src/GeneralConversion.cpp src/GeneralConversionValidation.cpp
    src/GeneralGeometry.cpp src/GeneralPCurve.cpp src/GeneralRefinement.cpp)
set_source_files_properties(src/IntervalBounds.cpp src/ConversionValidation.cpp src/GeneralConversionValidation.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
target_compile_features(BSpline2Analytic PUBLIC cxx_std_17)
target_compile_options(BSpline2Analytic PRIVATE /W4 /utf-8 /permissive-)
target_compile_definitions(BSpline2Analytic PUBLIC NOMINMAX)
target_include_directories(BSpline2Analytic
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_include_directories(BSpline2Analytic SYSTEM PUBLIC "${SGGK_SDK_DIR}/include")
set(_sgk_config "$<IF:$<CONFIG:Debug>,d,$<IF:$<CONFIG:RelWithDebInfo>,w,>>")
target_link_directories(BSpline2Analytic PUBLIC "${SGGK_SDK_DIR}/x64-win/lib${_sgk_config}")
target_link_libraries(BSpline2Analytic PUBLIC SGGK_Foundation SGGK_Math SGGK_GeomBase SGGK_Geometry SGGK_BSplineAlgo)
add_executable(BSpline2AnalyticSample examples/main.cpp)
target_link_libraries(BSpline2AnalyticSample PRIVATE BSpline2Analytic)
option(BSPLINE2ANALYTIC_BUILD_TESTS "Build module regression tests" ON)
if(BSPLINE2ANALYTIC_BUILD_TESTS)
    enable_testing()
    add_executable(BSpline2AnalyticGeneralDatasetTests tests/GeneralDataset.cpp
        tests/GeneralDatasetCertificate.cpp tests/GeneralDatasetTests.cpp)
    target_link_libraries(BSpline2AnalyticGeneralDatasetTests PRIVATE BSpline2Analytic)
    target_include_directories(BSpline2AnalyticGeneralDatasetTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticGeneralDatasetTests PRIVATE /W4 /utf-8)
    set_source_files_properties(tests/GeneralDatasetCertificate.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
    add_test(NAME BSpline2Analytic.General600 COMMAND BSpline2AnalyticGeneralDatasetTests --run
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/general600-$<CONFIG>.csv")
    add_test(NAME BSpline2Analytic.GeneralExtra60 COMMAND BSpline2AnalyticGeneralDatasetTests --extras
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/general-extra60-$<CONFIG>.csv")
    set_tests_properties(BSpline2Analytic.General600 BSpline2Analytic.GeneralExtra60 PROPERTIES
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    add_executable(BSpline2AnalyticDataset500Tests tests/Dataset500.cpp tests/PairedSurfaceBounds.cpp tests/Dataset500Tests.cpp)
    target_link_libraries(BSpline2AnalyticDataset500Tests PRIVATE BSpline2Analytic)
    target_compile_options(BSpline2AnalyticDataset500Tests PRIVATE /W4 /utf-8)
    set_source_files_properties(tests/PairedSurfaceBounds.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
    add_test(NAME BSpline2Analytic.Dataset500 COMMAND BSpline2AnalyticDataset500Tests
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/dataset500-$<CONFIG>.csv")
    set_tests_properties(BSpline2Analytic.Dataset500 PROPERTIES
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    add_test(NAME BSpline2Analytic.ConeDataset COMMAND BSpline2AnalyticDataset500Tests --cone-only
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/dataset500-cone-$<CONFIG>.csv")
    add_test(NAME BSpline2Analytic.SpherePolynomialDataset COMMAND BSpline2AnalyticDataset500Tests --sphere-polynomial-only
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/dataset500-sphere-polynomial-$<CONFIG>.csv")
    set_tests_properties(BSpline2Analytic.SpherePolynomialDataset PROPERTIES
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    add_executable(BSpline2AnalyticConeTests tests/ConeTests.cpp tests/Dataset500.cpp
        tests/PairedSurfaceBounds.cpp tests/ApproximateCases.cpp)
    target_link_libraries(BSpline2AnalyticConeTests PRIVATE BSpline2Analytic)
    target_compile_options(BSpline2AnalyticConeTests PRIVATE /W4 /utf-8)
    add_test(NAME BSpline2Analytic.Cone COMMAND BSpline2AnalyticConeTests
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/cone249-$<CONFIG>.csv")
    set_tests_properties(BSpline2Analytic.ConeDataset BSpline2Analytic.Cone PROPERTIES
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    add_executable(BSpline2AnalyticTests tests/RecognitionTests.cpp tests/ReferenceCases.cpp)
    target_link_libraries(BSpline2AnalyticTests PRIVATE BSpline2Analytic)
    target_compile_options(BSpline2AnalyticTests PRIVATE /W4 /utf-8)
    add_test(NAME BSpline2Analytic.Regression COMMAND BSpline2AnalyticTests)
    add_executable(BSpline2AnalyticApproximateTests tests/ApproximateConversionTests.cpp tests/ApproximateCases.cpp tests/ReferenceCases.cpp)
    target_link_libraries(BSpline2AnalyticApproximateTests PRIVATE BSpline2Analytic)
    target_compile_options(BSpline2AnalyticApproximateTests PRIVATE /W4 /utf-8)
    add_test(NAME BSpline2Analytic.Approximate COMMAND BSpline2AnalyticApproximateTests)
    add_executable(BSpline2AnalyticApproximateBenchmark tests/ApproximateBenchmark.cpp tests/ApproximateCases.cpp tests/ReferenceCases.cpp)
    target_link_libraries(BSpline2AnalyticApproximateBenchmark PRIVATE BSpline2Analytic)
    target_compile_options(BSpline2AnalyticApproximateBenchmark PRIVATE /W4 /utf-8)
    add_executable(BSpline2AnalyticPlateauTests tests/PlateauTests.cpp tests/ApproximateCases.cpp tests/ReferenceCases.cpp)
    target_link_libraries(BSpline2AnalyticPlateauTests PRIVATE BSpline2Analytic)
    target_compile_options(BSpline2AnalyticPlateauTests PRIVATE /W4 /utf-8)
    add_test(NAME BSpline2Analytic.Plateau COMMAND BSpline2AnalyticPlateauTests "${CMAKE_CURRENT_BINARY_DIR}/plateau-$<CONFIG>.csv")
    add_executable(BSpline2AnalyticJacobianTests tests/JacobianTests.cpp tests/ApproximateCases.cpp tests/ReferenceCases.cpp)
    target_link_libraries(BSpline2AnalyticJacobianTests PRIVATE BSpline2Analytic)
    target_compile_options(BSpline2AnalyticJacobianTests PRIVATE /W4 /utf-8)
    add_test(NAME BSpline2Analytic.Jacobian COMMAND BSpline2AnalyticJacobianTests "${CMAKE_CURRENT_BINARY_DIR}/perturbation-$<CONFIG>.csv")
    set_tests_properties(BSpline2Analytic.Jacobian PROPERTIES
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    set_tests_properties(BSpline2Analytic.Plateau PROPERTIES
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    set_tests_properties(BSpline2Analytic.Approximate PROPERTIES
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    set_tests_properties(BSpline2Analytic.Regression PROPERTIES
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    # White-box regression helpers use private implementation headers.
    foreach(_test_target IN ITEMS
        BSpline2AnalyticDataset500Tests BSpline2AnalyticConeTests BSpline2AnalyticTests
        BSpline2AnalyticApproximateTests BSpline2AnalyticApproximateBenchmark
        BSpline2AnalyticPlateauTests BSpline2AnalyticJacobianTests)
        target_include_directories(${_test_target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    endforeach()
endif()
