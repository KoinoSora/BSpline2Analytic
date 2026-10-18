# 文件职责：静态库、示例及当前测试目标的构建定义。
# 区间认证源文件使用 /fp:strict；保留历史源文件不表示允许运行其测试，测试开关以仓库约定为准。

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
    src/SequentialConversion.cpp
    src/BSpline2Analytic.cpp src/SurfaceConversion.cpp src/BoundaryAnalysis.cpp src/AnalyticReconstruction.cpp
    src/PatchRecovery.cpp src/RecognitionValidation.cpp src/Numerics.cpp
    src/IntervalBounds.cpp src/ApproximateGeometry.cpp src/ApproximateNumerics.cpp src/IsoCurveAnalysis.cpp
    src/PairedRationalChart.cpp src/PairedConversionValidation.cpp src/PairedSourcePreparation.cpp src/PairedPositionBounds.cpp src/PairedBoundaryValidation.cpp
    src/MotherDistanceValidation.cpp src/MotherDistanceBounds.cpp src/MotherDistanceBoundary.cpp src/MotherSourcePreparation.cpp src/MotherTangentValidation.cpp
    src/SurfaceSampling.cpp src/SurfaceSampleContext.cpp src/ApproximateEvaluation.cpp src/AnalyticJacobian.cpp src/PlaneOffsetRefinement.cpp src/CandidateRefinement.cpp src/ConversionValidation.cpp src/ApproximateConversion.cpp
    src/GeneralClassification.cpp src/GeneralSeedNumerics.cpp src/GeneralConversion.cpp src/GeneralConversionValidation.cpp
    src/GeneralGeometry.cpp src/GeneralPCurve.cpp src/GeneralRefinement.cpp src/GeneralPairedBounds.cpp src/GeneralSphereChart.cpp
    src/GeneralHybridValidation.cpp src/GeneralRationalBoundaryRecovery.cpp
    src/GeneralRationalBoundaryCurves.cpp src/GeneralBoundaryError.cpp src/GeneralBoundaryEnvelope.cpp src/GeneralRationalBoundaryValidation.cpp)
set_source_files_properties(src/GeneralRationalBoundaryRecovery.cpp src/GeneralRationalBoundaryCurves.cpp
    src/GeneralBoundaryError.cpp src/GeneralBoundaryEnvelope.cpp src/GeneralRationalBoundaryValidation.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
set_source_files_properties(src/GeneralHybridValidation.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
set_source_files_properties(src/IntervalBounds.cpp src/ConversionValidation.cpp src/GeneralConversionValidation.cpp src/GeneralPairedBounds.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
set_source_files_properties(src/PairedRationalChart.cpp src/PairedConversionValidation.cpp src/PairedSourcePreparation.cpp src/PairedPositionBounds.cpp src/PairedBoundaryValidation.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
set_source_files_properties(src/MotherDistanceValidation.cpp src/MotherDistanceBounds.cpp src/MotherDistanceBoundary.cpp src/MotherSourcePreparation.cpp src/MotherTangentValidation.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
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
    add_executable(BSpline2AnalyticMotherDistanceTests tests/MotherDistanceTests.cpp
        tests/Dataset500.cpp tests/PairedSurfaceBounds.cpp)
    target_link_libraries(BSpline2AnalyticMotherDistanceTests PRIVATE BSpline2Analytic)
    target_include_directories(BSpline2AnalyticMotherDistanceTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticMotherDistanceTests PRIVATE /W4 /utf-8 /fp:strict)
    add_executable(BSpline2AnalyticMotherDistanceSelectedTests tests/MotherDistanceSelectedTests.cpp
        tests/Dataset500.cpp tests/PairedSurfaceBounds.cpp)
    target_link_libraries(BSpline2AnalyticMotherDistanceSelectedTests PRIVATE BSpline2Analytic)
    target_include_directories(BSpline2AnalyticMotherDistanceSelectedTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticMotherDistanceSelectedTests PRIVATE /W4 /utf-8 /fp:strict)
    target_compile_definitions(BSpline2AnalyticMotherDistanceSelectedTests PRIVATE BSPLINE_MOTHER_STREAMING)
    add_executable(BSpline2AnalyticMotherDistanceDatasetTests tests/MotherDistanceDatasetTests.cpp
        tests/Dataset500.cpp tests/PairedSurfaceBounds.cpp tests/ExtremeDataset.cpp tests/BusinessExtremeDataset.cpp)
    target_link_libraries(BSpline2AnalyticMotherDistanceDatasetTests PRIVATE BSpline2Analytic)
    target_include_directories(BSpline2AnalyticMotherDistanceDatasetTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticMotherDistanceDatasetTests PRIVATE /W4 /utf-8 /fp:strict)
    add_test(NAME BSpline2Analytic.MotherDistanceUnit COMMAND BSpline2AnalyticMotherDistanceTests --unit
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500")
    add_test(NAME BSpline2Analytic.MotherDistanceSphereTorus COMMAND BSpline2AnalyticMotherDistanceTests --sphere-torus
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/mother-sphere-torus-$<CONFIG>.csv")
    set_tests_properties(BSpline2Analytic.MotherDistanceUnit BSpline2Analytic.MotherDistanceSphereTorus PROPERTIES TIMEOUT 120
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    add_library(BSpline2AnalyticRationalFixtures STATIC
        tests/GeneralDataset.cpp tests/GeneralDatasetCertificate.cpp
        tests/GeneralExtremeDataset.cpp tests/GeneralExtremeCertificate.cpp
        tests/GeneralRationalReference.cpp tests/GeneralRationalDataset.cpp
        tests/GeneralRationalDatasetCertificate.cpp)
    target_link_libraries(BSpline2AnalyticRationalFixtures PUBLIC BSpline2Analytic)
    target_include_directories(BSpline2AnalyticRationalFixtures PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticRationalFixtures PRIVATE /W4 /utf-8 /fp:strict)
    add_executable(BSpline2AnalyticGeneralRationalDatasetTests tests/GeneralRationalDatasetTests.cpp)
    add_executable(BSpline2AnalyticSequentialTests tests/SequentialConversionTests.cpp tests/Dataset500.cpp tests/PairedSurfaceBounds.cpp)
    target_link_libraries(BSpline2AnalyticSequentialTests PRIVATE BSpline2AnalyticRationalFixtures)
    target_include_directories(BSpline2AnalyticSequentialTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticSequentialTests PRIVATE /W4 /utf-8 /fp:strict)
    add_executable(BSpline2AnalyticIntroductionCases tests/IntroductionCases.cpp tests/Dataset500.cpp tests/PairedSurfaceBounds.cpp)
    target_link_libraries(BSpline2AnalyticIntroductionCases PRIVATE BSpline2AnalyticRationalFixtures)
    target_include_directories(BSpline2AnalyticIntroductionCases PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticIntroductionCases PRIVATE /W4 /utf-8 /fp:strict)
    add_test(NAME BSpline2Analytic.SequentialShortCircuit COMMAND BSpline2AnalyticSequentialTests
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data" iso 0 1 "${CMAKE_CURRENT_BINARY_DIR}/sequential-unit.csv" unit)
    add_test(NAME BSpline2Analytic.SequentialFallback COMMAND BSpline2AnalyticSequentialTests
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data" main 100 1 "${CMAKE_CURRENT_BINARY_DIR}/sequential-unit.csv" unit)
    set_tests_properties(BSpline2Analytic.SequentialShortCircuit BSpline2Analytic.SequentialFallback PROPERTIES TIMEOUT 30
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    target_link_libraries(BSpline2AnalyticGeneralRationalDatasetTests PRIVATE BSpline2AnalyticRationalFixtures)
    target_include_directories(BSpline2AnalyticGeneralRationalDatasetTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticGeneralRationalDatasetTests PRIVATE /W4 /utf-8 /fp:strict)
    add_executable(BSpline2AnalyticGeneralRationalBoundaryTests tests/GeneralRationalBoundaryTests.cpp)
    set_source_files_properties(tests/GeneralExtremeCertificate.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
    target_link_libraries(BSpline2AnalyticGeneralRationalBoundaryTests PRIVATE BSpline2AnalyticRationalFixtures)
    target_include_directories(BSpline2AnalyticGeneralRationalBoundaryTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticGeneralRationalBoundaryTests PRIVATE /W4 /utf-8 /fp:strict)
    add_executable(BSpline2AnalyticGeneralRationalRecoveryDiagnostics tests/GeneralRationalRecoveryDiagnostics.cpp)
    target_link_libraries(BSpline2AnalyticGeneralRationalRecoveryDiagnostics PRIVATE BSpline2AnalyticRationalFixtures)
    target_include_directories(BSpline2AnalyticGeneralRationalRecoveryDiagnostics PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticGeneralRationalRecoveryDiagnostics PRIVATE /W4 /utf-8 /fp:strict)
    add_test(NAME BSpline2Analytic.GeneralRationalBoundaryUnit COMMAND BSpline2AnalyticGeneralRationalBoundaryTests --rational-unit
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/general-rational-v2")
    set_tests_properties(BSpline2Analytic.GeneralRationalBoundaryUnit PROPERTIES TIMEOUT 120
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    # ValidateGeneralHybrid 已废弃：保留实现、共享工具及历史测试源码，但不构建或注册其测试。
    add_executable(BSpline2AnalyticPairedValidationTests tests/PairedValidationTests.cpp
        tests/Dataset500.cpp tests/PairedSurfaceBounds.cpp tests/ApproximateCases.cpp)
    target_link_libraries(BSpline2AnalyticPairedValidationTests PRIVATE BSpline2Analytic)
    target_include_directories(BSpline2AnalyticPairedValidationTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticPairedValidationTests PRIVATE /W4 /utf-8 /fp:strict)
    add_test(NAME BSpline2Analytic.PairedValidation COMMAND BSpline2AnalyticPairedValidationTests --unit)
    add_test(NAME BSpline2Analytic.PairedSphere COMMAND BSpline2AnalyticPairedValidationTests --sphere
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/paired-sphere-$<CONFIG>.csv")
    add_test(NAME BSpline2Analytic.PairedDataset500 COMMAND BSpline2AnalyticPairedValidationTests --dataset
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/paired-dataset500-$<CONFIG>.csv")
    set_tests_properties(BSpline2Analytic.PairedValidation BSpline2Analytic.PairedSphere BSpline2Analytic.PairedDataset500 PROPERTIES
        TIMEOUT 600 WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    add_executable(BSpline2AnalyticClassificationGateTests tests/ClassificationGateTests.cpp
        tests/Dataset500.cpp tests/PairedSurfaceBounds.cpp)
    target_link_libraries(BSpline2AnalyticClassificationGateTests PRIVATE BSpline2Analytic)
    target_compile_options(BSpline2AnalyticClassificationGateTests PRIVATE /W4 /utf-8)
    target_include_directories(BSpline2AnalyticClassificationGateTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    add_test(NAME BSpline2Analytic.ClassificationGate COMMAND BSpline2AnalyticClassificationGateTests
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/extreme200")
    set_tests_properties(BSpline2Analytic.ClassificationGate PROPERTIES
        TIMEOUT 60 WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    add_executable(BSpline2AnalyticExtremeTests tests/ExtremeDataset.cpp tests/BusinessExtremeDataset.cpp tests/ExtremeTests.cpp
        tests/Dataset500.cpp tests/PairedSurfaceBounds.cpp)
    target_link_libraries(BSpline2AnalyticExtremeTests PRIVATE BSpline2Analytic)
    target_include_directories(BSpline2AnalyticExtremeTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticExtremeTests PRIVATE /W4 /utf-8)
    find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter)
    add_test(NAME BSpline2Analytic.ExtremeData COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/AuditBusinessExtreme.py"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/business-extreme200")
    set_tests_properties(BSpline2Analytic.ExtremeData PROPERTIES TIMEOUT 30)
    add_test(NAME BSpline2Analytic.Extreme200 COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/RunExtremeDataset.py"
        --exe "$<TARGET_FILE:BSpline2AnalyticExtremeTests>"
        --data "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/business-extreme200"
        --output "${CMAKE_CURRENT_BINARY_DIR}/business-extreme200-$<CONFIG>.csv" --timeout 300)
    add_test(NAME BSpline2Analytic.ExtremeModules COMMAND BSpline2AnalyticExtremeTests --modules
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/business-extreme200" "${CMAKE_CURRENT_BINARY_DIR}/business-extreme-modules-$<CONFIG>.csv")
    # 批次有自己的逐例 300 秒进程看门狗；避免 CTest 默认套件超时丢弃部分结果。
    # 进程看门狗覆盖 SDK/审计等工作，与算法内部 50 ms 成功验收预算不是同一个口径。
    set_tests_properties(BSpline2Analytic.Extreme200 PROPERTIES TIMEOUT 61200)
    set_tests_properties(BSpline2Analytic.ExtremeModules PROPERTIES TIMEOUT 60)
    set_tests_properties(BSpline2Analytic.Extreme200 BSpline2Analytic.ExtremeModules PROPERTIES
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    add_executable(BSpline2AnalyticGeneralDatasetTests tests/GeneralDataset.cpp
        tests/GeneralDatasetCertificate.cpp tests/GeneralDatasetTests.cpp)
    target_link_libraries(BSpline2AnalyticGeneralDatasetTests PRIVATE BSpline2Analytic)
    target_include_directories(BSpline2AnalyticGeneralDatasetTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_compile_options(BSpline2AnalyticGeneralDatasetTests PRIVATE /W4 /utf-8)
    set_source_files_properties(tests/GeneralDatasetCertificate.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
    add_test(NAME BSpline2Analytic.General600 COMMAND BSpline2AnalyticGeneralDatasetTests --run
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/general600-$<CONFIG>.csv")
    add_test(NAME BSpline2Analytic.GeneralPositive500 COMMAND BSpline2AnalyticGeneralDatasetTests --run
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/general-positive500-$<CONFIG>.csv" 0 500)
    add_test(NAME BSpline2Analytic.GeneralExtra60 COMMAND BSpline2AnalyticGeneralDatasetTests --extras
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/general-extra60-$<CONFIG>.csv")
    set_tests_properties(BSpline2Analytic.General600 BSpline2Analytic.GeneralExtra60 BSpline2Analytic.GeneralPositive500 PROPERTIES
        WORKING_DIRECTORY "${SGGK_SDK_DIR}"
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${SGGK_SDK_DIR}/x64-win/bin${_sgk_config}")
    add_executable(BSpline2AnalyticDataset500Tests tests/Dataset500.cpp tests/PairedSurfaceBounds.cpp tests/Dataset500Tests.cpp)
    target_link_libraries(BSpline2AnalyticDataset500Tests PRIVATE BSpline2Analytic)
    target_compile_options(BSpline2AnalyticDataset500Tests PRIVATE /W4 /utf-8)
    set_source_files_properties(tests/PairedSurfaceBounds.cpp PROPERTIES COMPILE_OPTIONS "/fp:strict")
    add_test(NAME BSpline2Analytic.Dataset500 COMMAND BSpline2AnalyticDataset500Tests
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500" "${CMAKE_CURRENT_BINARY_DIR}/dataset500-$<CONFIG>.csv")
    add_test(NAME BSpline2Analytic.Dataset500Audit COMMAND BSpline2AnalyticDataset500Tests --audit
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dataset500")
    set_tests_properties(BSpline2Analytic.Dataset500 BSpline2Analytic.Dataset500Audit PROPERTIES
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
    # 白盒回归需要访问私有实现头文件，但这些路径不向库的公共使用者开放。
    foreach(_test_target IN ITEMS
        BSpline2AnalyticDataset500Tests BSpline2AnalyticConeTests BSpline2AnalyticTests
        BSpline2AnalyticApproximateTests BSpline2AnalyticApproximateBenchmark
        BSpline2AnalyticPlateauTests BSpline2AnalyticJacobianTests)
        target_include_directories(${_test_target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
    endforeach()
    # 保留历史实现及测试源码，但不得直接或通过 ConvertApproximate 执行废弃投影验证器。
    set_tests_properties(
        BSpline2Analytic.ClassificationGate BSpline2Analytic.Extreme200
        BSpline2Analytic.ExtremeModules BSpline2Analytic.Dataset500
        BSpline2Analytic.ConeDataset BSpline2Analytic.SpherePolynomialDataset
        BSpline2Analytic.Cone BSpline2Analytic.Approximate
        BSpline2Analytic.Plateau BSpline2Analytic.Jacobian
        PROPERTIES DISABLED TRUE)
endif()

# 已归档的解析 UV 兜底数据仅作历史证据，禁止运行这些测试。
if(BSPLINE2ANALYTIC_BUILD_TESTS)
    set_tests_properties(BSpline2Analytic.General600 BSpline2Analytic.GeneralPositive500 BSpline2Analytic.GeneralExtra60 PROPERTIES DISABLED TRUE)
endif()
