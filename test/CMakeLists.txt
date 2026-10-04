cmake_minimum_required(VERSION 3.14)

# Allow this directory to be configured by itself.  This keeps the unit-test
# workflow independent of SLsimLib's optional external dependencies.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
	project(SLsimLibTests LANGUAGES CXX)
	include(CTest)
endif()

if(NOT BUILD_TESTING)
	return()
endif()

if(NOT TARGET GTest::gtest_main)
	include(FetchContent)
	FetchContent_Declare(
		googletest
		URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.tar.gz
		DOWNLOAD_EXTRACT_TIMESTAMP TRUE
	)
	set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
	set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
	FetchContent_MakeAvailable(googletest)
endif()

add_executable(
	slsimlib_example_test
	example_test.cpp
)

target_compile_features(slsimlib_example_test PRIVATE cxx_std_11)
target_include_directories(
	slsimlib_example_test
	PRIVATE "${CMAKE_CURRENT_LIST_DIR}/../include"
)
target_link_libraries(slsimlib_example_test PRIVATE GTest::gtest_main)

if(TARGET SLsimLib)
	add_executable(
		slsimlib_lenshalohalos_test
		test_lenshalohalos.cpp
	)
	target_compile_features(slsimlib_lenshalohalos_test PRIVATE cxx_std_14)
	target_link_libraries(
		slsimlib_lenshalohalos_test
		PRIVATE
			SLsimLib
			CosmoLib
			NR
			GTest::gtest_main
			${CMAKE_THREAD_LIBS_INIT}
			${CFITSIO_LIBRARIES}
			${FFTW3_LIBRARIES}
			${GSL_LIBRARIES}
			${HEALPIX_LIBRARIES}
			${HDF5_LIBRARIES}
			${CERF_LIBRARIES}
	)
endif()

include(GoogleTest)
gtest_discover_tests(slsimlib_example_test)

if(TARGET slsimlib_lenshalohalos_test)
	gtest_discover_tests(slsimlib_lenshalohalos_test)
endif()
