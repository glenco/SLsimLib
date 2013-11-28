cmake_minimum_required(VERSION 2.8)

####
# project settings
####

project(SLsimLib)

set(SLsimLib_VERSION_MAJOR 0)
set(SLsimLib_VERSION_MINOR 1)


####
# macro for easy subdirectory adding
####

macro(add_sources)
	file(RELATIVE_PATH _rel "${PROJECT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
	foreach(_src ${ARGN})
		if(_rel)
			list(APPEND _group "${_rel}/${_src}")
		else()
			list(APPEND _group "${_src}")
		endif()
	endforeach()
	if(_rel)
		list(APPEND SOURCES ${_group})
		set(SOURCES ${SOURCES} PARENT_SCOPE)
		list(APPEND GROUPS ${_rel})
		set(GROUPS ${GROUPS} PARENT_SCOPE)
		set("GROUP_${_rel}" ${_group} PARENT_SCOPE)
	endif()
endmacro()


####
# options
####

set(N_THREADS 1 CACHE STRING "Number of threads for rayshooting.")
add_definitions(-DN_THREADS=${N_THREADS})
message(STATUS "Number of threads: ${N_THREADS}")

####
# sources
####

add_subdirectory(include)
add_subdirectory(AnalyticNSIE)
add_subdirectory(BLR)
add_subdirectory(Fitlens)
add_subdirectory(FullRange)
add_subdirectory(ImageProcessing)
add_subdirectory(InputParameters)
add_subdirectory(MultiPlane)
add_subdirectory(Source)
add_subdirectory(TreeCode)
add_subdirectory(TreeCode_link)


####
# dependencies
####

if(ENABLE_FITS)
	include_directories(${CCFITS_INCLUDE_DIRS})
endif()

if(ENABLE_FFTW)
	include_directories(${FFTW3_INCLUDE_DIRS})
endif()

include_directories(
	"${NR_SOURCE_DIR}/include"
	"${CosmoLib_SOURCE_DIR}/include"
	"include"
)


####
# building
####

add_library(SLsimLib ${HEADERS} ${SOURCES})

add_dependencies(SLsimLib CosmoLib NR)


####
# group sources
####

foreach(_grp ${GROUPS})
	source_group(${_grp} FILES ${GROUP_${_grp}})
endforeach()

