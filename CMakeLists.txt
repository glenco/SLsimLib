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
# headers
####

set(HEADERS
	include/analytic_lens.h
	include/base_analens.h
	include/causticdata.h
	include/divide_images.h
	include/doc.h
	include/fitlens.h
	include/forceTree.h
	include/grid_maintenance.h
	include/image_info.h
	include/image_processing.h
	include/InputParams.h
	include/Kist.h
	include/KistDriver.h
	include/lens.h
	include/lens_halos.h
	include/List1.h
	include/List_old.h
	include/map_images.h
	include/mcmc.h
	include/model.h
	include/MOKAlens.h
	include/nestedsampler.h
	include/nsie.h
	include/overzier_source.h
	include/parameters.h
	include/peak_refinement.h
	include/planes.h
	include/point.h
	include/pointlist.h
	include/profile.h
	include/quadTree.h
	include/sersic_source.h
	include/simpleTree.h
	include/sky.h
	include/slsimlib.h
	include/source.h
	include/source_models.h
	include/sourceAnaGalaxy.h
	include/standard.h
	include/Tree.h
	include/tree_maintenance.h
	include/TreeNB.h
	include/uniform_lens.h
	include/utilities_slsim.h
)


####
# sources
####

add_subdirectory(AnalyticNSIE)
add_subdirectory(BLR)
add_subdirectory(Fitlens)
add_subdirectory(FullRange)
add_subdirectory(ImageProcessing)
add_subdirectory(InputParameters)
add_subdirectory(MCMC)
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

source_group("include" FILES ${HEADERS})

foreach(_grp ${GROUPS})
	source_group(${_grp} FILES ${GROUP_${_grp}})
endforeach()

