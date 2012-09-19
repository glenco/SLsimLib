NR_DIR = /Users/cgiocoli/Dropbox/Astro/Codes/EclipseWorkspace/NR
COSMOLIB_DIR = /Users/cgiocoli/Dropbox/Astro/Codes/EclipseWorkspace/CosmoLib

OPT = -DENABLE_FITS

CFITS_DIR = /Users/cgiocoli/Library/CppLib/cfitsio
CCFITS_DIR = /Users/cgiocoli/Library/CppLib/CCfits

INCL = -I./include/ \
	-I$(NR_DIR)/include \
	-I$(COSMOLIB_DIR)/include \
	-I$(CFITS_DIR)/include \
	-I$(CCFITS_DIR)/include

OUT = libSLsimLib.a

SRC = AnalyticNSIE/expanded_powerlaw.cpp \
AnalyticNSIE/lens.cpp \
AnalyticNSIE/lens_expand.cpp \
AnalyticNSIE/nfw_lens.cpp \
AnalyticNSIE/nsie.cpp \
AnalyticNSIE/powerlaw.cpp \
AnalyticNSIE/randomize_lens.cpp \
AnalyticNSIE/readfiles_ana.cpp \
AnalyticNSIE/source.cpp \
BLR/blr_surface_brightness.cpp \
BLR/blr_surface_brightness2.0.cpp \
Fitlens/fitlens.cpp \
FullRange/implant_stars.cpp \
FullRange/internal_rayshooter_multi.cpp \
FullRange/internal_rayshooter_nfw.cpp \
Galaxies/create_sersic.cpp \
Galaxies/overzier_galaxy.cpp \
ImageProcessing/pixelize.cpp \
Kist/Kist.cpp \
MultiPlane/MOKAlens.cpp \
MultiPlane/MOKAfits.cpp \
MultiPlane/model.cpp \
MultiPlane/multiplane.cpp \
MultiPlane/singlelens.cpp \
MultiPlane/profile.cpp \
TreeCode/forceTree.cpp \
TreeCode/profiles.cpp \
TreeCode/qTreeNB.cpp \
TreeCode/quadForceTrees.cpp \
TreeCode/quadTree.cpp \
TreeCode/readfiles.cpp \
TreeCode/simpleTree.cpp \
TreeCode_link/List/List.cpp \
TreeCode_link/KistDriver.cpp \
TreeCode_link/Tree.cpp \
TreeCode_link/TreeDriver.cpp \
TreeCode_link/change_redshifts.cpp \
TreeCode_link/change_redshits.cpp \
TreeCode_link/curve_routines.cpp \
TreeCode_link/dirtycode.cpp \
TreeCode_link/divide_images.cpp \
TreeCode_link/double_sort.cpp \
TreeCode_link/find_crit.cpp \
TreeCode_link/grid_maintenance.cpp \
TreeCode_link/image_finder.cpp \
TreeCode_link/image_finder_kist.cpp \
TreeCode_link/map_images.cpp \
TreeCode_link/peak_refinement.cpp \
TreeCode_link/tree_maintenance.cpp \
TreeCode_link/utilities.cpp \
Utilities/mcmc.cpp 

#
OBJ = $(SRC:.cpp=.o)
# 
OPTIMIZE = -O2
# 
DEBUG = -g 
# compiler  
CC = g++
#
RM = rm -fr
#

CFLAGS = $(INCL) $(DEBUG) $(OPTIMIZE)  -Wno-deprecated -fopenmp
#
CLEAR = clear

.cpp.o:
	$(CC) $(OPT) $(CFLAGS) -c $< -o $@

all: $(OBJ)
	ar rcs $(OUT) $(OBJ)

clean:
	$(RM) $(OBJ) $(OUT) *~

