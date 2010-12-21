/*
 * map_images.h
 *
 *  Created on: Oct 6, 2010
 *      Author: bmetcalf
 */
#include "Tree.h"

void map_images(AnaLens *lens,TreeHndl s_tree,TreeHndl i_tree,int *Nimages,ImageInfo *imageinfo
		,int Nimagesmax,double initial_size,Boolean splitimages
		,ExitCriterion criterion,Boolean kappa_off);
int refine_grid_on_image(AnaLens *lens,TreeHndl i_tree,TreeHndl s_tree,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,ExitCriterion criterion,Boolean kappa_off);
