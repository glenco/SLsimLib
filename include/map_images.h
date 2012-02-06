/*
 * map_images.h
 *
 *  Created on: Oct 6, 2010
 *      Author: bmetcalf
 */
#include "Tree.h"
#include <model.h>

void map_images(ModelHndl model,GridHndl grid,int *Nimages,ImageInfo *imageinfo
		,int Nimagesmax,double initial_size,bool splitimages
		,ExitCriterion criterion,bool kappa_off);
int refine_grid_on_image(ModelHndl model,GridHndl grid,ImageInfo *imageinfo
		,unsigned long Nimages,double res_target,ExitCriterion criterion,bool kappa_off);
