/*
 * divide_images.h
 *
 *  Created on: Nov 11, 2010
 *      Author: bmetcalf
 */

#ifndef DIVIDE_IMAGES_H_
#define DIVIDE_IMAGES_H_

#include "Tree.h"

void find_divide_images(TreeHndl i_tree,TreeHndl s_tree,double *source_x,double source_r
		,ImageInfo *imageinfo,int *Nimages,int Nimagesmax);
int ImageGenus(TreeHndl i_tree,ImageInfo *imageinfo);
void divide_images(TreeHndl i_tree,ImageInfo *imageinfo,int *Nimages,int Nimagesmax);
void partition_images(Point *point,unsigned long *N_in_image,TreeHndl i_tree);
void divide_images_kist(TreeHndl i_tree,ImageInfo *imageinfo,int *Nimages,int Nimagesmax);
double partition_images_kist(Point *point,KistHndl imagekist,TreeHndl i_tree);

#endif
