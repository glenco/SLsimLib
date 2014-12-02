/*
 * divide_images.h
 *
 *  Created on: Nov 11, 2010
 *      Author: bmetcalf
 */

#ifndef DIVIDE_IMAGES_H_
#define DIVIDE_IMAGES_H_

#include "Tree.h"

void find_divide_images(TreeHndl i_tree,TreeHndl s_tree,PosType *source_x,PosType source_r
		,ImageInfo *imageinfo,int *Nimages,int Nimagesmax);
int ImageGenus(TreeHndl i_tree,ImageInfo *imageinfo);
void divide_images(TreeHndl i_tree,ImageInfo *imageinfo,int *Nimages,int Nimagesmax);
void partition_images(Point *point,unsigned long *N_in_image,TreeHndl i_tree);
void divide_images_kist(TreeHndl i_tree,std::vector<ImageInfo> &imageinfo,int *Nimages);
PosType partition_images_kist(Point *point,Kist<Point> * imagekist,TreeHndl i_tree);

#endif
