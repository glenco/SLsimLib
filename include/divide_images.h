/*
 * divide_images.h
 *
 *  Created on: Nov 11, 2010
 *      Author: bmetcalf
 */

void find_divide_images(TreeHndl i_tree,TreeHndl s_tree,double *source_x,double source_r
		,ImageInfo *imageinfo,int *Nimages,int Nimagesmax);
void divide_images(TreeHndl i_tree,ImageInfo *imageinfo,int *Nimages,int Nimagesmax);
void divide_images2(TreeHndl i_tree,ImageInfo *imageinfo,int *Nimages,int Nimagesmax);
void partition_images(Point *point,unsigned long *N_in_image,TreeHndl i_tree);
double partition_images2(Point *point,KistHndl imagekist,TreeHndl i_tree);
