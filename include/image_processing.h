/*
 * image_processing.h
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */


void pixelize(double *map,long Npixels,double range,double *center
		,ImageInfo *imageinfo,int Nimages,Boolean constant_sb,Boolean cleanmap);
void _SplitFluxIntoPixels(TreeHndl ptree,Branch *leaf,double *leaf_sb);
void smoothmap(double *map_out,double *map_in,long Npixels,double range,double sigma);
