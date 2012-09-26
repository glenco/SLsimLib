/*
 * image_processing.h
 *
 *  Created on: Feb 28, 2010
 *      Author: R.B. Metcalf
 */

#ifndef IMAGE_PROCESSING_H_
#define IMAGE_PROCESSING_H_

#include <slsimlib.h>
#include <image_info.h>
#include <Tree.h>
#include <point.h>
#include <valarray>

/** \ingroup Image
 * \brief Takes images and pixelizes the flux into regular pixel grid.
 *
 * The routine constructs an internal tree structure of the pixel points for fast calculation of location
 * of each point in the images.
 */

class PixelMap{
public:
	PixelMap(unsigned long Npixels,double range,double *center);
	~PixelMap();

	unsigned long getNpixels(){return Npixels;}
	double getRange(){return range;}
	double getResolution(){return resolution;}

	void Clean();
	void AddImages(ImageInfo *imageinfo,int Nimages,bool constant_sb);
	void AddImages(ImageInfo *imageinfo,int Nimages,double sigma);
	void printASCII();
	void printASCIItoFile(std::string filename);
	void printFITS(std::string filename);
	void smooth(double *map_out,double sigma);

private:
	std::valarray<float> map;
	unsigned long Npixels;
	TreeHndl ptree;
	double resolution,range,center[2];

	void Convert();
	void _SplitFluxIntoPixels(TreeHndl ptree,Branch *leaf,double *leaf_sb);
};

void pixelize(double *map,long Npixels,double range,double *center
		,ImageInfo *imageinfo,int Nimages,bool constant_sb,bool cleanmap
		,bool write_for_skymaker = false, std::string filename="");
void _SplitFluxIntoPixels(TreeHndl ptree,Branch *leaf,double *leaf_sb);
void smoothmap(double *map_out,double *map_in,long Npixels,double range,double sigma);

#endif
