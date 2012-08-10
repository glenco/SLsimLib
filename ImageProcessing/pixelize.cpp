/*
 * pixelize.c
 *
 *  Created on: Feb 27, 2010
 *      Author: R.B. Metcalf
 */

#include <slsimlib.h>


PixelMap::PixelMap(
		unsigned long my_Npixels  /// Number of pixels in one dimension of map.
		,double my_range          /// One dimensional range of map in whatever units the point positions are in
		,double *my_center        /// The location of the center of the map
		): Npixels(my_Npixels), range(my_range)
		{

	if( (Npixels & (Npixels-1)) != 0){
		ERROR_MESSAGE();
		std::printf("ERROR: pixelsize, Npixels is not a power of 2\n");
		exit(1);
	}

	center[0] = my_center[0];
	center[1] = my_center[1];

	resolution=range/(Npixels-1);

	// initialize pixel tree
	Point *points=NewPointArray(Npixels*Npixels,true);
	xygridpoints(points,range,center,Npixels,false);
	ptree=BuildTree(points,Npixels*Npixels);

	MoveToTopList(ptree->pointlist);
	for(unsigned long i=0 ; i < ptree->pointlist->Npoints ; ++i){
		ptree->pointlist->current->surface_brightness = 0.0;
		MoveDownList(ptree->pointlist);
	}

	map = NULL;

	return;
}
PixelMap::~PixelMap(){

	freeTree(ptree);
	if(map != NULL) delete[] map;
}
/// Zero the whole map
void PixelMap::Clean(){
	unsigned long i;

	if(map != NULL) for(i=0;i < Npixels*Npixels; ++i) map[i] = 0.0;

	MoveToTopList(ptree->pointlist);
	for(i=0 ; i < ptree->pointlist->Npoints ; ++i){
		ptree->pointlist->current->surface_brightness = 0.0;
		MoveDownList(ptree->pointlist);
	}

	return;
}
/// Add an image to the map
void PixelMap::AddImages(
		ImageInfo *imageinfo   /// An array of ImageInfo-s.  There is no reason to separate images for this routine
		,int Nimages           /// Number of images on input.
		,bool constant_sb      /// true - all images will have surface brightness = 1,
		                       /// false - surface brightness is taken from surface_brighness in  the image points
		){

	if(Nimages <= 0) return;
	if(imageinfo->imagekist->Nunits() == 0) return;

	double sb = 1;
	for(long ii=0;ii<Nimages;++ii){
		if(imageinfo->imagekist->Nunits() > 0){
			MoveToTopKist(imageinfo[ii].imagekist);
			do{

				if(!constant_sb) sb = getCurrentKist(imageinfo[ii].imagekist)->surface_brightness;

				assert(getCurrentKist(imageinfo[ii].imagekist)->leaf);
				moveTop(ptree);
				_SplitFluxIntoPixels(ptree,getCurrentKist(imageinfo[ii].imagekist)->leaf,&sb);

			}while(MoveDownKist(imageinfo[ii].imagekist));
		}
	}

	return;
}

void PixelMap::Convert(){
	long ix;

	unsigned long i;
	if(map == NULL) map = new double[Npixels*Npixels];
	for(i=0 ; i < Npixels*Npixels ; ++i) map[i]=0.0;

	MoveToTopList(ptree->pointlist);
	for(i=0;i<ptree->pointlist->Npoints;++i){
		if(ptree->pointlist->current->surface_brightness > 0.0){
			ix = IndexFromPosition(ptree->pointlist->current->x,Npixels,range,center);
			if(ix > -1)	map[ix] =  ptree->pointlist->current->surface_brightness/resolution/resolution;
		}
		MoveDownList(ptree->pointlist);
	}

	return;
}

void PixelMap::print(){

	Convert();
	std::cout << Npixels << "  " << range << std::endl;
	for(unsigned long i=0;i < Npixels*Npixels; ++i) std::cout << map[i] << std::endl;
	std::cout << Npixels << "  " << range << std::endl;

	return;
}


/** \ingroup LowLevel
 *  Recursively determine which pixels a leaf intersects with
 *  and add flux to that pixel.  Used in pixelizer().
 */
void PixelMap::_SplitFluxIntoPixels(TreeHndl ptree,Branch *leaf,double *leaf_sb){

	double area = BoxIntersection(leaf,ptree->current);

	if(area > 0.0){
		if(atLeaf(ptree)){
		
			assert(ptree->current->npoints == 1);
			ptree->current->points->surface_brightness += (*leaf_sb)*area;
			//ptree->current->points->surface_brightness = *leaf_sb;

			return;
		}

		if(ptree->current->child1 != NULL){
			moveToChild(ptree,1);
			_SplitFluxIntoPixels(ptree,leaf,leaf_sb);
			moveUp(ptree);
		}

		if(ptree->current->child2 != NULL){
			moveToChild(ptree,2);
			_SplitFluxIntoPixels(ptree,leaf,leaf_sb);
			moveUp(ptree);
		}
	}

	return;
}
/** \ingroup Image
 *
 * \brief Smoothes a map with a Gaussian kernel. Needs to be tested.
 */
//TODO Ben, This should be improved and made a member of PixelMap
void smoothmap(double *map_out,double *map_in,long Npixels,double range,double sigma){
	double sum=0,**mask;
	unsigned long i=0;
	long j,k,ix,iy;
	int Nmask;

	Nmask=2*(int)(3*sigma*Npixels/range + 1);
	if(Nmask < 4 ) std::printf("WARNING: pixels are large compare to psf Nmask=%i\n",Nmask);

	// set up mask
	mask=dmatrix(0,Nmask,0,Nmask);
	for(j=-Nmask/2,sum=0;j<=Nmask/2;++j){
		for(k=-Nmask/2;k<=Nmask/2;++k){
			mask[j+Nmask/2][k+Nmask/2]= exp(-(pow(j*range/(Npixels-1),2)
					                        + pow(k*range/(Npixels-1),2))/2/pow(sigma,2) );
			sum+=mask[j+Nmask/2][k+Nmask/2];
		}
	}
	for(j=-Nmask/2;j<=Nmask/2;++j) for(k=-Nmask/2;k<=Nmask/2;++k) mask[j+Nmask/2][k+Nmask/2]/=sum;

	for(i=0;i<Npixels*Npixels;++i){
		for(j=0;j<=Nmask;++j){
			ix=i%Npixels + j-Nmask/2;
			if( (ix>-1)*(ix<Npixels) ){
				for(k=0;k<=Nmask;++k){
					iy=i/Npixels + k-Nmask/2;
					if( (iy>-1)*(iy<Npixels) ){

						map_out[ix+Npixels*iy] += mask[ix][iy]*map_in[i];
					}
				}
			}
		}
	}

	free_dmatrix(mask,0,Nmask,0,Nmask);

	return ;
}


void pixelize(
		double *map    /// Output map in one dimensional array. It is always square. map[0...Npixels*Npixels-1]
		,long Npixels   /// Number of pixels in one dimension of map.
		,double range   /// One dimensional range of map in whatever units the point positions are in (generally Mpc on the lens plane.)
		,double *center /// The location of the center of the map
		,ImageInfo *imageinfo  /// An array of ImageInfo-s.  There is no reason to separate images for this routine
		,int Nimages           /// Number of images on input.
		,bool constant_sb  /// true - all images will have surface brightness = 1,
		                      /// false - surface brightness is taken from surface_brighness in  the image points
		,bool cleanmap     ///  true - erases previous pixel map, false - adds new flux to map
		){

	if(imageinfo->imagekist->Nunits() == 0) return;

	long ix;
	double sb,resolution=0;
	unsigned long i,ii;
	Point *points;
	TreeHndl ptree;

	//printf("%d %g %g %g\n", Npixels, range, center[0], center[1]);

	if( (Npixels & (Npixels-1)) != 0){
		ERROR_MESSAGE();
		std::printf("ERROR: pixelsize, Npixels is not a power of 2\n");
		exit(1);
	}

	resolution=range/(Npixels-1);

	// initialize pixel tree
	points=NewPointArray(Npixels*Npixels,true);
	xygridpoints(points,range,center,Npixels,false);
	ptree=BuildTree(points,Npixels*Npixels);

	MoveToTopList(ptree->pointlist);
	for(i=0 ; i < ptree->pointlist->Npoints ; ++i){
		ptree->pointlist->current->surface_brightness = 0.0;
		MoveDownList(ptree->pointlist);
	}

	if(cleanmap)
		for(i=0 ; i < Npixels*Npixels ; ++i) map[i]=0.0;

	sb = 1;
	for(ii=0;ii<Nimages;++ii){
		MoveToTopKist(imageinfo[ii].imagekist);
		do{

			if(!constant_sb) sb = getCurrentKist(imageinfo[ii].imagekist)->surface_brightness;

			assert(getCurrentKist(imageinfo[ii].imagekist)->leaf);
			moveTop(ptree);
			//_SplitFluxIntoPixels(ptree,getCurrentKist(imageinfo[ii].imagekist)->leaf,&sb);

		}while(MoveDownKist(imageinfo[ii].imagekist));
	}

	int mycount;
	MoveToTopList(ptree->pointlist);
	for(i=0,mycount=0;i<ptree->pointlist->Npoints;++i){
		if(ptree->pointlist->current->surface_brightness > 0.0){
			ix = IndexFromPosition(ptree->pointlist->current->x,Npixels,range,center);
			if(ix > -1){
				mycount++;
				map[ix] =  ptree->pointlist->current->surface_brightness/resolution/resolution;
			}
		}
		MoveDownList(ptree->pointlist);
	}

	FreePointArray(points);

	std::cout << "Found " << mycount << " pixels!" << std::endl;

	return;
}
