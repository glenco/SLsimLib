/*
 * pixelize.c
 *
 *  Created on: Feb 27, 2010
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
#endif

PixelMap::PixelMap()
: isvalid(false), Npixels(0), range(0), center({0,0}), resolution(0),
  map_boundary_p1({0,0}), map_boundary_p2({0,0})
{
}

PixelMap::PixelMap(
		unsigned long my_Npixels  /// Number of pixels in one dimension of map.
		,double my_range          /// One dimensional range of map in whatever units the point positions are in
		,double *my_center        /// The location of the center of the map
		): isvalid(true), Npixels(my_Npixels), range(my_range)
		{

	center[0] = my_center[0];
	center[1] = my_center[1];

	resolution=range/Npixels;
	map_boundary_p1[0] = center[0]-range/2.;
	map_boundary_p1[1] = center[1]-range/2.;
	map_boundary_p2[0] = center[0]+range/2.;
	map_boundary_p2[1] = center[1]+range/2.;

	range = range - resolution;

	map.resize(Npixels*Npixels, 0);
	return;
}

PixelMap::PixelMap(std::string filename)
: isvalid(false)
		{
#ifdef ENABLE_FITS

		if(filename == ""){
			std::cout << "Please enter a valid filename for the FITS file input" << std::endl;
			exit(1);
		}

		std::auto_ptr<CCfits::FITS> fp (new CCfits::FITS (filename, CCfits::Read));
		CCfits::PHDU *h0=&fp->pHDU();
		const CCfits::ExtMap *h1=&fp->extension();
		Npixels = h0->axis(0);
		int ny=h0->axis(1);
		if(Npixels != ny){
			std::cout << "Only squared maps are allowed!" << std::endl;
			exit(1);
		}
		map.resize(Npixels*Npixels, 0);
		h0->readKey("CRVAL1",center[0]);
		h0->readKey("CRVAL2",center[1]);
		h0->readKey("CDELT1",resolution);
		std::cout << "Resolution is " << resolution << std::endl;
		resolution = fabs(resolution)*pi/180.;
		std::cout << "Resolution is " << resolution << std::endl;
		range = resolution*Npixels;
		map_boundary_p1[0] = center[0] - range/2.;
		map_boundary_p1[1] = center[1] - range/2.;
		map_boundary_p2[0] = center[0] + range/2.;
		map_boundary_p2[1] = center[1] + range/2.;
		range = range - resolution;
		h0->read(map);
		std::cout << "Resolution is " << resolution << std::endl;
		isvalid = true;

#else
		std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
		exit(1);
#endif

		}

PixelMap::~PixelMap(){

	map.resize(0);
}

bool PixelMap::valid() const
{
	return isvalid;
}

/// Zero the whole map
void PixelMap::Clean(){
	map = 0.0;

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
	std::list <unsigned long> neighborlist;
	std::list<unsigned long>::iterator it;
	for(long ii=0;ii<Nimages;++ii){

		if(imageinfo->imagekist->Nunits() > 0){
			MoveToTopKist(imageinfo[ii].imagekist);
			do{
				if(!constant_sb) sb = getCurrentKist(imageinfo[ii].imagekist)->surface_brightness;

				assert(getCurrentKist(imageinfo[ii].imagekist)->leaf);

				if ((inMapBox(getCurrentKist(imageinfo[ii].imagekist)->leaf)) == true){
					PointsWithinLeaf(getCurrentKist(imageinfo[ii].imagekist)->leaf,neighborlist);
					for(it = neighborlist.begin();it != neighborlist.end();it++){
						float area = LeafPixelArea(*it,getCurrentKist(imageinfo[ii].imagekist)->leaf);
						map[*it] += sb*area;
					}
				}
			}while(MoveDownKist(imageinfo[ii].imagekist));
		}
	}

	return;
}

void PixelMap::PointsWithinLeaf(Branch * branch1, std::list <unsigned long> &neighborlist){

	neighborlist.clear();

	int line_s,line_e,col_s,col_e;

	line_s = std::max(0,IndexFromPosition(branch1->boundary_p1[0],Npixels,range,center[0]));
	col_s = std::max(0,IndexFromPosition(branch1->boundary_p1[1],Npixels,range,center[1]));
	line_e = IndexFromPosition(branch1->boundary_p2[0],Npixels,range,center[0]);
	col_e = IndexFromPosition(branch1->boundary_p2[1],Npixels,range,center[1]);
	if (line_e < 0) line_e = Npixels-1;
	if (col_e < 0) col_e = Npixels-1;

	for (int iy = col_s; iy<= col_e; ++iy)
	{
		for (int ix = line_s; ix <= line_e; ++ix)
			{
				neighborlist.push_back(ix+Npixels*iy);
			}
		}
}

bool PixelMap::inMapBox(Branch * branch1){
	if (branch1->boundary_p1[0] > map_boundary_p2[0] || branch1->boundary_p2[0] < map_boundary_p1[0]) return false;
	if (branch1->boundary_p1[1] > map_boundary_p2[1] || branch1->boundary_p2[1] < map_boundary_p1[1]) return false;
	return true;
}

double PixelMap::LeafPixelArea(IndexType i,Branch * branch1){
	double area=0;
	PosType p[2],p1[2],p2[2];

	PositionFromIndex(i,p,Npixels,range,center);
	p1[0] = p[0] - .5*resolution;
	p1[1] = p[1] - .5*resolution;
	p2[0] = p[0] + .5*resolution;
	p2[1] = p[1] + .5*resolution;
	area = MIN(p2[0],branch1->boundary_p2[0])
	     - MAX(p1[0],branch1->boundary_p1[0]);
	if(area < 0) return 0.0;

	area *= MIN(p2[1],branch1->boundary_p2[1])
	      - MAX(p1[1],branch1->boundary_p1[1]);
	if(area < 0) return 0.0;

	return area;

}

/// Add an image to the map with Gaussian smoothing
void PixelMap::AddImages(
		ImageInfo *imageinfo   /// An array of ImageInfo-s.  There is no reason to separate images for this routine
		,int Nimages           /// Number of images on input.
		,double sigma          /// Gaussion width of smoothing kernal
		){

	if(sigma < resolution){
		ERROR_MESSAGE();
		std::cout << "ERROR in PixelMap::AddImages(), Smoothing scale must be larger than resolution of final image." << std::endl;
		exit(1);
	}
	if(Nimages <= 0) return;
	if(imageinfo->imagekist->Nunits() == 0) return;

	double sb,r[2],res,norm=0;
	KistHndl kist = new Kist();

	// find numerical normalization of mask on grid
//	PointsWithinKist(ptree,center,3*sigma,kist,0);
	kist->MoveToTop();
	do{
		r[0] = kist->getCurrent()->x[0] - center[0];
		r[1] = kist->getCurrent()->x[1] - center[1];
		norm += exp(-0.5*(r[0]*r[0] + r[1]*r[1] )/sigma/sigma);
	}while(kist->Down());


	for(long ii=0;ii<Nimages;++ii){
		if(imageinfo->imagekist->Nunits() > 0){
			MoveToTopKist(imageinfo[ii].imagekist);
			do{

				sb = getCurrentKist(imageinfo[ii].imagekist)->surface_brightness;
				res = getCurrentKist(imageinfo[ii].imagekist)->gridsize;

				if(res >= resolution){
					ERROR_MESSAGE();
					std::cout << "ERROR in PixelMap::AddImages(), Resolution of simulation must be higher than resolution of final image." << std::endl;
					exit(1);
				}
//				PointsWithinKist(ptree,imageinfo[ii].imagekist->getCurrent()->x,3*sigma,kist,0);
				kist->MoveToTop();
				do{
					r[0] = kist->getCurrent()->x[0] - imageinfo[ii].imagekist->getCurrent()->x[0];
					r[1] = kist->getCurrent()->x[1] - imageinfo[ii].imagekist->getCurrent()->x[1];
					kist->getCurrent()->surface_brightness += sb*res*res*exp(-0.5*(r[0]*r[0] + r[1]*r[1] )/sigma/sigma)/norm;
				}while(kist->Down());

			}while(MoveDownKist(imageinfo[ii].imagekist));
		}
	}

	delete kist;
	return;
}

/// Print an ASCII table of all the pixel values.
void PixelMap::printASCII(){

	std::cout << Npixels << "  " << range << std::endl;
	for(unsigned long i=0;i < Npixels*Npixels; ++i) std::cout << map[i] << std::endl;
	std::cout << Npixels << "  " << range << std::endl;

	//map.resize(0);
	return;
}
/// Print an ASCII table of all the pixel values.
void PixelMap::printASCIItoFile(std::string filename){
	std::ofstream file_map(filename.c_str());

	if(!file_map){
		std::cout << "unable to open file " << filename << std::endl;
		exit(0);
	}

	file_map << Npixels << "  " << range << std::endl;
	for(unsigned long i=0;i < Npixels*Npixels; ++i) file_map << std::scientific << map[i] << std::endl;
	file_map << Npixels << "  " << range << std::endl;

	//map.resize(0);

	file_map.close();

	return;
}
// Output the pixel map as a fits file.
void PixelMap::printFITS(std::string filename){

#ifdef ENABLE_FITS
		if(filename == ""){
			std::cout << "Please enter a valid filename for the FITS file output" << std::endl;
			exit(1);
		}

		//int Np = (int)Npixels;
		//writeImage(filename,map,Np,Np);

		long naxis=2;
		long naxes[2] = {Npixels,Npixels};

		std::auto_ptr<CCfits::FITS> fout(0);

		try{
			fout.reset(new CCfits::FITS(filename,FLOAT_IMG,naxis,naxes));
		}
		catch(CCfits::FITS::CantCreate){
			std::cout << "Unable to open fits file " << filename << std::endl;
			ERROR_MESSAGE();
			exit(1);
		}

		long seed;

		ran2(&seed)*1000;

		std::vector<long> naxex(2);
		naxex[0]=Npixels;
		naxex[1]=Npixels;

		CCfits::PHDU *phout = &fout->pHDU();

		phout->write( 1,Npixels*Npixels,map );
std::cout<< range << "  " << resolution << "  " << Npixels << std::endl;
		phout->addKey ("CRPIX1",naxex[0]/2,"");
		phout->addKey ("CRPIX2",naxex[1]/2,"");
		phout->addKey ("CRVAL1",0.0,"");
		phout->addKey ("CRVAL2",0.0,"");
		phout->addKey ("CDELT1",-180*range/(Npixels-1)/pi,"degrees");
		phout->addKey ("CDELT2", 180*range/(Npixels-1)/pi,"degrees");
		phout->addKey ("CTYPE1","RA--TAN","");
		phout->addKey ("CTYPE2","RA-TAN","");
		phout->addKey ("CROTA2",0.0,"");
		phout->addKey ("CD1_1",-180*range/(Npixels-1)/pi,"degrees");
		phout->addKey ("CD1_2",0.0,"");
		phout->addKey ("CD2_1",0.0,"");
		phout->addKey ("CD2_2", 180*range/(Npixels-1)/pi,"degrees");

		phout->addKey ("Npixels", Npixels,"");
		phout->addKey ("range ", range," radians");

		std::cout << *phout << std::endl;

		//map.resize(0);
#else
		std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
		exit(1);
#endif
}


/** \ingroup Image
 *
 * \brief Smoothes a map with a Gaussian kernel. Needs to be tested.
 *  This smoothes onto a map of the same resolution
 */
void PixelMap::smooth(double *map_out,double sigma){
	double sum=0,**mask;
	unsigned long i=0;
	long j,k,ix,iy;
	int Nmask;

	Nmask=2*(int)(3*sigma*Npixels/range + 1);
	if(Nmask < 4 ) std::printf("WARNING: pixels are large compare to psf Nmask=%i\n",Nmask);

	// set up mask
	Matrix(mask,Nmask,Nmask);
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

						map_out[ix+Npixels*iy] += mask[ix][iy]*map[i];
					}
				}
			}
		}
	}

	free_Matrix(mask,Nmask,Nmask);

	return ;
}

/*
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
		,bool write_for_skymaker /// true -- produces a fits map in the proper Skymaker format
		,std::string filename /// the filename for the FITS file
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

	if(write_for_skymaker == true){
#ifdef ENABLE_FITS
		if(filename == ""){
			std::cout << "Please enter a valid filename for the FITS file output" << std::endl;
			exit(1);
		}

		std::valarray<float> quantity;
		quantity.resize(Npixels*Npixels, 0);

		int j;
		for(i=0; i<Npixels; i++)
			for(j=0; j<Npixels; j++)
				quantity[j+i*Npixels] = map[j+i*Npixels];

		int Np = (int)Npixels;

		writeImage(filename,quantity,Np,Np);

		quantity.resize(0);
#else
		std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
		exit(1);
#endif
	}

	return;
}
*/
