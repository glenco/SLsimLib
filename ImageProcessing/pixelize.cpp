/*
 * pixelize.c
 *
 *  Created on: Feb 27, 2010
 *      Author: R.B. Metcalf
 */

#include "slsimlib.h"

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
//#include <CCfits>
#endif

#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

#include <fstream>
#include <algorithm>
#include <utility>

void swap(PixelMap& x, PixelMap& y)
{
	using std::swap;
	
	swap(x.map_size, y.map_size);
	swap(x.map, y.map);
	
	swap(x.Npixels, y.Npixels);
	swap(x.resolution, y.resolution);
	swap(x.range, y.range);
	
	swap(x.center[0], y.center[0]);
	swap(x.center[1], y.center[1]);
	
	swap(x.map_boundary_p1[0], y.map_boundary_p1[0]);
	swap(x.map_boundary_p1[1], y.map_boundary_p1[1]);
	swap(x.map_boundary_p2[0], y.map_boundary_p2[0]);
	swap(x.map_boundary_p2[1], y.map_boundary_p2[1]);
}
// TODO comment
bool agree(const PixelMap& a, const PixelMap& b)
{
	return (a.Npixels == b.Npixels) && (a.resolution == b.resolution) &&
		(a.center[0] == b.center[0] && a.center[1] == b.center[1]);
}

PixelMap::PixelMap()
: map_size(0), map(0), Npixels(0), resolution(0), range(0)
{
	center[0] = 0;
	center[1] = 0;
	
	map_boundary_p1[0] = 0;
	map_boundary_p1[1] = 0;
	map_boundary_p2[0] = 0;
	map_boundary_p2[1] = 0;
}

PixelMap::PixelMap(const PixelMap& other)
: map_size(other.map_size), map(0),
  Npixels(other.Npixels), resolution(other.resolution), range(other.range)
{
	if(map_size)
	{
		map = new float[map_size];
		std::copy(other.map, other.map + map_size, map);
	}
	
	std::copy(other.center, other.center + 2, center);
	
	std::copy(other.map_boundary_p1, other.map_boundary_p1 + 2, map_boundary_p1);
	std::copy(other.map_boundary_p2, other.map_boundary_p2 + 2, map_boundary_p2);
}

PixelMap::PixelMap(
		const double* center,  /// The location of the center of the map
		std::size_t Npixels,  /// Number of pixels in one dimension of map.
		double resolution        /// One dimensional range of map in whatever units the point positions are in
		): Npixels(Npixels), resolution(resolution)
		{

	std::copy(center, center + 2, this->center);
	range = resolution*(Npixels-1);
	
	map_boundary_p1[0] = center[0]-(Npixels*resolution)/2.;
	map_boundary_p1[1] = center[1]-(Npixels*resolution)/2.;
	map_boundary_p2[0] = center[0]+(Npixels*resolution)/2.;
	map_boundary_p2[1] = center[1]+(Npixels*resolution)/2.;

	map_size = Npixels*Npixels;
	map = new float[map_size];
	std::fill(map, map + map_size, 0);
}

/// Constructs a PixelMap reading in a fits file
/// Infos about resolution, Npixels and center are read from the header.
PixelMap::PixelMap(std::string filename)
{
#ifdef ENABLE_FITS

		if(filename == ""){
			std::cout << "Please enter a valid filename for the FITS file input" << std::endl;
			exit(1);
		}

		std::auto_ptr<CCfits::FITS> fp (new CCfits::FITS (filename, CCfits::Read));
		CCfits::PHDU *h0=&fp->pHDU();
		//const CCfits::ExtMap *h1=&fp->extension();
		Npixels = h0->axis(0);
		if(Npixels != (std::size_t)h0->axis(1))
		{
			std::cout << "Only squared maps are allowed!" << std::endl;
			exit(1);
		}
		h0->readKey("CRVAL1",center[0]);
		h0->readKey("CRVAL2",center[1]);
		h0->readKey("CDELT1",resolution);
		std::cout << "Resolution is " << resolution << std::endl;
		resolution = fabs(resolution)*pi/180.;
		std::cout << "Resolution is " << resolution << std::endl;
		range = resolution*(Npixels-1);
		map_boundary_p1[0] = center[0] - (Npixels*resolution)/2.;
		map_boundary_p1[1] = center[1] - (Npixels*resolution)/2.;
		map_boundary_p2[0] = center[0] + (Npixels*resolution)/2.;
		map_boundary_p2[1] = center[1] + (Npixels*resolution)/2.;
		
		std::valarray<float> image;
		h0->read(image);
		
		map_size = Npixels*Npixels;
		map = new float[map_size];
		std::copy(&image[0], &image[0] + map_size, map);
		
		std::cout << "Resolution is " << resolution << std::endl;
		
#else
		std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
		exit(1);
#endif

}

/// Creates a new PixelMap from a region of a PixelMap.
/// If the region exceeds the boundaries of the original map, the new map is completed with zeros.
PixelMap::PixelMap(const PixelMap& pmap,  /// Input PixelMap (from which the stamp is taken)
		const double* center, /// center of the region to be duplicated (in rads)
		std::size_t Npixels /// size of the region to be duplicated (in pixels)
		): resolution(pmap.resolution), Npixels(Npixels)
	{
		std::copy(center, center + 2, this->center);
		range = resolution*(Npixels-1);

		map_boundary_p1[0] = center[0]-(Npixels*resolution)/2.;
		map_boundary_p1[1] = center[1]-(Npixels*resolution)/2.;
		map_boundary_p2[0] = center[0]+(Npixels*resolution)/2.;
		map_boundary_p2[1] = center[1]+(Npixels*resolution)/2.;

		map_size = Npixels*Npixels;
		map = new float[map_size];

		int * edge = new int[2];
		edge[0] = (center[0]-pmap.map_boundary_p1[0])/resolution - Npixels/2;
		edge[1] = (center[1]-pmap.map_boundary_p1[1])/resolution - Npixels/2;
		if (edge[0] > pmap.Npixels || edge[1] > pmap.Npixels || edge[0]+Npixels < 0 || edge[1]+Npixels < 0)
		{
			std::cout << "The region you selected is completely outside PixelMap!" << std::endl;
			exit(1);
		}
		for (unsigned long i=0; i < map_size; ++i)
		{
			int ix = i%Npixels;
			int iy = i/Npixels;
			map[i] = 0;
			if (ix+edge[0] > 0 && ix+edge[0] < pmap.Npixels && iy+edge[1] > 0 && iy+edge[1] < pmap.Npixels)
				map[i] = pmap.map[ix+edge[0]+(iy+edge[1])*pmap.Npixels];
		}
	}

/// Creates a PixelMap at a different resolution.
/// The new counts are calculated integrating over the input pixels.
/// No interpolation or smoothing is performed.
PixelMap::PixelMap(const PixelMap& pmap, double res_ratio)
	{
	resolution = res_ratio*pmap.resolution;
	Npixels = pmap.Npixels/res_ratio + .5;
	range = resolution*(Npixels-1);
	center[0] = pmap.center[0];
	center[1] = pmap.center[1];
	map_boundary_p1[0] = center[0] - (Npixels*resolution)/2.;
	map_boundary_p1[1] = center[1] - (Npixels*resolution)/2.;
	map_boundary_p2[0] = center[0] + (Npixels*resolution)/2.;
	map_boundary_p2[1] = center[1] + (Npixels*resolution)/2.;
	
	map_size = Npixels*Npixels;
	map = new float[map_size];
	
	int old_Npixels = pmap.Npixels;
	int ix, iy;
	double area;
	double* old_p1 = new double[2];
	double* old_p2 = new double[2];

	for(unsigned long i=0;i < map_size; ++i)
	{
		ix = i%Npixels;
		iy = i/Npixels;
		map[ix+Npixels*iy] = 0.;
		old_p1[0] = std::max(0,int(ix*res_ratio));
		old_p1[1] = std::max(0,int(iy*res_ratio));
		old_p2[0] = std::min(old_Npixels-1,int((ix+1.)*res_ratio));
		old_p2[1] = std::min(old_Npixels-1,int((iy+1.)*res_ratio));
		for (int old_iy = old_p1[1]; old_iy<= old_p2[1]; ++old_iy)
		{
			for (int old_ix = old_p1[0]; old_ix<= old_p2[0]; ++old_ix)
				{
					area = MIN(old_ix+0.5,(ix+1.)*res_ratio-0.5) - MAX(old_ix-0.5,ix*res_ratio-0.5);
					area *= MIN(old_iy+0.5,(iy+1.)*res_ratio-0.5) - MAX(old_iy-0.5,iy*res_ratio-0.5);
					map[ix+Npixels*iy] += area*pmap.map[old_ix+old_Npixels*old_iy];
				}
			}
		}
}

PixelMap::~PixelMap()
{
	delete[] map;

}
/* TODO I don't like this because it is not what you would expect the "=" operator to do.
 * This is a swap and not a copy.
*/
PixelMap& PixelMap::operator=(PixelMap other)
{
	swap(*this, other);
	return *this;
}

/// Zero the whole map
void PixelMap::Clean()
{
	std::fill(map, map + map_size, 0);
}

/// Multiplies the whole map by a scalar factor
void PixelMap::Renormalize(double factor)
{
	for(std::size_t i=0;i < map_size; ++i) map[i] *= factor;
}

/// Adds a value to the i-th pixel
void PixelMap::AddValue(std::size_t i, double value)
{
	map[i] += value;
}

/// Assigns a value to the i-th pixel
void PixelMap::AssignValue(std::size_t i, double value)
{
	map[i] = value;
}

/// Add an image to the map
void PixelMap::AddImages(
		ImageInfo *imageinfo   /// An array of ImageInfo-s.  There is no reason to separate images for this routine
		,int Nimages           /// Number of images on input.
		,float rescale         /// rescales the surface brightness while leaving the image unchanged
			                   /// , rescale==0 gives constant surface brightness (default: 1)
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
				if(rescale != 0.0) sb = rescale*getCurrentKist(imageinfo[ii].imagekist)->surface_brightness;

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

/// returns the grid points within the branch
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
/// checks if the branch is within map boundaries
bool PixelMap::inMapBox(Branch * branch1){
	if (branch1->boundary_p1[0] > map_boundary_p2[0] || branch1->boundary_p2[0] < map_boundary_p1[0]) return false;
	if (branch1->boundary_p1[1] > map_boundary_p2[1] || branch1->boundary_p2[1] < map_boundary_p1[1]) return false;
	return true;
}

//// Finds the area of the intersection between pixel i and branch1
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

/*// Add an image to the map with Gaussian smoothing
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
	Kist<Point> * kist = new Kist<Point>();

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
*/
/// Print an ASCII table of all the pixel values.
void PixelMap::printASCII(){

	std::cout << Npixels << "  " << range << std::endl;
	for(std::size_t i=0;i < map_size; ++i) std::cout << map[i] << std::endl;
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
	for(std::size_t i=0;i < map_size; ++i) file_map << std::scientific << map[i] << std::endl;
	file_map << Npixels << "  " << range << std::endl;

	//map.resize(0);

	file_map.close();

	return;
}
/// Output the pixel map as a fits file.
void PixelMap::printFITS(std::string filename){

#ifdef ENABLE_FITS
		if(filename == ""){
			std::cout << "Please enter a valid filename for the FITS file output" << std::endl;
			exit(1);
		}

		//int Np = (int)Npixels;
		//writeImage(filename,map,Np,Np);

		long naxis=2;
		long naxes[2] = {(long)Npixels,(long)Npixels};

		std::auto_ptr<CCfits::FITS> fout(0);

		try{
			fout.reset(new CCfits::FITS(filename,FLOAT_IMG,naxis,naxes));
		}
		catch(CCfits::FITS::CantCreate&){
			std::cout << "Unable to open fits file " << filename << std::endl;
			ERROR_MESSAGE();
			exit(1);
		}

		//long seed;

		//ran2(&seed)*1000;

		std::vector<long> naxex(2);
		naxex[0]=Npixels;
		naxex[1]=Npixels;

		CCfits::PHDU *phout = &fout->pHDU();
		phout->write(1, map_size, std::valarray<float>(map, map_size));
std::cout<< range << "  " << resolution << "  " << Npixels << std::endl;		
		phout->addKey ("CRPIX1",naxex[0]/2,"");
		phout->addKey ("CRPIX2",naxex[1]/2,"");
		phout->addKey ("CRVAL1",0.0,"");
		phout->addKey ("CRVAL2",0.0,"");
		phout->addKey ("CDELT1",-180*resolution/pi,"degrees");
		phout->addKey ("CDELT2", 180*resolution/pi,"degrees");
		phout->addKey ("CTYPE1","RA--TAN","");
		phout->addKey ("CTYPE2","RA-TAN","");
		phout->addKey ("CROTA2",0.0,"");
		phout->addKey ("CD1_1",-180*resolution/pi,"degrees");
		phout->addKey ("CD1_2",0.0,"");
		phout->addKey ("CD2_1",0.0,"");
		phout->addKey ("CD2_2", 180*resolution/pi,"degrees");

		phout->addKey ("Npixels", Npixels,"");
		phout->addKey ("range ", map_boundary_p2[0]-map_boundary_p1[0]," radians");

		std::cout << *phout << std::endl;

		//map.resize(0);
#else
		std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
		exit(1);
#endif
}


/** \ingroup Image
 *
 * \brief Smoothes a map with a Gaussian kernel of width sigma (in arcseconds)
 */
void PixelMap::smooth(double sigma){
	double sum=0,**mask;
	int ix,iy;
	int Nmask,Nmask_half;
	double *map_out;
	int j_cen, k_cen;

	sigma /= 3600.*180/pi;
	map_out = new double[map_size];
	Nmask=2*(int)(3*sigma/resolution + 1);
	std::cout << Nmask << std::endl;
	if(Nmask < 4 ) std::cout << "WARNING: pixels are large compare to psf Nmask=" << Nmask << std::endl;

	Nmask_half = int(Nmask/2);
	mask = new double*[Nmask];
	for (int j = 0; j <Nmask; j++)
		mask[j] = new double[Nmask];

	for(int j=0;j<Nmask;j++)
	{
		for(int k=0;k<Nmask;k++)
		{
			j_cen = j - Nmask_half;
			k_cen = k - Nmask_half;
			mask[j][k]= exp(-(pow(j_cen*resolution,2) + pow(k_cen*resolution,2))/2/pow(sigma,2) );
			sum+=mask[j][k];
		}
	}
	for(int j=0;j<Nmask;j++)
	{
		for(int k=0;k<Nmask;k++)
		{
			mask[j][k]/=sum;
			std::cout << mask[j][k] << std::endl;
		}
	}
	for(long i=0;i<map_size;i++){
		map_out[i] = 0.;
	}

	for(long i=0;i<map_size;i++){
		for(int j=0;j<Nmask;j++){
			ix=i%Npixels + j-Nmask_half;
			if( (ix>-1) && (ix<Npixels) ){
				for(int k=0;k<Nmask;k++){
					iy=i/Npixels + k-Nmask_half;
					if( (iy>-1) && (iy<Npixels) ){
						map_out[ix+Npixels*iy] += mask[j][k]*map[i];
					}
				}
			}
		}
	}
	for(long i=0;i<map_size;i++){
		map[i] = map_out[i];
	}

	for (int j = 0; j <Nmask; j++)
		delete mask[j];
	delete [] mask;
	delete [] map_out;

}

/** \brief Smooths the image with a PSF map.
*
*/
PixelMap Observation::ApplyPSF(PixelMap &pmap)
{
	if (map_psf.size() == 0)
	{
		if (seeing > 0.)
		{
			PixelMap outmap(pmap);
			outmap.smooth(seeing);
			return outmap;
		}
		else
		{
			return pmap;
		}
	}
	else
	{
#ifdef ENABLE_FITS
#ifdef ENABLE_FFTW

	PixelMap outmap(pmap);
	// creates plane for fft of map, sets properly input and output data, then performs fft
	fftw_plan p;
	long Npix = outmap.getNpixels();
	double* in = new double[Npix*Npix];
	std::complex<double>* out=new std::complex<double> [Npix*(Npix/2+1)];
	for (unsigned long i = 0; i < Npix*Npix; i++)
	{
		in[i] = outmap[i];
	}
	p = fftw_plan_dft_r2c_2d(Npix,Npix,in, reinterpret_cast<fftw_complex*>(out), FFTW_ESTIMATE);
	fftw_execute(p);

	// creates plane for perform backward fft after convolution, sets output data
	fftw_plan p2;
	double* out2 = new double[Npix*Npix];
	p2 = fftw_plan_dft_c2r_2d(Npix,Npix,reinterpret_cast<fftw_complex*>(out), out2, FFTW_ESTIMATE);

	// calculates normalisation of psf
	int N_psf = map_psf.size();
	int side_psf = sqrt(N_psf);
	double map_norm = 0.;
	for (int i = 0; i < N_psf; i++)
	{
		map_norm += map_psf[i];
	}
	fftw_plan p_psf;

	// arrange psf data for fft, creates plane, then performs fft
	int psf_big_Npixels = static_cast<int>(Npix*oversample);
	double* psf_big = new double[psf_big_Npixels*psf_big_Npixels];
	std::complex<double>* out_psf=new std::complex<double> [psf_big_Npixels*(psf_big_Npixels/2+1)];
	p_psf = fftw_plan_dft_r2c_2d(psf_big_Npixels,psf_big_Npixels,psf_big, reinterpret_cast<fftw_complex*>(out_psf), FFTW_ESTIMATE);
	long ix, iy;
	for (int i = 0; i < psf_big_Npixels*psf_big_Npixels; i++)
	{
		ix = i/psf_big_Npixels;
		iy = i%psf_big_Npixels;
		if(ix<side_psf/2 && iy<side_psf/2)
			psf_big[i] = map_psf[(ix+side_psf/2)*side_psf+(iy+side_psf/2)]/map_norm;
		else if(ix<side_psf/2 && iy>=psf_big_Npixels-side_psf/2)
			psf_big[i] = map_psf[(ix+side_psf/2)*side_psf+(iy-(psf_big_Npixels-side_psf/2))]/map_norm;
		else if(ix>=psf_big_Npixels-side_psf/2 && iy<side_psf/2)
			psf_big[i] = map_psf[(ix-(psf_big_Npixels-side_psf/2))*side_psf+(iy+side_psf/2)]/map_norm;
		else if(ix>=psf_big_Npixels-side_psf/2 && iy>=psf_big_Npixels-side_psf/2)
			psf_big[i] = map_psf[(ix-(psf_big_Npixels-side_psf/2))*side_psf+(iy-(psf_big_Npixels-side_psf/2))]/map_norm;
		else
			psf_big[i] = 0.;
	}
	fftw_execute(p_psf);

	// performs convolution in Fourier space, and transforms back to real space
	for (unsigned long i = 0; i < Npix*(Npix/2+1); i++)
	{
		ix = i/(Npix/2+1);
		iy = i%(Npix/2+1);
		if (ix>Npix/2)
			out[i] *= out_psf[(psf_big_Npixels-(Npix-ix))*(psf_big_Npixels/2+1)+iy];
		else
			out[i] *= out_psf[ix*(psf_big_Npixels/2+1)+iy];
	}
	fftw_execute(p2);

	// translates array of data in (normalised) counts map
	for (unsigned long i = 0; i < Npix*Npix; i++)
	{
		ix = i/Npix;
		iy = i%Npix;
		outmap.AssignValue(i,out2[i]/double(Npix*Npix));
	}
	return outmap;
#else
		std::cout << "Please enable the preprocessor flag ENABLE_FFTW !" << std::endl;
		exit(1);
#endif

#else
		std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
		exit(1);
#endif
	}
}

/// Applies realistic noise (read-out + Poisson) on an image
PixelMap Observation::AddNoise(PixelMap &pmap)
{
	PixelMap outmap(pmap);
	double Q = pow(10,0.4*(mag_zeropoint+48.6));
	double res_in_arcsec = outmap.getResolution()*180.*60.*60/pi;
	double back_mean = pow(10,-0.4*(48.6+back_mag))*res_in_arcsec*res_in_arcsec*Q*exp_time;
	double rms, noise;
	long seed = 24;
	double norm_map, noised_map;
	for (unsigned long i = 0; i < outmap.getNpixels()*outmap.getNpixels(); i++)
	{
		norm_map = outmap[i]*exp_time;
		rms = sqrt(pow(exp_num*ron,2)+norm_map+back_mean);
		noise = gasdev(&seed)*rms;
		outmap.AddValue(i,noise/exp_time);
		}
	return outmap;
	}

/// Translates photon flux (in photons/(cm^2*Hz)) into telescope pixel counts
PixelMap Observation::PhotonToCounts(PixelMap &pmap)
{
	PixelMap outmap(pmap);
	double Q = pow(10,0.4*(mag_zeropoint+48.6))*hplanck;
	outmap.Renormalize(Q);
	return outmap;
}

// TODO comment (First reorganize the constructors, probably we do not need so many of them)
Observation::Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, float seeing):
		exp_time(exp_time), exp_num(exp_num), back_mag(back_mag), diameter(diameter), transmission(transmission), ron(ron), seeing(seeing)
		{
			mag_zeropoint = 2.5*log10(diameter*diameter*transmission*pi/4./hplanck) - 48.6;
		}
// TODO comment
Observation::Observation(float diameter, float transmission, float exp_time, int exp_num, float back_mag, float ron, std::string psf_file, float oversample):
		exp_time(exp_time), exp_num(exp_num), back_mag(back_mag), diameter(diameter), transmission(transmission), ron(ron), oversample(oversample)
		{
	mag_zeropoint = 2.5*log10(diameter*diameter*transmission*pi/4./hplanck) - 48.6;

	std::auto_ptr<CCfits::FITS> fp (new CCfits::FITS (psf_file.c_str(), CCfits::Read));
	CCfits::PHDU *h0=&fp->pHDU();
	int side_psf = h0->axis(0);
	int N_psf = side_psf*side_psf;
	map_psf.resize(N_psf);
	h0->read(map_psf);
		}
// TODO comment
Observation::Observation(float diameter, float transmission):
		diameter(diameter), transmission(transmission)
		{
	mag_zeropoint = 2.5*log10(diameter*diameter*transmission*pi/4./hplanck) - 48.6;
		}
// TODO comment
Observation::Observation(float zeropoint):
		mag_zeropoint(zeropoint)
{
}

void swap(PixelMask& x, PixelMask& y)
	{

	using std::swap;
	
	swap(x.map_size, y.map_size);
	swap(x.mask_size, y.mask_size);
	swap(x.pixels, y.pixels);
}

PixelMask::PixelMask()
: map_size(0), mask_size(0)
{
}

PixelMask::PixelMask(std::size_t map_size)
: map_size(map_size), mask_size(map_size)
{
	pixels.resize(map_size);
	for(std::size_t i = 0; i < map_size; ++i)
		pixels[i] = true;
}

PixelMask::PixelMask(const PixelMap& base, double threshold, ThresholdType type)
: map_size(base.size())
{
	pixels.reserve(map_size);
	
	switch(type)
	{
	case Greater:
		for(std::size_t i = 0; i < map_size; ++i)
			if(base[i] > threshold)
				pixels.push_back(i);
		break;
	case GreaterOrEqual:
		for(std::size_t i = 0; i < map_size; ++i)
			if(base[i] >= threshold)
				pixels.push_back(i);
		break;
	case Less:
		for(std::size_t i = 0; i < map_size; ++i)
			if(base[i] < threshold)
				pixels.push_back(i);
		break;
	case LessOrEqual:
		for(std::size_t i = 0; i < map_size; ++i)
			if(base[i] <= threshold)
				pixels.push_back(i);
		break;
	}
	
	mask_size = pixels.size();
	pixels.resize(mask_size);
}

PixelMask::PixelMask(std::string file, double threshold, ThresholdType type)
{
#ifdef ENABLE_FITS
	// get PixelMap for file
	PixelMap map(file);
	
	// create PixelMask for map
	PixelMask mask(map);
	
	// swap contents of mask into this
	swap(*this, mask);
#else
	// warn about ENABLE_FITS
	SLSIMLIB_DEBUG("FITS support disabled, use ENABLE_FITS flag");
	
	// default mask
	PixelMask mask;
	
	// swap defaults into this
	swap(*this, mask);
#endif
}

PixelMask& PixelMask::operator=(PixelMask other)
{
	swap(*this, other);
	return *this;
}

std::size_t PixelMask::operator[](std::size_t i) const
{
	return pixels[i];
}

bool PixelMask::valid() const
{
	return map_size;
}

bool PixelMask::empty() const
{
	return !mask_size;
}

std::size_t PixelMask::size() const
{
	return mask_size;
}

std::size_t PixelMask::base_size() const
{
	return map_size;
}

void swap(PixelData& a, PixelData& b)
{
	using std::swap;
	
	swap(a.img, b.img);
	swap(a.noi, b.noi);
}

PixelData::PixelData(const PixelMap& image, const PixelMap& noise, double zp_mag)
: img(image), noi(noise), zp(zp_mag)
{
	// must be a valid image
	assert(img.valid());
	
	// must be a valid sigma map
	assert(noi.valid());
	
	// image and sigma need to agree
	assert(agree(img, noi));
	
	// renormalize pixmaps by amount given from zero-point magnitude.
	double norm = std::pow(10., -0.4*(zp_mag+48.6))*inv_hplanck;
	img.Renormalize(norm);
	noi.Renormalize(norm);
}

PixelData::PixelData(const PixelData& other)
: img(other.img), noi(other.noi), zp(other.zp)
{
}

PixelData& PixelData::operator=(PixelData rhs)
{
	swap(rhs, *this);
	return *this;
}

double PixelData::chi_square(const PixelMap &model) const
{
	assert(model.valid());
	assert(agree(model, img));
	
	double chi2 = 0;
	for(std::size_t i = 0, n = model.size(); i < n; ++i)
		chi2 += std::pow(img[i] - model[i], 2)/(model[i] + noi[i]);
	return chi2;
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
