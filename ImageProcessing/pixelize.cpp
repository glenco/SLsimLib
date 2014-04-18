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

#include <fstream>
#include <algorithm>
#include <utility>
#include <stdexcept>

#if __cplusplus < 201103L
template<typename T>
void swap(std::valarray<T>& x, std::valarray<T>& y)
{
	std::valarray<T> z(x);
	
	x.resize(y.size());
	x = y;
	
	y.resize(z.size());
	y = z;
}
#endif

void swap(PixelMap& x, PixelMap& y)
{
	using std::swap;
	
	swap(x.map,y.map);
	
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

PixelMap::PixelMap()
: map(), Npixels(0), resolution(0), range(0)
{
	center[0] = 0;
	center[1] = 0;
	
	map_boundary_p1[0] = 0;
	map_boundary_p1[1] = 0;
	map_boundary_p2[0] = 0;
	map_boundary_p2[1] = 0;
}

PixelMap::PixelMap(const PixelMap& other)
: map(other.map),
  Npixels(other.Npixels), resolution(other.resolution), range(other.range)
{
	std::copy(other.center, other.center + 2, center);
	
	std::copy(other.map_boundary_p1, other.map_boundary_p1 + 2, map_boundary_p1);
	std::copy(other.map_boundary_p2, other.map_boundary_p2 + 2, map_boundary_p2);
}

PixelMap::PixelMap(
		const PosType* center,  /// The location of the center of the map
		std::size_t Npixels,  /// Number of pixels in one dimension of map.
		PosType resolution        /// One dimensional range of map in whatever units the point positions are in
		)
: map(0.0, Npixels*Npixels),
  Npixels(Npixels), resolution(resolution)
{
	std::copy(center, center + 2, this->center);
	range = resolution*(Npixels-1);
	
	map_boundary_p1[0] = center[0]-(Npixels*resolution)/2.;
	map_boundary_p1[1] = center[1]-(Npixels*resolution)/2.;
	map_boundary_p2[0] = center[0]+(Npixels*resolution)/2.;
	map_boundary_p2[1] = center[1]+(Npixels*resolution)/2.;
}

/** \brief Constructs a PixelMap reading in a fits file
* Infos about resolution, Npixels and center are read from the header.
 */
PixelMap::PixelMap(std::string fitsfilename)
{
#ifdef ENABLE_FITS
	if(fitsfilename.empty())
		throw std::invalid_argument("Please enter a valid filename for the FITS file input");
	
	std::auto_ptr<CCfits::FITS> fp(new CCfits::FITS(fitsfilename, CCfits::Read));
	
	CCfits::PHDU& h0 = fp->pHDU();
	
	//const CCfits::ExtMap *h1=&fp->extension();
	
	Npixels = h0.axis(0);
	if(Npixels != (std::size_t)h0.axis(1))
		throw std::runtime_error("Only square maps are allowed!");
	
	h0.readKey("CRVAL1", center[0]);
	h0.readKey("CRVAL2", center[1]);
	
	// read the resolution
	try
	{
		double cdelt2;
		h0.readKey("CDELT1", resolution);
		h0.readKey("CDELT2", cdelt2);
		if(std::abs(resolution) - std::abs(cdelt2) > 1e-6)
			throw std::runtime_error("non-square pixels in FITS file " + fitsfilename);
	}
	catch(CCfits::HDU::NoSuchKeyword&)
	{
		double cd12, cd21, cd22;
		h0.readKey("CD1_1", resolution);
		h0.readKey("CD1_2", cd12);
		h0.readKey("CD2_1", cd21);
		h0.readKey("CD2_2", cd22);
		if(std::abs(resolution) - std::abs(cd22) > 1e-6)
			throw std::runtime_error("non-square pixels in FITS file " + fitsfilename);
		if(cd12 || cd21)
			throw std::runtime_error("pixels not aligned with coordingate in FITS file " + fitsfilename);
	}
	
	resolution = fabs(resolution)*pi/180.;
	range = resolution*(Npixels-1);
	map_boundary_p1[0] = center[0] - (Npixels*resolution)/2.;
	map_boundary_p1[1] = center[1] - (Npixels*resolution)/2.;
	map_boundary_p2[0] = center[0] + (Npixels*resolution)/2.;
	map_boundary_p2[1] = center[1] + (Npixels*resolution)/2.;
	
	h0.read(map);
#else
	std::cerr << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif
}

/** \brief Creates a new PixelMap from a region of a PixelMap.
 * If the region exceeds the boundaries of the original map, the new map is completed with zeros.
 */
PixelMap::PixelMap(const PixelMap& pmap,  /// Input PixelMap (from which the stamp is taken)
		const PosType* center, /// center of the region to be duplicated (in rads)
		std::size_t Npixels /// size of the region to be duplicated (in pixels)
		)
: map(0.0, Npixels*Npixels),
  Npixels(Npixels), resolution(pmap.resolution)
	{
		std::copy(center, center + 2, this->center);
		range = resolution*(Npixels-1);

		map_boundary_p1[0] = center[0]-(Npixels*resolution)/2.;
		map_boundary_p1[1] = center[1]-(Npixels*resolution)/2.;
		map_boundary_p2[0] = center[0]+(Npixels*resolution)/2.;
		map_boundary_p2[1] = center[1]+(Npixels*resolution)/2.;

		int * edge = new int[2];
		edge[0] = (center[0]-pmap.map_boundary_p1[0])/resolution - Npixels/2;
		edge[1] = (center[1]-pmap.map_boundary_p1[1])/resolution - Npixels/2;
		if (edge[0] > int(pmap.Npixels) || edge[1] > int(pmap.Npixels) || edge[0]+int(Npixels) < 0 || edge[1]+int(Npixels) < 0)
		{
			std::cout << "The region you selected is completely outside PixelMap!" << std::endl;
			throw std::runtime_error("Attempting to make Sub-PixelMap outside of parent PixelMap!");
		}
		for (unsigned long i=0; i < map.size(); ++i)
		{
			int ix = i%Npixels;
			int iy = i/Npixels;
			map[i] = 0;
			if (ix+edge[0] > 0 && ix+edge[0] < pmap.Npixels && iy+edge[1] > 0 && iy+edge[1] < pmap.Npixels)
				map[i] = pmap.map[ix+edge[0]+(iy+edge[1])*pmap.Npixels];
		}
	}

/** \brief Creates a PixelMap at a different resolution.
* The new counts are calculated integrating over the input pixels.
* No interpolation or smoothing is performed.
 */
PixelMap::PixelMap(
                   const PixelMap& pmap
                   , PosType res_ratio     /// resolution of map is res_ratio times the resolution of the input map
                   )
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
	
	map.resize(Npixels*Npixels);
	
	int old_Npixels = pmap.Npixels;
	int ix, iy;
	PosType area;
	PosType* old_p1 = new PosType[2];
	PosType* old_p2 = new PosType[2];

	for(unsigned long i=0;i < map.size(); ++i)
	{
		ix = i%Npixels;
		iy = i/Npixels;
		map[ix+Npixels*iy] = 0.;
		old_p1[0] = std::max(0,int(ix*res_ratio));
		old_p1[1] = std::max(0,int(iy*res_ratio));
		old_p2[0] = std::min(old_Npixels-1,int((ix+1.)*res_ratio));
		old_p2[1] = std::min(old_Npixels-1,int((iy+1.)*res_ratio));
		for (int old_iy = old_p1[1]; old_iy <= old_p2[1]; ++old_iy)
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
}

PixelMap& PixelMap::operator=(PixelMap other)
{
	swap(*this, other);
	return *this;
}

/// Zero the whole map
void PixelMap::Clean()
{
	map = 0;
}

/// Multiplies the whole map by a scalar factor
void PixelMap::Renormalize(PosType factor)
{
	map *= factor;
}

/// Adds a value to the i-th pixel
void PixelMap::AddValue(std::size_t i, PosType value)
{
	map[i] += value;
}

/// Assigns a value to the i-th pixel
void PixelMap::AssignValue(std::size_t i, PosType value)
{
	map[i] = value;
}

/// Check whether two PixelMaps agree in their physical dimensions.
bool PixelMap::agrees(const PixelMap& other) const
{
	return
		(Npixels == other.Npixels) &&
		(resolution == other.resolution) &&
		(center[0] == other.center[0]) &&
		(center[1] == other.center[1]);
}

/// Add the values of another PixelMap to this one.
PixelMap& PixelMap::operator+=(const PixelMap& rhs)
{
	// TODO: maybe check if PixelMaps agree, but is slower
    if(Npixels != rhs.getNpixels()) throw std::runtime_error("Dimensions of maps are not compatible");
	for(size_t i=0;i<map.size();++i) map[i] += rhs.map[i];
	return *this;
}

/// Add two PixelMaps.
PixelMap operator+(const PixelMap& a, const PixelMap& b)
{
	PixelMap sum(a);
	sum += b;
	return sum;
}

/// Subtract the values of another PixelMap from this one.
PixelMap& PixelMap::operator-=(const PixelMap& rhs)
{
	// TODO: maybe check if PixelMaps agree, but is slower
	map -= rhs.map;
	return *this;
}

/// Subtract two PixelMaps.
PixelMap operator-(const PixelMap& a, const PixelMap& b)
{
	PixelMap diff(a);
	diff -= b;
	return diff;
}


/** \brief Add an image to the map
 *
 *  If rescale==0 gives constant surface brightness, if < 0
 *  the surface brightness is not scaled by the pixel area as for the flux (default: 1).
 *  Negative values are good for mapping some quantity independant of the pixel size
 */
void PixelMap::AddImages(
		ImageInfo *imageinfo   /// An array of ImageInfo-s.  There is no reason to separate images for this routine
		,int Nimages           /// Number of images on input.
		,float rescale         /// rescales the surface brightness while leaving the image unchanged,
                               ///  see full notes
		){

	if(Nimages <= 0) return;
	if(imageinfo->imagekist->Nunits() == 0) return;

	PosType sb = 1;
    float area = 1;
	std::list <unsigned long> neighborlist;
	std::list<unsigned long>::iterator it;
	for(long ii=0;ii<Nimages;++ii){

		if(imageinfo->imagekist->Nunits() > 0){
			imageinfo[ii].imagekist->MoveToTop();
			do{
				if(rescale != 0.0) sb = fabs(rescale)*imageinfo[ii].imagekist->getCurrent()->surface_brightness;

				assert(imageinfo[ii].imagekist->getCurrent()->leaf);

				if ((inMapBox(imageinfo[ii].imagekist->getCurrent()->leaf)) == true){
					PointsWithinLeaf(imageinfo[ii].imagekist->getCurrent()->leaf,neighborlist);
					for(it = neighborlist.begin();it != neighborlist.end();it++){
						area = LeafPixelArea(*it,imageinfo[ii].imagekist->getCurrent()->leaf);
            map[*it] += sb*area;
					}
				}
			}while(imageinfo[ii].imagekist->Down());
		}
	}
    
    if(rescale < 0){
        for(size_t i=0; i< Npixels*Npixels ;++i) map[i] /= resolution*resolution;
    }

	return;
}

/// returns the grid points within the branch
void PixelMap::PointsWithinLeaf(Branch * branch1, std::list <unsigned long> &neighborlist){

	neighborlist.clear();

	int line_s,line_e,col_s,col_e;

	line_s = std::max(0,int(Utilities::IndexFromPosition(branch1->boundary_p1[0],Npixels,range,center[0])));
	col_s = std::max(0,int(Utilities::IndexFromPosition(branch1->boundary_p1[1],Npixels,range,center[1])));
	line_e = Utilities::IndexFromPosition(branch1->boundary_p2[0],Npixels,range,center[0]);
	col_e = Utilities::IndexFromPosition(branch1->boundary_p2[1],Npixels,range,center[1]);
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
/// checks if point is within map boundaries
bool PixelMap::inMapBox(PosType * x){
	if (x[0] > map_boundary_p2[0] || x[0] < map_boundary_p1[0]) return false;
	if (x[1] > map_boundary_p2[1] || x[1] < map_boundary_p1[1]) return false;
	return true;
}

//// Finds the area of the intersection between pixel i and branch1
PosType PixelMap::LeafPixelArea(IndexType i,Branch * branch1){
	PosType area=0;
	PosType p[2],p1[2],p2[2];

	Utilities::PositionFromIndex(i,p,Npixels,range,center);
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
		,PosType sigma          /// Gaussion width of smoothing kernal
		){

	if(sigma < resolution){
		ERROR_MESSAGE();
		std::cout << "ERROR in PixelMap::AddImages(), Smoothing scale must be larger than resolution of final image." << std::endl;
		exit(1);
	}
	if(Nimages <= 0) return;
	if(imageinfo->imagekist->Nunits() == 0) return;

	PosType sb,r[2],res,norm=0;
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
void PixelMap::printASCII() const
{
	std::cout << Npixels << "  " << range << std::endl;
	for(std::size_t i=0;i < map.size(); ++i) std::cout << map[i] << std::endl;
	std::cout << Npixels << "  " << range << std::endl;

	//map.resize(0);
	return;
}
/// Print an ASCII table of all the pixel values.
void PixelMap::printASCIItoFile(std::string filename) const
{
	std::ofstream file_map(filename.c_str());

	if(!file_map){
		std::cout << "unable to open file " << filename << std::endl;
		exit(0);
	}

	file_map << Npixels << "  " << range << std::endl;
	for(std::size_t i=0;i < map.size(); ++i) file_map << std::scientific << map[i] << std::endl;
	file_map << Npixels << "  " << range << std::endl;

	//map.resize(0);

	file_map.close();

	return;
}
/// Output the pixel map as a fits file.
void PixelMap::printFITS(std::string filename, bool verbose) const
{
#ifdef ENABLE_FITS
	if(filename.empty())
		throw std::invalid_argument("Please enter a valid filename for the FITS file output");
	
	long naxis = 2;
	long naxes[2] = {(long)Npixels, (long)Npixels};
	
	// might throw CCfits::FITS::CantCreate
	std::auto_ptr<CCfits::FITS> fout(new CCfits::FITS(filename, FLOAT_IMG, naxis, naxes));
	
	std::vector<long> naxex(2);
	naxex[0] = Npixels;
	naxex[1] = Npixels;
	
	CCfits::PHDU& phout = fout->pHDU();
	
	phout.write(1, map.size(), map);
	
	phout.addKey("WCSAXES", 2, "number of World Coordinate System axes");
	phout.addKey("CRPIX1", 0.5*(naxex[0]+1), "x-coordinate of reference pixel");
	phout.addKey("CRPIX2", 0.5*(naxex[1]+1), "y-coordinate of reference pixel");
	phout.addKey("CRVAL1", 0.0, "first axis value at reference pixel");
	phout.addKey("CRVAL2", 0.0, "second axis value at reference pixel");
	phout.addKey("CTYPE1", "RA---TAN", "the coordinate type for the first axis");
	phout.addKey("CTYPE2", "DEC--TAN", "the coordinate type for the second axis");
	phout.addKey("CUNIT1", "deg     ", "the coordinate unit for the first axis");
	phout.addKey("CUNIT2", "deg     ", "the coordinate unit for the second axis");
	phout.addKey("CDELT1", -180*resolution/pi, "partial of first axis coordinate w.r.t. x");
	phout.addKey("CDELT2", 180*resolution/pi, "partial of second axis coordinate w.r.t. y");
	phout.addKey("CROTA2", 0.0, "");
	phout.addKey("CD1_1", -180*resolution/pi, "partial of first axis coordinate w.r.t. x");
	phout.addKey("CD1_2", 0.0, "partial of first axis coordinate w.r.t. y");
	phout.addKey("CD2_1", 0.0, "partial of second axis coordinate w.r.t. x");
	phout.addKey("CD2_2", 180*resolution/pi, "partial of second axis coordinate w.r.t. y");
	
	phout.addKey("Npixels", Npixels, "");
	phout.addKey("range", map_boundary_p2[0]-map_boundary_p1[0], "radians");
	phout.addKey("RA", center[0], "radians");
	phout.addKey("DEC", center[1], "radians");
    
	if(verbose)
		std::cout << phout << std::endl;
#else
	std::cerr << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif
}

/** \ingroup Image
 *
 * \brief Smoothes a map with a Gaussian kernel of width sigma (in arcseconds)
 */
void PixelMap::smooth(PosType sigma){
	PosType sum=0,**mask;
	int ix,iy;
	int Nmask,Nmask_half;
	int j_cen, k_cen;

	sigma /= 3600.*180/pi;
	Nmask=2*(int)(3*sigma/resolution + 1);
	std::cout << Nmask << std::endl;
	if(Nmask < 4 ) std::cout << "WARNING: pixels are large compare to psf Nmask=" << Nmask << std::endl;

	Nmask_half = int(Nmask/2);
	mask = new PosType*[Nmask];
	for (int j = 0; j <Nmask; j++)
		mask[j] = new PosType[Nmask];

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
		}
	}

	std::valarray<PosType> map_out(0.0, map.size());
	
	for(long i=0;i<map.size();i++){
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
	
	using std::swap;
	swap(map, map_out);

	for (int j = 0; j <Nmask; j++)
		delete[] mask[j];
	delete[] mask;
}

/**
 * \brief Draws a line between two points on the image by setting
 * the pixels equal to value.
 *
 * TODO: Could be improved by detecting if the line passes through the map
 * at all before starting.  Could also be improved by making the line fatter
 * by including neighbor points.
 */
void PixelMap::drawline(
        PosType x1[]     /// one end point of line
        ,PosType x2[]    /// other end point of line
        ,PosType value   /// value that it is set to on the map
                        ){

	PosType x[2],s1,s2,r;
	size_t index;
	PosType d = 0;

	r = sqrt( (x2[0] - x1[0])*(x2[0] - x1[0]) + (x2[1] - x1[1])*(x2[1] - x1[1]) );

	if(r==0.0){
		if(inMapBox(x1)){
			index = Utilities::IndexFromPosition(x1,Npixels,range,center);
			map[index] = value;
		}
		return;
	}

	s1 = (x2[0] - x1[0])/r;
	s2 = (x2[1] - x1[1])/r;

	x[0] = x1[0];
	x[1] = x1[1];
	while(d <= r){
		if(inMapBox(x)){
			index = Utilities::IndexFromPosition(x,Npixels,range,center);
			map[index] = value;
		}
		x[0] += s1*resolution;
		x[1] += s2*resolution;
		d += resolution;
	}

	return;
}

/**
 * \brief Draws a circle
 */
void PixelMap::drawcircle(
                        PosType r_center[]    /// center of circle
                        ,PosType radius       ///  radius of circle
                        ,PosType value        /// value that it is set to on the map
                        ){
  
	PosType x1[2],x2[2];
  PosType dtheta = resolution/fabs(radius);

  for(float theta = 0; theta < 2*pi; theta += dtheta){
    x1[0] = r_center[0] + radius*cos(theta);
    x1[1] = r_center[1] + radius*sin(theta);
    x2[0] = r_center[0] + radius*cos(theta+dtheta);
    x2[1] = r_center[1] + radius*sin(theta+dtheta);
    drawline(x1,x2,value);
  }
    
	return;
}
/**
 * \brief Draws a closed curve through the points in curve->imagekist
 *
 * This differs form PixelMap::AddImage() in that it draws lines between the points
 * and takes no account of the cells that the points are in or the surface brightness.
 * The points must be ordered already.  Particularly useful for drawing the caustics
 * that may have irregular cell sizes.  The last point will be connected to the first point.
 */
void PixelMap::AddCurve(ImageInfo *curve,PosType value){
	PosType x[2];

  if(curve->imagekist->Nunits() == 0 ) return;

	curve->imagekist->MoveToTop();
	x[0] = curve->imagekist->getCurrent()->x[0];
	x[1] = curve->imagekist->getCurrent()->x[1];
	curve->imagekist->Down();
	for(;!(curve->imagekist->OffBottom());curve->imagekist->Down()){
		drawline(x,curve->imagekist->getCurrent()->x,value);
		x[0] = curve->imagekist->getCurrent()->x[0];
		x[1] = curve->imagekist->getCurrent()->x[1];
	}
	curve->imagekist->MoveToTop();
	drawline(x,curve->imagekist->getCurrent()->x,value);

	return;
}

/**
 *  \brief Fills in pixels where the image plane points in the grid are located with the value given
 */
void PixelMap::AddGrid(Grid &grid,PosType value){
  PointList* list = grid.i_tree->pointlist;
  size_t index;

  list->current = list->top;
  do{
    if(inMapBox(list->current->x)){
      index = Utilities::IndexFromPosition(list->current->x,Npixels,range,center);
      map[index] = value;
    }
  }while(MoveDownList(list));
}


/// Find arcs in image
void PixelMap::FindArc(
                       PosType &radius
                       ,PosType *xc
                       ,PosType *arc_center
                       ,PosType &arclength
                       ,PosType &width
                       ,PosType threshold    // threshold in pixal value
                       ){
  
  std::vector<size_t> mask(Npixels*Npixels);
  size_t j=0;
  long k=0;
  PosType const tmp_center[2] = {0,0};
    
  // mask pixels below threshhold
  PosType maxval = map[0],minval = map[0];
  for(size_t i=0;i<Npixels*Npixels;i++){
    if(map[i] > threshold){
      if(j==0) minval = map[i];
      else minval = MIN(minval,map[i]);
      mask[j++]=i;
    }
    maxval = MAX(maxval,map[i]);
  }
  mask.resize(j);
  minval *= 0.99;
  
  if(j == 0){
    std::cout << "PixelMap::FindArc() - No pixels above surface brighness limit" << std::endl;
    radius = arclength = width = 0.0;
    xc[0] = xc[1] = 0.0;
    return;
  }
  
  PosType Rmax,Rmin,r2;
  Rmax = Npixels;
  Rmin = 2;
  
  size_t Nc,Nr;
  
  Nc = (size_t)(2*Rmax);
  Nr = (size_t)(Rmax-Rmin)/2;
  
  std::vector<PosType> x(Nc),y(Nc),R2(Nr);
  Utilities::D3Matrix<float> votes(Nc,Nc,Nr);
  for(size_t i = 0;i<Nc*Nc*Nr;++i) votes(i) = 0;
  
  for(size_t i = 0;i<Nc;++i) x[i] = i*2*Rmax/(Nc-1) - Rmax;
  for(size_t i = 0;i<Nc;++i) y[i] = i*2*Rmax/(Nc-1) - Rmax;
  for(size_t i = 0;i<Nr;++i) R2[i] = pow(Rmin + i*(Rmax-Rmin)/(Nr-1),2);
  
  const PosType range = 1.0*Npixels;
  PosType RmaxSqr = Rmax*Rmax;
  PosType rminSqr = RmaxSqr,rmax2=0;
  for(size_t m=0;m<mask.size();++m){
    
    Utilities::PositionFromIndex(mask[m], xc, Npixels, range, tmp_center);
    
    for(size_t ii=0;ii<Nc;++ii){
      for(size_t jj=0;jj<Nc;++jj){
        
        r2 = (xc[0]-x[ii])*(xc[0]-x[ii]) + (xc[1]-y[jj])*(xc[1]-y[jj]);
        
        rminSqr = MIN(rminSqr,r2);
        rmax2 = MAX(rmax2,r2);
        
        if(r2 < RmaxSqr){
          k = locate<PosType>(R2,r2);
          if(k > -1 && k < Nr){
            //votes(ii,jj,k) += log(map[mask[m]]/minval);
            votes(ii,jj,k) += map[mask[m]];
            //std::cout << "vote = " << votes(ii,jj,k) << std::endl;
          }
        }
      }
    }
  }
  
  printFITS("!fit_test.fits");
  
  // find maximum votes
  size_t kmax=0,ksecond=0;
  PosType maxvotes,secondvotes;
  maxvotes = votes(0);
  for(size_t kk=0;kk < Nc*Nc*Nr;++kk){
        
    if(votes(kk) >= maxvotes ){
      
      secondvotes = maxvotes;
      ksecond = kmax;
      
      maxvotes = votes(kk);
      kmax = kk;
    }
  }
  
  xc[0] = x[votes.xindex(kmax)];
  xc[1] = y[votes.yindex(kmax)];
  radius = sqrt(R2[votes.zindex(kmax)]);
  
  PosType x_tmp[2],r_tmp,rmax,rmin;
  double xave[2] = {0,0};
  rmax = 0.0;
  rmin = radius;
  arclength = 0;
  
  // find arc length, width and center
  for(size_t m=0;m<mask.size();++m){
    Utilities::PositionFromIndex(mask[m], x_tmp, Npixels, range, tmp_center);
    r_tmp = sqrt( (xc[0] - x_tmp[0])*(xc[0] - x_tmp[0]) + (xc[1] - x_tmp[1])*(xc[1] - x_tmp[1]) );
    
    if(fabs(r_tmp-radius) < 1.0){
      arclength += 1;
      xave[0] += x_tmp[0];
      xave[1] += x_tmp[1];
    }
    rmax = MAX(rmax,r_tmp);
    rmin = MIN(rmin,r_tmp);
  }
  xave[0] /= arclength;
  xave[1] /= arclength;

  double tmp = sqrt( (xave[0] - xc[0])*(xave[0] - xc[0]) + (xave[1] - xc[1])*(xave[1] - xc[1]) );
  arc_center[0] = radius*(xave[0] - xc[0])/tmp + xc[0];
  arc_center[1] = radius*(xave[1] - xc[1])/tmp + xc[1];

  // convert from pixel units to angular units
  width = MAX(rmax-rmin,1.0)*resolution;
  arclength *= resolution;
  radius *= resolution;

  xc[0] = xc[0]*resolution + center[0];
  xc[1] = xc[1]*resolution + center[1];

  arc_center[0] = arc_center[0]*resolution + center[0];
  arc_center[1] = arc_center[1]*resolution + center[1];

}
