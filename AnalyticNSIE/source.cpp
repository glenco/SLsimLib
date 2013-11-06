/*
 * source.cpp
 *
 *  Created on: Feb 8, 2012
 *      Author: mpetkova
 */

#include "slsimlib.h"
#include <typeinfo>

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
//#include <CCfits>
#endif

using namespace std;

Source::Source()
{
	source_r = 0;
	source_x[0] = 0;
	source_x[1] = 0;
	zsource = 0;
	setSBlimit_magarcsec(30.);
}

std::size_t Source::Nparams() const
{
	return 2;
}

double Source::getParam(std::size_t p) const
{
	switch(p)
	{
		case 0:
			// source x position
			return (source_x[0]/pi*180*60*60);
		case 1:
			// source y position
			return (source_x[1]/pi*180*60*60);
		default:
			throw std::invalid_argument("bad parameter index for getParam()");
	}
}

double Source::setParam(std::size_t p, double val)
{
	switch(p)
	{
		case 0:
			// source x position
			return (source_x[0] = val*pi/180/60/60);
		case 1:
			// source y position
			return (source_x[1] = val*pi/180/60/60);
		default:
			throw std::invalid_argument("bad parameter index for setParam()");
	}
}

void Source::printCSV(std::ostream&, bool header) const
{
	const std::type_info& type = typeid(*this);
	std::cerr << "Source subclass " << type.name() << " does not implement printCSV()" << std::endl;
	std::exit(1);
}

SourceUniform::SourceUniform(InputParams& params) : Source(){
	assignParams(params);
}

SourceGaussian::SourceGaussian(InputParams& params) : Source(){
	assignParams(params);
}

SourceBLR::SourceBLR(InputParams& params) : Source(){
	assignParams(params);
}

SourceBLRDisk::SourceBLRDisk(InputParams& params) : SourceBLR(params){

}

SourceBLRSph1::SourceBLRSph1(InputParams& params) : SourceBLR(params){

}

SourceBLRSph2::SourceBLRSph2(InputParams& params) : SourceBLR(params){

}

Source::~Source(){
}

SourceUniform::~SourceUniform(){
}

SourceGaussian::~SourceGaussian(){
}

SourceBLR::~SourceBLR(){
}

SourceBLRDisk::~SourceBLRDisk(){
}

SourceBLRSph1::~SourceBLRSph1(){
}

SourceBLRSph2::~SourceBLRSph2(){
}

void SourceUniform::assignParams(InputParams& params){

	if(!params.get("z_source",zsource)){
		  ERROR_MESSAGE();
		  cout << "parameter z_source needs to be set in parameter file " << params.filename() << endl;
		  exit(0);
	  }

	  printSource();
}

void SourceGaussian::assignParams(InputParams& params){


	if(!params.get("z_source",zsource)){
		  ERROR_MESSAGE();
		  cout << "parameter z_source needs to be set in parameter file " << params.filename() << endl;
		  exit(0);
	  }

	if(!params.get("gauss_r2",source_gauss_r2)){
		  ERROR_MESSAGE();
		  cout << "parameter gauss_r2 needs to be set in parameter file " << params.filename() << endl;
		  exit(0);
	  }

	  printSource();
}

void SourceBLR::assignParams(InputParams& params){

	bool fail = false;
	if(!params.get("source_z",zsource)){
		  ERROR_MESSAGE();
		  cout << "parameter source_z_source needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
	  }
	if(!params.get("source_BHmass",source_BHmass)){
		  ERROR_MESSAGE();
		  cout << "parameter source_BHmass needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
	  }

	if(!params.get("source_gamma",source_gamma)){
		  ERROR_MESSAGE();
		  cout << "parameter source_gamma needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
	  }
	if(!params.get("source_inclin",source_inclination)){
		  ERROR_MESSAGE();
		  cout << "parameter source_inclin needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
	  }
	if(!params.get("source_opening_ang",source_opening_angle)){
		  ERROR_MESSAGE();
		  cout << "parameter source_opening_ang needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
	  }
	if(!params.get("source_r_in",source_r_in)){
		  ERROR_MESSAGE();
		  cout << "parameter source_r_in needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
	  }
	if(!params.get("source_r_out",source_r_out)){
		  ERROR_MESSAGE();
		  cout << "parameter source_r_out needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
	  }
	if(!params.get("source_nuo",source_nuo)){
		  ERROR_MESSAGE();
		  cout << "parameter source_nuo needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
	  }
	if(!params.get("source_fK",source_fK)){
		  ERROR_MESSAGE();
		  cout << "parameter source_fK needs to be set in parameter file " << params.filename() << endl;
		  fail = true;
	  }

	if(fail) exit(1);

	source_inclination *= pi/180;
	source_opening_angle *= pi/180;
	source_monocrome = false;

	printSource();
}

void SourceUniform::printSource(){
	cout << endl << "**Source model**" << endl;

	cout << "z_source " << zsource << endl << endl;
}

void SourceGaussian::printSource(){
	cout << endl << "**Source model**" << endl;

	cout << "z_source " << zsource << endl;
	cout << "gauss_r2 " << source_gauss_r2 << endl << endl;
}

void SourceBLR::printSource(){
	cout << endl << "**Source model**" << endl;

	cout << "z_source " << zsource << endl;
	cout << "BHmass " << source_BHmass << endl;
	cout << "gamma " << source_gamma << endl;
	cout << "incl " << source_inclination*180/pi << endl;
	cout << "opening angl " << source_opening_angle*180/pi << endl;
	cout << "r_in " << source_r_in << endl;
	cout << "r_out " << source_r_out << endl;
	cout << "nuo " << source_nuo << endl;
	cout << "source_fK " << source_fK << endl << endl;
}

double SourceUniform::SurfaceBrightness(double *y){
	return (double)( (pow(y[0]-getX()[0],2) + pow(y[1]-getX()[1],2)) < source_r*source_r );
}

double SourceGaussian::SurfaceBrightness(double *y){
	return exp( -(pow(y[0]-getX()[0],2) + pow(y[1]-getX()[1],2))/source_gauss_r2 );
}
// surface brightness for models of the Broad Line Region
double SourceBLRDisk::SurfaceBrightness(double *y){
	double x[2] = {y[0]-getX()[0],y[1]-getX()[1]};
	return blr_surface_brightness_disk(x,this);
}

double SourceBLRSph1::SurfaceBrightness(double *y){
	return blr_surface_brightness_spherical_circular_motions(sqrt((pow(y[0]-getX()[0],2) + pow(y[1]-getX()[1],2))),this);
}
double SourceBLRSph2::SurfaceBrightness(double *y){
	return blr_surface_brightness_spherical_random_motions(sqrt((pow(y[0]-getX()[0],2) + pow(y[1]-getX()[1],2))),this);
}

//void in_source(double *y_source,ListHndl sourcelist){
//  return;
//}
SourcePixelled::SourcePixelled(InputParams& params)
{}

SourcePixelled::SourcePixelled(
		double my_z            /// redshift of the source
		, double* my_center  /// center (in rad)
		, int my_Npixels           /// number of pixels per side
		, double my_resolution  /// resolution (in rad)
		, double* arr_val          /// array of pixel values (must be of size = Npixels*Npixels)
		)
	:Source(), resolution(my_resolution), Npixels (my_Npixels){
	zsource = my_z;

	range = resolution*(Npixels-1);

	values.resize(Npixels*Npixels);
	for (int i = 0; i < Npixels*Npixels; i++)
			values[i] = arr_val[i];
	source_r =  range/sqrt(2.);
	source_x[0] = my_center[0];
	source_x[1] = my_center[1];
	calcTotalFlux();
	calcCentroid();
	calcEll();
	calcSize();
}

/** \brief Creates a SourcePixelled from a PixelMap image.
 *  The idea is to use stamps of observed galaxies as input sources for simuations.
 *  Surface brightness of the source is conserved, taking into account the input pixel size.
 *  Factor allows for rescaling of the flux, in case one wants to simulate a different observation.
 */
SourcePixelled::SourcePixelled(
		const PixelMap& gal_map  /// Input image and information
		, double my_z                 /// redshift of the source
		, double factor                /// optional rescaling factor for the flux
		)
	:Source(){

	zsource = my_z;
	resolution = gal_map.getResolution();
	Npixels = gal_map.getNpixels();
	range = resolution*(Npixels-1);
	source_x[0] = gal_map.getCenter()[0];
	source_x[1] = gal_map.getCenter()[1];
	source_r =  range/sqrt(2.);
	values.resize(Npixels*Npixels);
	for (int i = 0; i < Npixels*Npixels; i++)
			values[i] = gal_map[i]/resolution/resolution*factor;

	calcTotalFlux();
	calcCentroid();
	calcEll();
	calcSize();
}

SourcePixelled::~SourcePixelled(){
}

void SourcePixelled::calcCentroid(){
	double x_sum = 0;
	double y_sum = 0;
	double sum = 0;
	double x[2];
	for (unsigned long i = 0; i < Npixels*Npixels; i++)
	{
		Utilities::PositionFromIndex(i,x,Npixels,range,source_x);
		x_sum += x[0]*values[i];
		y_sum += x[1]*values[i];
		sum += values[i];
	}
	centroid[0] = x_sum/sum;
	centroid[1] = y_sum/sum;
}

void SourcePixelled::calcEll(){
	for (int j = 0; j < 2; j++)
	{
		for (int k = 0; k < 2; k++)
		{
			quad[j][k] = 0.;
		}
	}

	double x[2];
	for (unsigned long i = 0; i < Npixels*Npixels; i++)
	{
		Utilities::PositionFromIndex(i,x,Npixels,range,source_x);
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				quad[j][k] += values[i]*(x[j]-centroid[j])*(x[k]-centroid[k]);
			}
		}
	}

	ell[0] = (quad[0][0] - quad[1][1])/(quad[0][0] + quad[1][1]);
	ell[1] = 2*quad[0][1]/(quad[0][0] + quad[1][1]);
}

void SourcePixelled::calcSize(){
	double r_sum = 0.;
	double sum = 0.;
	double rad;
	double x[2];

	for (unsigned long i = 0; i < Npixels*Npixels; i++)
	{
		Utilities::PositionFromIndex(i,x,Npixels,range,source_x);
		rad = sqrt((x[0]-centroid[0])*(x[0]-centroid[0])+(x[1]-centroid[1])*(x[1]-centroid[1]));
		r_sum += rad*values[i];
		sum += values[i];
	}
	size = r_sum/sum;
}

double SourcePixelled::SurfaceBrightness(double *y){
	long ix = Utilities::IndexFromPosition(y[0],Npixels,range,source_x[0]);
	long iy = Utilities::IndexFromPosition(y[1],Npixels,range,source_x[1]);
	if (ix>-1 && iy>-1)
	{
		return values[iy*Npixels+ix];
	}
	else
		return 0.;
}

void SourcePixelled::calcTotalFlux(){
	double val_tot = 0.;
	for (int i = 0; i < Npixels*Npixels; i++)
		val_tot += values[i];
	flux = val_tot*resolution*resolution;
}

void SourcePixelled::printSource(){}
void SourcePixelled::assignParams(InputParams& params){}

/**  \brief Calculates the difference in magnitude when changing the observing filter
 *
 */
double Source::changeFilter(
		std::string filter_in  				/// file with the old observing filter
		, std::string filter_out			/// file with the new observing filter
		, std::string sed					/// file with the galaxy spectral energy distribution
		)
{

	// reads in the input filter
	std::ifstream fin(filter_in.c_str());
	std::vector<double> wavel_in, ampl_in;
	fin.seekg(0,fin.beg);
	double x,y;
	for (int i = 0; ; i++)
	{
		if( fin.eof() ) break;
		if (fin.peek() == '#')
		{
			fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			continue;
		}
		fin >> x;
		wavel_in.push_back(x);
		fin >> y;
		ampl_in.push_back(y);
	}

	// reads in the output filter
	std::ifstream fout(filter_out.c_str());
	std::vector<double> wavel_out, ampl_out;
	fout.seekg(0,fout.beg);
	for (int i = 0; ; i++)
	{
		if (fout.eof() ) break;
		if (fout.peek() == '#')
		{
			fout.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			continue;
		}
		fout >> x;
		wavel_out.push_back(x);
		fout >> y;
		ampl_out.push_back(y);
	}

	// reads in the source sed
	std::ifstream sed_input(sed.c_str());
	std::vector<double> wavel_sed, ampl_sed;
	sed_input.seekg(0,sed_input.beg);
	for (int i = 0; ; i++)
	{
		if(sed_input.eof() ) break;
		if (sed_input.peek() == '#')
		{
			sed_input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			continue;
		}
		sed_input >> x;
		wavel_sed.push_back(x);
		sed_input >> y;
		ampl_sed.push_back(y);
		wavel_sed[i] *= 1+zsource;
	}

	// Applies Delta m = int (sed*fout) / int (sed*fin) * int fin / int fout
 	double fin_int = integrateFilter(wavel_in,ampl_in);
	double fout_int = integrateFilter(wavel_out,ampl_out);
	double sed_in = integrateFilterSed(wavel_in, ampl_in, wavel_sed, ampl_sed);
	double sed_out = integrateFilterSed(wavel_out, ampl_out, wavel_sed, ampl_sed);
	double delta_mag = -2.5 * log10(sed_out/sed_in*fin_int/fout_int);
	return delta_mag;
}

/**  \brief Calculates the integral of the filter curve given as an array of (x,y) values.
 *
 */
double Source::integrateFilter(std::vector<double> wavel_fil, std::vector<double> fil)
{
	double integr = 0.;
	for (int i = 0; i < wavel_fil.size()-1; i++)
		integr += (wavel_fil[i+1] - wavel_fil[i])*(fil[i+1]+fil[i]);
	integr /= 2.;
	return integr;
}

/**  \brief Calculates the integral of the sed mutliplied by the filter curve.
 *
 */
double Source::integrateFilterSed(std::vector<double> wavel_fil, std::vector<double> fil, std::vector<double> wavel_sed, std::vector<double> sed)
{
	int wavel_new_size = 10000;
	vector<double> wavel_new(wavel_new_size), fil_val(wavel_new_size), sed_val(wavel_new_size);

	double integr = 0.;
	for (int i = 0; i < wavel_new_size; i++)
	{
		wavel_new[i] = max(wavel_fil[0],wavel_sed[0]) + double(i)/double(wavel_new_size-1)*(min(wavel_fil[wavel_fil.size()-1],wavel_sed[wavel_sed.size()-1])-max(wavel_fil[0],wavel_sed[0]) );
		vector<double>::iterator it_f = lower_bound(wavel_fil.begin(), wavel_fil.end(), wavel_new[i]);
		int p = it_f-wavel_fil.begin()-1;
		fil_val[i] = fil[p] + (fil[p+1]-fil[p])*(wavel_new[i]-wavel_fil[p])/(wavel_fil[p+1]-wavel_fil[p]);
		vector<double>::iterator it_s = lower_bound(wavel_sed.begin(), wavel_sed.end(), wavel_new[i]);
		int q = it_s-wavel_sed.begin()-1;
		sed_val[i] = sed[q] + (sed[q+1]-sed[q])*(wavel_new[i]-wavel_sed[q])/(wavel_sed[q+1]-wavel_sed[q]);
	}
	for (int i = 0; i < wavel_new_size-1; i++)
	{
		integr += (wavel_new[i+1] - wavel_new[i])*(fil_val[i+1]*sed_val[i+1]*wavel_new[i+1]*wavel_new[i+1]+fil_val[i]*sed_val[i]*wavel_new[i]*wavel_new[i]);
	}
	integr /= 2.;
	return integr;
}

SourceShapelets::SourceShapelets(
		double my_z                              /// redshift of the source
		, double* my_center           			/// center (in rad)
		, double my_mag							/// magnitude
		, double my_scale						/// scale of the shapelets decomposition
		, std::valarray<double> my_coeff  	/// coefficients of the shapelets decomposition
		, double my_ang					/// rotation angle (in rad)
		)
		:Source()
{
	zsource = my_z;
	mag = my_mag;
	source_x[0] = my_center[0];
	source_x[1] = my_center[1];
	ang = my_ang;
	n1 = sqrt(my_coeff.size());
	n2 = n1;
	coeff = my_coeff;
	source_r = my_scale;

	NormalizeFlux();
}

SourceShapelets::SourceShapelets(
		double my_z							/// redshift of the source
		, double* my_center					/// center (in rad)
		, double my_mag						/// magnitude
		, std::string shap_file				/// fits file with coefficients in a square array
		, double my_ang			/// rotation angle (in rad)
		)
		:Source()
{
	//TODO Fabio: add parameter for rotation

	zsource = my_z;
	mag = my_mag;
	source_x[0] = my_center[0];
	source_x[1] = my_center[1];
	ang = my_ang;

#ifdef ENABLE_FITS
	if(shap_file.empty())
		throw std::invalid_argument("Please enter a valid filename for the FITS file input");

	std::auto_ptr<CCfits::FITS> fp(new CCfits::FITS(shap_file.c_str(), CCfits::Read));

	CCfits::PHDU& h0 = fp->pHDU();

	h0.readKey("BETA", source_r);
	source_r *= 0.03/180./60./60.*pi;
	h0.readKey("DIM", n1);
	n2 = n1;
	h0.read(coeff);

#else
	std::cerr << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif

	NormalizeFlux();
}

SourceShapelets::SourceShapelets(
		double* my_center  					/// center (in rad)
		, std::string shap_file				/// fits file with coefficients in a square array. Mag and redshift are read from the header.
		, double my_ang				 /// rotation angle (in rad)
		)
		:Source()
{
	//TODO Fabio: add parameter for rotation

	source_x[0] = my_center[0];
	source_x[1] = my_center[1];
	ang = my_ang;

#ifdef ENABLE_FITS
	if(shap_file.empty())
		throw std::invalid_argument("Please enter a valid filename for the FITS file input");

	std::auto_ptr<CCfits::FITS> fp(new CCfits::FITS(shap_file.c_str(), CCfits::Read));

	CCfits::PHDU& h0 = fp->pHDU();

	h0.readKey("BETA", source_r);
	source_r *= 0.03/180./60./60.*pi;
	h0.readKey("MAG", mag);
	h0.readKey("REDSHIFT", zsource);
	h0.readKey("DIM", n1);
	n2 = n1;
	h0.read(coeff);

#else
	std::cerr << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif

	NormalizeFlux();
}


double SourceShapelets::SurfaceBrightness(double *y)
{
	double sb = 0.;
	double y_norm[2];
	y_norm[0] = ((y[0]-source_x[0])*cos(ang)+(y[1]-source_x[1])*sin(ang))/source_r;
	y_norm[1] = ((y[0]-source_x[0])*sin(ang)-(y[1]-source_x[1 ])*cos(ang))/source_r;
	double dist = sqrt(y_norm[0]*y_norm[0]+y_norm[1]*y_norm[1]);
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++)
		{
			double norm = 1./sqrt(pow(2,i+j)*pi*factrl(i)*factrl(j));
			sb += norm*coeff[j*n1+i]*Hermite(i,y_norm[0])*Hermite(j,y_norm[1]);
		}
	}
	sb *= exp(-dist*dist/2.)/source_r;
	return max(sb,std::numeric_limits<double>::epsilon());
}

/// Returns the value of the Hermite polynomial of degree n at position x
double SourceShapelets::Hermite(int n, double x)
{
	double hg[n];
	hg[0] = 1.;
	for (int i = 1; i <= n; i++)
	{
		if (i==1)
			hg[1] = 2.*x;
		else
			hg[i] = 2.*x*hg[i-1]-2.*(i-1)*hg[i-2];
	}
	return hg[n];
}

void SourceShapelets::printSource(){};
void SourceShapelets::assignParams(InputParams& params){};

/// Rescales the coefficients to make the source bright as we want.
void SourceShapelets::NormalizeFlux()
{
	double shap_flux = 0.;
	for (int i = 0; i < n1; i=i+2)
	{
		for (int j = 0; j < n2; j=j+2)
		{
			shap_flux += pow(2,0.5*(2-i-j))*sqrt(factrl(i))/factrl(i/2.)*sqrt(factrl(j))/factrl(j/2.)*coeff[j*n1+i];
		}
	}
	shap_flux *= sqrt(pi)*source_r;

	flux = pow(10,-0.4*(mag+48.6))*inv_hplanck;

	coeff *= flux/shap_flux;
}
