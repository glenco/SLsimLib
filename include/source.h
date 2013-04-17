/*
 * source.h
 *
 *  Created on: Feb 6, 2012
 *      Author: mpetkova
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include "standard.h"
#include "galaxies.h"
#include "InputParams.h"

class Source{
public:

	  Source();
	  virtual ~Source();

	  /// names of clump and sb models
	  typedef enum {Uniform,Gaussian,BLR_Disk,BLR_Sph1,BLR_Sph2,MultiAnaSource,Pixelled,Sersic} SBModel;

	  // in lens.cpp
	  /// Surface brightness of source in grid coordinates not source centered coordinates.
	  virtual double SurfaceBrightness(double *y) = 0;
	  virtual double getTotalFlux() = 0;
	  virtual void printSource() = 0;
	  // TODO Fabio: What are the units?
	  double getSBlimit(){return sb_limit;}


	  // accessor functions that will sometimes be over ridden in class derivatives
	  /// Redshift of source
	  virtual inline double getZ(){return zsource;}
	  virtual void setZ(double my_z){zsource = my_z;}
	  /// Radius of source TODO units?
	  virtual inline double getRadius(){return source_r;}
	  virtual void setRadius(double my_radius){source_r = my_radius;}
	  /// position of source TODO units?
	  virtual inline double* getX(){return source_x;}
	  virtual inline void setX(double *xx){source_x[0] = xx[0]; source_x[1] = xx[1];}
	  void setX(double my_x,double my_y){source_x[0] = my_x; source_x[1] = my_y;}
	  /// In the case of a single plane lens, the ratio of angular size distances
	  virtual inline double getDlDs(){return DlDs;}
	  //TODO BEN I think this need only be in the BLR source models
	  virtual void setDlDs(double my_DlDs){DlDs = my_DlDs;}
	  void setSBlimit(float limit) {sb_limit = limit;}

protected:
	  SBModel source_sb_type;
	  virtual void assignParams(InputParams& params) = 0;

	  // source parameters
	  /// total source size, ie no flux outside this radius
	  double source_r;
	  /// center of source
	  double source_x[2];

	  /// redshift of source
	  double zsource;
	  /// Dl / Ds -- needed for the blr source models
	  //TODO Could this be moved into the BLR classes because they are the only ones that use it.
	  double DlDs;
	  double sb_limit;
};

typedef Source *SourceHndl;

class SersicSource: public Source{
public:
	SersicSource();
	SersicSource(double mag,double Reff,double PA,double my_index,double my_q,double my_z=0,const double *theta=0);
	~SersicSource();

	void setInternals(double mag,double Reff,double PA,double my_index,double my_q,double my_z=0,const double *theta=0);

	inline double getMag() { return mag; }
	// approximation of the radius that encloses 99% of the flux
	inline double getRadius() { return (3.73 - 0.926*index + 1.164*index*index)*Reff; }
	inline double getPA() { return PA; }

	double SurfaceBrightness(double *x);
	inline double getTotalFlux() {return flux;}
	void printSource();

private:
	void assignParams(InputParams& params);
	double Reff;
	double mag;
	double PA;
	double index;
	double bn;
	double q;
	double Ieff;
	double flux;
};


class PixelledSource: public Source{
public:
	PixelledSource(double my_z, int Npixels, double range, double* center, double* arr_val);
	PixelledSource(InputParams& params);
	~PixelledSource();
	double SurfaceBrightness(double *y);
	void printSource();
	inline double getTotalFlux(){return flux;}
	inline double getRadius(){return source_r;}
	inline double* getEll(){return ell;};
	inline double getQuad(int i, int j){return quad[i][j];};
	inline double getSize(){return size;};
	inline double* getCentroid(){return centroid;};
	inline double getMag(){return -2.5*log10(flux)-48.6;};
private:
	void assignParams(InputParams& params);
	void calcEll();
	void calcSize();
	void calcCentroid();
	void calcTotalFlux();
	double resolution;
	double range;
	long Npixels;
	double flux;
	double quad[2][2];
	double ell[2];
	double size;
	double centroid[2];
	std::valarray<float> values;
};

class SourceUniform : public Source{
public:
	SourceUniform(InputParams& params);
	~SourceUniform();

	double SurfaceBrightness(double *y);
	void assignParams(InputParams& params);
	void printSource();
	double getTotalFlux(){return pi*source_r*source_r;}
};

class SourceGaussian : public Source{
public:
	  SourceGaussian(InputParams& params);
	  ~SourceGaussian();

	  /// internal scale parameter
	  double source_gauss_r2;

	  double SurfaceBrightness(double *y);
	  void assignParams(InputParams& params);
	  void printSource();
	  double getTotalFlux(){std::cout << "No total flux in SourceGaussian yet" << std::endl; exit(1);}
};

class SourceBLR : public Source{
public:

	  SourceBLR(InputParams& params);
	  ~SourceBLR();

	  void printSource();
	  double getTotalFlux(){std::cout << "No total flux in SourceBLR yet" << std::endl; exit(1);}

	  virtual inline double getRadius(){return source_r_out;}

	/// lag time
	  double source_tau;
	  /// frequency
	  double source_nu;
	  float source_nuo;
	  /// inner radius of BLR
	  float source_r_in;
	  /// outer radius of BLR
	  float source_r_out;
	  ///inclination of BLR in radians, face on is
	  float source_inclination;
	  float source_opening_angle;
	  float source_gamma;
	  float source_BHmass;
	  /// fraction of Keplerian velocity in random motions
	  float source_fK;
	  /// set to true to integrate over frequency
	  bool source_monocrome;
private:
	  void assignParams(InputParams& params);
};

class SourceBLRDisk : public SourceBLR{
public:
	double SurfaceBrightness(double *y);
	double getTotalFlux(){std::cout << "No total flux in SourceBLRDisk yet" << std::endl; exit(1);}

	SourceBLRDisk(InputParams&);
	~SourceBLRDisk();
};

class SourceBLRSph1 : public SourceBLR{
public:
	double SurfaceBrightness(double *y);
	double getTotalFlux(){std::cout << "No total flux in SourceBLRSph1 yet" << std::endl; exit(1);}

	SourceBLRSph1(InputParams&);
	~SourceBLRSph1();
};

class SourceBLRSph2 : public SourceBLR{
public:
	double SurfaceBrightness(double *y);
	double getTotalFlux(){std::cout << "No total flux in SourceBLRSph2 yet" << std::endl; exit(1);}

	SourceBLRSph2(InputParams&);
	~SourceBLRSph2();
};

/// pointer to surface brightness function
//double (Source::*SurfaceBrightness)(double *y);

#endif /* SOURCE_H_ */
