/*
 * source.h
 *
 *  Created on: Feb 6, 2012
 *      Author: mpetkova
 */

#ifndef SOURCE_H_
#define SOURCE_H_

#include <galaxies.h>

class Source{
public:

	  Source();
	  ~Source();

	  /// names of clump and sb models
	  typedef enum {Uniform,Gaussian,BLR_Disk,BLR_Sph1,BLR_Sph2,MultiAnaSource} SBModel;

	  // in lens.cpp
	  virtual double SurfaceBrightness(double *y) = 0;
	  virtual void readParamfile(std::string) = 0;
	  virtual void printSource() = 0;

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
	  /// In the case of a single plane lens, the ratio of angular size distances
	  virtual inline double getDlDs(){return DlDs;}
	  //TODO BEN I think this need only be in the BLR source models
	  virtual void setDlDs(double my_DlDs){DlDs = my_DlDs;}

protected:
	  SBModel source_sb_type;

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
};

typedef Source *SourceHndl;

class SourceUniform : public Source{
public:
	double SurfaceBrightness(double *y);
	void readParamfile(std::string);
	void printSource();

	SourceUniform(std::string);
	~SourceUniform();
};

class SourceGaussian : public Source{
public:
	  /// internal scale parameter
	  double source_gauss_r2;

	  double SurfaceBrightness(double *y);
	  void readParamfile(std::string);
	  void printSource();

	  SourceGaussian(std::string);
	  ~SourceGaussian();
};

class SourceBLR : public Source{
public:
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

	  void readParamfile(std::string);
	  void printSource();

	  SourceBLR(std::string);
	  ~SourceBLR();
};

class SourceBLRDisk : public SourceBLR{
public:
	double SurfaceBrightness(double *y);

	SourceBLRDisk(std::string);
	~SourceBLRDisk();
};

class SourceBLRSph1 : public SourceBLR{
public:
	double SurfaceBrightness(double *y);

	SourceBLRSph1(std::string);
	~SourceBLRSph1();
};

class SourceBLRSph2 : public SourceBLR{
public:
	double SurfaceBrightness(double *y);

	SourceBLRSph2(std::string);
	~SourceBLRSph2();
};

/// pointer to surface brightness function
//double (Source::*SurfaceBrightness)(double *y);

#endif /* SOURCE_H_ */
