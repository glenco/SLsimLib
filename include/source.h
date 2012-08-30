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

	  // source parameters
	  /// total source size, ie no flux outside this radius
	  double source_r;
	  /// center of source
	  double source_x[2];

	  /// redshift of source
	  double zsource;

	  // Dl / Ds -- needed for the blr source models
	  double DlDs;

	  // in lens.cpp
	  virtual double SurfaceBrightness(double *y) = 0;

	  Source();
	  ~Source();

	  virtual void readParamfile(std::string) = 0;
	  virtual void printSource() = 0;
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
