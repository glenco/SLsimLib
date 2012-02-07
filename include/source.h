/*
 * source.h
 *
 *  Created on: Feb 6, 2012
 *      Author: mpetkova
 */

#ifndef SOURCE_H_
#define SOURCE_H_

class Source{
public:
	 /// names of clump and sb models
	  typedef enum {Uniform,Gaussian,BLR_Disk,BLR_Sph1,BLR_Sph2} SBModel;

	  // source parameters
	  /// lag time
	  double source_tau;
	  /// frequency
	  double source_nu;
	  /// total source size, ie no flux outside this radius
	  double source_r;
	  /// center of source
	  double source_x[2];

	  SBModel source_sb_type;

	  /// redshift of source
	  double zsource;

	  // Dl / Ds -- needed for the blr source models
	  double DlDs;

	  // in lens.cpp
	  virtual double source_sb_func(double *y) = 0;

	/*  double uniform_SB(double *y);
	  double gaussian_SB(double *y);
	  double BLR_Disk_SB(double *y);
	  double BLR_Sph1_SB(double *y);
	  double BLR_Sph2_SB(double *y);*/

	  Source();
	  ~Source();

	  virtual void readParamfile(string) = 0;
	  virtual void printSource() = 0;
};

typedef Source *SourceHndl;

class SourceUniform : public Source{
public:
	double source_sb_func(double *y);
	void readParamfile(string);
	void printSource();

	SourceUniform(string);
	~SourceUniform();
};

class SourceGaussian : public Source{
public:
	  /// internal scale parameter
	  double source_gauss_r2;

	  double source_sb_func(double *y);
	  void readParamfile(string);
	  void printSource();

	  SourceGaussian(string);
	  ~SourceGaussian();
};

class SourceBLR : public Source{
public:
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

	  void readParamfile(string);
	  void printSource();

	  SourceBLR(string);
	  ~SourceBLR();
};

class SourceBLRDisk : public SourceBLR{
public:
	double source_sb_func(double *y);

	SourceBLRDisk(string);
	~SourceBLRDisk();
};

class SourceBLRSph1 : public SourceBLR{
public:
	double source_sb_func(double *y);

	SourceBLRSph1(string);
	~SourceBLRSph1();
};

class SourceBLRSph2 : public SourceBLR{
public:
	double source_sb_func(double *y);

	SourceBLRSph2(string);
	~SourceBLRSph2();
};

/// pointer to surface brightness function
//double (Source::*source_sb_func)(double *y);

#endif /* SOURCE_H_ */
