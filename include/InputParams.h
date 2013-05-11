/*
 * InputParams.h
 *
 *  Created on: Oct 17, 2012
 *      Author: bmetcalf
 */


#ifndef INPUTPARAMS_H_
#define INPUTPARAMS_H_

/**
 * This class is for reading in and writing parameter files and storing the parameters to be passed around the code.
 *
 * In the parameter file comments can be added with the comment character is #.  Any information
 * on a line after a # is stripped out.
 *
 * A valid parameter must have a label that starts a new line followed by at least one space
 * and then the value of the parameter followed by a comment (with #) if desired.
 *
 * The get() methods are used to retrieve parameter values and test if they exist.
 */

#include "standard.h"

typedef enum {PowerLaw, NFW, PseudoNFW, NSIE, NFW_NSIE} IntProfType;
typedef enum {PS, ST, PL} MassFuncType;
typedef enum {null, ana_lens, moka_lens} InputLens;
/// names of clump and sb models
typedef enum {nfw,powerlaw,pointmass} ClumpInternal;
enum IMFtype {One,Mono,BrokenPowerLaw,Salpeter,SinglePowerLaw,Kroupa,Chabrier};
enum Band {SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,i1,i2};


class InputParams {
public:
	InputParams(std::string paramfile);
	~InputParams();

	bool get(std::string label,bool& value);
	bool get(std::string label,std::string& value);
	bool get(std::string label,int& value);
	bool get(std::string label,double& value);
	bool get(std::string label,float& value);
	bool get(std::string label,IntProfType& value);
	bool get(std::string label,MassFuncType& value);
	bool get(std::string label,InputLens& value);
	bool get(std::string label,ClumpInternal& value);
	bool get(std::string label,Band& value);
	bool get(std::string label,IMFtype& value);

	bool exist(std::string label);
	void print();
	void print_used();
	void print_unused();
	/// Returns total number of parameters.
	unsigned int Nparams(){return labels.size();}
	void PrintToFile(std::string filename);
	void PrintUsedToFile(std::string filename);
	/// Return name of the parameter file.
	std::string filename(){return paramfile_name;}

	void put(std::string label,std::string value,std::string comment = "");
	void put(std::string label,int value,std::string comment = "");
	void put(std::string label,float value,std::string comment = "");
	void put(std::string label,double value,std::string comment = "");

private:
	std::string paramfile_name;

	std::vector<std::string> labels;
	std::vector<std::string> char_values;
	std::vector<std::string> comments;
	// The number of times a parameter was retrieved by get
	std::vector<int> use_number;

	double string_to_double( const std::string& s,bool& numeric);
	int string_to_int( const std::string& s,bool& numeric);
};

#endif /* INPUTPARAMS_H_ */
