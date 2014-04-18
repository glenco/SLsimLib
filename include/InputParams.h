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
#include <map>

#if __cplusplus >= 201103L
#include <mutex>
#endif

enum MassFuncType
{
	PS,
	ST,
	PL
};

enum LensHaloType
{
	null_lens,
	nfw_lens,
	pnfw_lens,
	pl_lens,
	nsie_lens,
	ana_lens,
	uni_lens,
	moka_lens,
	dummy_lens,
	hern_lens,
	jaffe_lens,
	multi_dark_lens
};

enum GalaxyLensHaloType
{
	null_gal,
	nsie_gal
};

/// names of clump and sb models
typedef enum {nfw,powerlaw,pointmass} ClumpInternal;
/// Initial mass function type
enum IMFtype {One,Mono,BrokenPowerLaw,Salpeter,SinglePowerLaw,Kroupa,Chabrier};
/// Photometric bands
enum Band {SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,i1,i2};
/// Photometric bands for shapelets galaxies
enum shap_band {MAG_B, MAG_V, MAG_I, MAG_Z, MAG_J, MAG_H, MAG_u_KIDS, MAG_g_KIDS, MAG_r_KIDS, MAG_i_KIDS};

/** \brief Structure for reading and writing parameters to and from a parameter file as well as a container 
 * for passing the parameters to other classes and functions.
 *
 *  The constructor reads in all parameters in the parameter file.  They can then be accessed with the 
 *  get functions.  There should be no interaction with the parameter file except through the InputParams 
 *  structure.
 */
class InputParams
{
public:
	InputParams();
	InputParams(std::string paramfile);
	~InputParams();

	bool get(std::string label,bool& value) const;
	bool get(std::string label,std::string& value) const;
	bool get(std::string label,MassFuncType& value) const;
	bool get(std::string label,LensHaloType& value) const;
	bool get(std::string label,GalaxyLensHaloType& value) const;
	bool get(std::string label,ClumpInternal& value) const;
	bool get(std::string label,Band& value) const;
	bool get(std::string label,shap_band& value) const;
	bool get(std::string label,IMFtype& value) const;
	
	// get numbers
	template<typename Number>
	bool get(std::string label, Number& value) const;
	
	template<typename T>
	T get(std::string label, const T& def = T()) const;
	
	template<typename T>
	T require(std::string label) const;
	
	bool exist(std::string label) const;
	void print() const;
	void print_used() const;
	void print_unused() const;
	/// Returns total number of parameters.
	std::size_t Nparams() const { return params.size(); }
	void PrintToFile(std::string filename, bool strip_unused = false) const;
	/// Return name of the parameter file.
	std::string filename() const{return paramfile_name;}

	void put(std::string label,std::string value,std::string comment = std::string());
	
	template<typename Number>
	void put(std::string label,Number value,std::string comment = std::string());
	
	// read parameters from a MOKA FITS header
	void readMOKA();
	
	/** return sample input parameters */
	static InputParams sample();

	/** add a new parameter with default value and comment */
	static void add(std::string label, std::string value = std::string(), std::string comment = std::string());
	
private:
	// thread-safe counter
	class counter
	{
	public:
		counter();
		counter(const counter& other);
		
		counter& operator=(const counter& rhs);
		
		void use(const std::string& label);
		bool is_used(const std::string& label) const;
		
		friend void swap(counter&, counter&);
		
	private:
		std::map<std::string, std::size_t> c;
		
		mutable std::mutex mutex;
	};
	
	typedef std::map<std::string, std::string>::iterator iterator;
	typedef std::map<std::string, std::string>::const_iterator const_iterator;
	
	static bool is_known(const std::string& name);
	
	std::string paramfile_name;
	
	std::map<std::string, std::string> params;
	std::map<std::string, std::string> comments;
	
	// The number of times a parameter was retrieved by get
	mutable counter use_counter;
};

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned and the value is set to the default value of the numerical type.
 * 
 * The entry in the parameter file must have a numerical value.  If it
 * does not, an exception message is printed and false is returned.
 */
template<typename Number>
bool InputParams::get(std::string label, Number& value) const
{
	const_iterator it = params.find(label);
	
	if(it == params.end())
	{
		value = Number();
		return false;
	}
	
	std::istringstream sstr(it->second);
	
	sstr >> value;
	if(sstr.fail())
	{
		std::cout << "Expected a numerical value for " << it->first << " in parameter file " << paramfile_name << std::endl;
		return false;
	}
	
	use_counter.use(it->first);
	
	return true;
}

/**
 * Add a new parameter to the parameter list.
 */
template<typename Number>
void InputParams::put(std::string label, Number value, std::string comment)
{
	// convert to string
	std::ostringstream sstr;
	sstr << value;
	params.insert(std::make_pair(label, sstr.str()));
	
	if(comment.empty())
		comments.erase(label);
	else
		comments.insert(std::make_pair(label, comment));
}

/**
 * Return a parameter value.
 */
template<typename T>
T InputParams::get(std::string label, const T& def) const
{
	T value;
	if(!get(label, value))
		return def;
	return value;
}

/**
 * Return a parameter value or throw an error.
 */
template<typename T>
T InputParams::require(std::string label) const
{
	// TODO: add getting of default value
	T value;
	if(!get(label, value))
		throw std::runtime_error("Parameter " + label + " is required in " + paramfile_name);
	return value;
}

#endif /* INPUTPARAMS_H_ */
