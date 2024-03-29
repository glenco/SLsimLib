/*
 * InputParams.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: bmetcalf
 */
//#include "slsimlib.h"
#include <iomanip>

#include "cpfits.h"
#include "InputParams.h"

namespace
{
	template<typename T>
	std::string to_str(const T& v)
	{
		std::stringstream sstr;
		sstr << v;
		return sstr.str();
	}
	
	template<typename Stream, typename ValueType>
	void printrow(Stream& str, const std::string& label, const ValueType& value)
	{
		str << std::left << std::setw(23) << label << " " << value << std::endl;
	}
	
	template<typename Stream, typename ValueType>
	void printrow(Stream& str, const std::string& label, const ValueType& value, const std::string& comment)
	{
		str << std::left << std::setw(23) << label << " " << std::setw(11) << value << " # " << comment << std::endl;
	}
	
	// whitespace characters
	const std::string WS = " \t\r\n";
}

InputParams::counter::counter()
{
}

InputParams::counter::counter(const InputParams::counter& other)
: c(other.c)
{
}

InputParams::counter& InputParams::counter::operator=(const InputParams::counter& rhs)
{
	c = rhs.c;
	return *this;
}

void InputParams::counter::use(const std::string& label)
{
	std::lock_guard<std::mutex> hold(mutex);
	
	++c[label];
}

bool InputParams::counter::is_used(const std::string& label) const
{
	std::map<std::string, std::size_t>::const_iterator it = c.find(label);
	if(it == c.end())
		return false;
	return (it->second > 0);
}

/**
 * \brief The constructor creates an empty map of input params.
 */
InputParams::InputParams()
{
}

/**
 * \brief The constructor reads in all the lines from the parameter file and stores the labels and values of each parameter.
 * 
 * If a parameter is provided that should not exist an error message is printed
 */
InputParams::InputParams(std::string paramfile)
{
	std::cout << "Reading parameters from file: " << paramfile << std::endl;
	std::ifstream file_in(paramfile.c_str());
	if(!file_in)
	{
		std::cout << "Can't open parameter file " << paramfile << std::endl;
		exit(1);
	}
	paramfile_name = paramfile;

	for(std::string myline; std::getline(file_in, myline);)
	{
		// string leading whitespace if present
		std::size_t np = myline.find_first_not_of(WS);
		if(np != std::string::npos)
			myline = myline.substr(np);
		
		// look for comment and store it
		std::string comment;
		np = myline.find_first_of('#');
		if(np != std::string::npos)
			comment = myline.substr(np);
		
		// strip comment
		myline = myline.substr(0, np);
		
		// get label
		np = myline.find_first_of(WS);
		std::string label = myline.substr(0, np);
		
		// skip empty lines
		if(label.empty())
			continue;
		
		// check if label is known
		if(!is_known(label))
			std::cerr << "WARNING: unknown parameter '" << label << "' in file " << paramfile << "!" << std::endl;
		
		// skip over label
		myline = myline.substr(np);
		
		// skip whitespace after label if present
		np = myline.find_first_not_of(WS);
		if(np != std::string::npos)
			myline = myline.substr(np);
		
		// store value
		np = myline.find_first_of(WS);
		myline = myline.substr(0, np);
		
		// value is necessary
		if(myline.empty())
		{
			std::cout << "ERROR: Paramters " << label << " does not have a valid value in parameter file " << paramfile_name << std::endl;
			exit(1);
		}
		
		// now store parameter, value and optional comment
		params.insert(std::make_pair(label, myline));
		if(!comment.empty())
			comments.insert(std::make_pair(label, comment));
	}

	std::cout << "number of lines read: " << params.size() << std::endl;
	//print();
	
	// check if MOKA parameters should be read
	bool read_moka = false;
	if(get("MOKA_input_params", read_moka) && read_moka)
		readMOKA();
}

InputParams::~InputParams()
{
	print_unused();
	//print_used();
}

/// Print all parameters and values to stdout.
void InputParams::print() const
{
	std::cout << "number of lines read: " << params.size() << std::endl;
	for(const_iterator it = params.begin(); it != params.end(); ++it)
	{
		const_iterator comment = comments.find(it->first);
		if(comment != comments.end())
			printrow(std::cout, it->first, it->second, comment->second);
		else
      printrow(std::cout, it->first, it->second);
	}
}

/// Print parameters and values that have been accessed within the code to stdout.
void InputParams::print_used() const
{
	std::cout << "###### Used Parameters #######" << std::endl;
	std::size_t n = 0;
	for(const_iterator it = params.begin(); it != params.end(); ++it)
	{
		if(use_counter.is_used(it->first))
		{
			const_iterator comment = comments.find(it->first);
			if(comment != comments.end())
				printrow(std::cout, it->first, it->second, comment->second);
			else
				printrow(std::cout, it->first, it->second);
			++n;
		}
	}
	std::cout << std::endl << n << " parameters where USED out of a total of " << params.size() << " paramaters read from the parameter file." << std::endl;
}

/// Print parameters and values that where read in but not accessed within the code to stdout.
void InputParams::print_unused() const
{
	std::cout << "###### Unused Parameters #######" << std::endl;
	std::cout << "number of lines read: " << params.size() << std::endl;
	std::size_t n = 0;
	for(const_iterator it = params.begin(); it != params.end(); ++it)
	{
		if(!use_counter.is_used(it->first))
		{
			const_iterator comment = comments.find(it->first);
			if(comment != comments.end())
				printrow(std::cout, it->first, it->second, comment->second);
			else
				printrow(std::cout, it->first, it->second);
			++n;
		}
	}
	std::cout << std::endl << n << " parameters where UNUSED out of a total of " << params.size() << " paramaters read from the parameter file." << std::endl;
}

std::ostream &operator<<(std::ostream &os, InputParams const &p) {
  std::size_t n = 0;
  os << "Used :" << std::endl;
  for(InputParams::const_iterator it = p.params.begin(); it != p.params.end(); ++it)
  {
    if(p.use_counter.is_used(it->first))
    {
      InputParams::const_iterator comment = p.comments.find(it->first);
      if(comment != p.comments.end())
        printrow(os, it->first, it->second, comment->second);
      else
        printrow(os, it->first, it->second);
      ++n;
    }
  }


  os << "Unused :" << std::endl;
 for(InputParams::const_iterator it = p.params.begin(); it != p.params.end(); ++it)
  {
    if(!p.use_counter.is_used(it->first))
    {
      InputParams::const_iterator comment = p.comments.find(it->first);
      if(comment != p.comments.end())
        printrow(os, it->first, it->second, comment->second);
      else
        printrow(os, it->first, it->second);
      ++n;
    }
  }
  return os;
}

/// Print all parameters to a file in the format needed for an input parameter file. Unused parameters can be stripped with the optional second argument.
void InputParams::PrintToFile(std::string filename, bool strip_unused) const
{
	std::ofstream file_out(filename.c_str());

	std::cout << "Creating parameter file: " << filename;
	file_out << "# This parameter file was created by GLAMER."<< std::endl;
	file_out << "# It can be used as an input parameter file."<< std::endl;
	file_out << "# number of parameters: " << params.size() << std::endl << std::endl;
	
	for(const_iterator it = params.begin(); it != params.end(); ++it)
	{
		if(strip_unused && !use_counter.is_used(it->first))
			continue;
		
		const_iterator comment = comments.find(it->first);
		if(comment != comments.end())
			printrow(file_out, it->first, it->second, comment->second);
		else
			printrow(file_out, it->first, it->second);
	}
}

/** \brief Read input parameters from a MOKA FITS header.
 * This function reads the FITS file given by the MOKA_input_file parameter and
 * uses the values found in its header to set various input parameters such as
 * z_lens, z_source, Omega_matter, ...
 *
 * This method needs ENABLE_FITS to be defined.
 */
void InputParams::readMOKA()
{
	std::cout << "Reading MOKA FITS parameters...\n" << std::endl;
	
	std::string MOKA_input_file;
	if(!get("MOKA_input_file", MOKA_input_file))
		throw new std::runtime_error("Parameter MOKA_input_file must be set for MOKA_input_params to work!");
	
  CPFITS_READ cpfits("MOKA_input_file");
  
		double sidel; // box side length in arc seconds
		double zlens; // redshift of lens
		double zsource; // redshift of source
		double omega_m; // omega matter
		double omega_l; // omega_lambda
		double hubble; // hubble constant H/100
		
    
    int error = 0;
    
    error += cpfits.readKey("SIDEL", sidel);
    error += cpfits.readKey("ZLENS", zlens);
    error += cpfits.readKey("ZSOURCE", zsource);
    error += cpfits.readKey("OMEGA", omega_m);
    error += cpfits.readKey("LAMBDA", omega_l);
    error += cpfits.readKey("H", hubble);
    
    if(error != 0){
      std::cerr << "MOKA fits file must have header keywords:" << std::endl
      << " SIDEL - length on a side" << std::endl
      << " ZLENS - redshift of lens" << std::endl
      << " ZSOURCE - redshift of source" << std::endl
      << " OMEGA - Omega matter" << std::endl
      << " LAMBDA - Omega lambda" << std::endl
      << " H - hubble constant" << std::endl;
      exit(1);
    }

		double fov = 1.4142*1.4142*sidel*sidel;
		
		params["field_fov"] = to_str(fov);
		comments["field_fov"] = "# [MOKA]";
		std::cout << std::left << std::setw(24) << "field_fov" << std::setw(12) << fov << "# [MOKA]" << std::endl;
		
		params["z_lens"] = to_str(zlens);
		comments["z_lens"] = "# [MOKA]";
		std::cout << std::left << std::setw(24) << "z_lens" << std::setw(12) << zlens << "# [MOKA]" << std::endl;
		
		params["z_source"] = to_str(zsource);
		comments["z_source"] = "# [MOKA]";
		std::cout << std::left << std::setw(24) << "z_source" << std::setw(12) << zsource << "# [MOKA]" << std::endl;
		
		params["Omega_matter"] = to_str(omega_m);
		comments["Omega_matter"] = "# [MOKA]";
		std::cout << std::left << std::setw(24) << "Omega_matter" << std::setw(12) << omega_m << "# [MOKA]" << std::endl;
		
		params["Omega_lambda"] = to_str(omega_l);
		comments["Omega_lambda"] = "# [MOKA]";
		std::cout << std::left << std::setw(24) << "Omega_lambda" << std::setw(12) << omega_l << "# [MOKA]" << std::endl;
		
		params["hubble"] = to_str(hubble);
		comments["hubble"] = "# [MOKA]";
		std::cout << std::left << std::setw(24) << "hubble" << std::setw(12) << hubble << "# [MOKA]" << std::endl;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * bool entries in the parameter file must be 0,1,true or false.
 */
bool InputParams::get(std::string label, bool& value) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	use_counter.use(it->first);
	
	if(!it->second.compare("0") || !it->second.compare("false"))
	{
		value = false;
		return true;
	}
	
	if(!it->second.compare("1") || !it->second.compare("true"))
	{
		value = true;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be 1 or 0 representing true or false!" << std::endl;
	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MassFuncType entries in the parameter file must be 0 through 2 or PS (Press & Schechter), ST (Sheth & Torman) or PowerLaw (Power-law).
 */
bool InputParams::get(std::string label, MassFuncType& value) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;

	use_counter.use(it->first);
	
	if(!it->second.compare("0") || !it->second.compare("PS"))
	{
		value = MassFuncType::PressSchechter;
		return true;
	}
	if(!it->second.compare("1") || !it->second.compare("ST"))
	{
		value = MassFuncType::ShethTormen;
		return true;
	}
	if(!it->second.compare("2") || !it->second.compare("PowerLaw"))
	{
		value = MassFuncType::PowerLaw;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be 0, 1 or 2 or PS, ST or PowerLaw!" << std::endl;
	return false;
}
/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MainLensType entries in the parameter file must needs to be 0 or nolens, 1 or NFW, 2 or PseudoNFW,
 * 3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens
 */

bool InputParams::get(std::string label, LensHaloType& value) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	use_counter.use(it->first);
	
	if(!it->second.compare("0") || !it->second.compare("nolens"))
	{
		value = LensHaloType::null_lens;
		return true;
	}
	if(!it->second.compare("1") || !it->second.compare("NFW"))
	{
		value = LensHaloType::nfw_lens;
		return true;
	}
	if(!it->second.compare("2") || !it->second.compare("PseudoNFW"))
	{
		value = LensHaloType::pnfw_lens;
		return true;
	}
	if(!it->second.compare("3") || !it->second.compare("PowerLaw"))
	{
		value = LensHaloType::pl_lens;
		return true;
	}
	if(!it->second.compare("4") || !it->second.compare("NSIE"))
	{
		value = LensHaloType::nsie_lens;
		return true;
	}
	if(!it->second.compare("5") || !it->second.compare("AnaLens"))
	{
		value = LensHaloType::ana_lens;
		return true;
	}
	if(!it->second.compare("6") || !it->second.compare("UniLens")){
		
		value = LensHaloType::uni_lens;
		return true;
	}
	if(!it->second.compare("7") || !it->second.compare("MOKALens"))
	{
		value = LensHaloType::moka_lens;
		return true;
	}
	if(!it->second.compare("8") || !it->second.compare("DummyLens"))
	{
		value = LensHaloType::dummy_lens;
		return true;
	}
	if(!it->second.compare("9") || !it->second.compare("Hernquist"))
	{
		value = LensHaloType::hern_lens;
		return true;
	}
	if(!it->second.compare("10") || !it->second.compare("Jaffe"))
	{
		value = LensHaloType::jaffe_lens;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be 0 or nolens, 1 or NFW, 2 or PseudoNFW,"
		"3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens, 9 or Hernquist, 10 or Jaffe" << std::endl;
	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * GalaxyLensType entries in the parameter file must be 0 or none, 1 or NSIE
 */

bool InputParams::get(std::string label, GalaxyLensHaloType& value) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	use_counter.use(it->first);
	
	if(!it->second.compare("0") || !it->second.compare("none"))
	{
		value = GalaxyLensHaloType::null_gal;
		return true;
	}
	if(!it->second.compare("1") || !it->second.compare("NSIE"))
	{
		value = GalaxyLensHaloType::nsie_gal;
		return true;
	}
	if(!it->second.compare("2") || !it->second.compare("PowerLaw"))
	{
		value = GalaxyLensHaloType::pl_gal;
		return true;
	}
	if(!it->second.compare("3") || !it->second.compare("Hernquist"))
	{
		value = GalaxyLensHaloType::hern_gal;
		return true;
	}
  if(!it->second.compare("4") || !it->second.compare("Jaffe"))
	{
		value = GalaxyLensHaloType::jaffe_gal;
		return true;
	}
  
  
	std::cout << label << " in parameter file " << paramfile_name << " needs to be 0 or none, 1 or NSIE, 2 or PowerLaw, 3 or Hernquist, 4 or Jaffe" << std::endl;

	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MainLensType entries in the parameter file must be 0 through 2 or NFW, PowerLaw, or PointMass.
 */

bool InputParams::get(std::string label, ClumpInternal& value) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;

	use_counter.use(it->first);
	
	if(!it->second.compare("0") || !it->second.compare("NFW"))
	{
		value = ClumpInternal::nfw;
		return true;
	}
	if(!it->second.compare("1") || !it->second.compare("PowerLaw"))
	{
		value = ClumpInternal::powerlaw;
		return true;
	}
	if(!it->second.compare("2") || !it->second.compare("PointMass"))
	{
		value = ClumpInternal::pointmass;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be 0, 1 or 2 or NFW, PowerLaw, or PointMass!" << std::endl;
	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MainLensType entries in the parameter file must be 0 through 2 or NFW, PowerLaw, or PointMass.
 */

bool InputParams::get(std::string label, HaloCatFormats& value) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;
  
	use_counter.use(it->first);
	
	if(!it->second.compare("MillenniumObs"))
	{
		value = HaloCatFormats::MillenniumObs;
		return true;
	}
  if(!it->second.compare("MultiDarkHalos"))
  {
    value = HaloCatFormats::MultiDarkHalos;
    return true;
  }
  if(!it->second.compare("ObservedData"))
  {
    value = HaloCatFormats::ObservedData;
    return true;
  }
  
	std::cout << label << " in parameter file " << paramfile_name << " needs to be MillenniumObs, MultiDarkHalos or ObservedData!" << std::endl;
	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MainLensType entries in the parameter file must be One,Mono,BrokenPowerLaw,Salpeter,SinglePowerLaw,Kroupa or Chabrier
 */

bool InputParams::get(std::string label, IMFtype& value) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	use_counter.use(it->first);
	
	if(!it->second.compare("One"))
	{
		value = IMFtype::One;
		return true;
	}
	if(!it->second.compare("Mono"))
	{
		value = IMFtype::Mono;
		return true;
	}
	if(!it->second.compare("BrokenPowerLaw"))
	{
		value = IMFtype::BrokenPowerLaw;
		return true;
	}
	if(!it->second.compare("Salpeter"))
	{
		value = IMFtype::Salpeter;
		return true;
	}
	if(!it->second.compare("SinglePowerLaw"))
	{
		value = IMFtype::SinglePowerLaw;
		return true;
	}
	if(!it->second.compare("Kroupa"))
	{
		value = IMFtype::Kroupa;
		return true;
	}
	if(!it->second.compare("Chabrier"))
	{
		value = IMFtype::Chabrier;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be One, Mono, BrokenPowerLaw, Salpeter, SinglePowerLaw, Kroupa or Chabrier!" << std::endl;
	return false;
}


/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MainLensType entries in the parameter file must be Fourier,Pseudo,Schramm or Keeton
 */

bool InputParams::get(std::string label, EllipMethod& value) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	use_counter.use(it->first);
	
  if(!it->second.compare("0") || !it->second.compare("Fourier"))
	{
		value = EllipMethod::Fourier;
		return true;
	}
	if(!it->second.compare("1") || !it->second.compare("Pseudo"))
	{
		value = EllipMethod::Pseudo;
		return true;
	}
	if(!it->second.compare("2") || !it->second.compare("Schramm"))
	{
		value = EllipMethod::Schramm;
		return true;
	}
	if(!it->second.compare("3") || !it->second.compare("Keeton"))
	{
		value = EllipMethod::Keeton;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be Fourier, Pseudo, Schramm or Keeton!" << std::endl;
	return false;
}


/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MainLensType entries in the parameter file must be EUC_VIS,EUC_Y,EUC_J,EUC_H,SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,IRAC1,IRAC2,F435W,F606W,F775W,F850LP,F814W,F110W,F160W
 */

bool InputParams::get(std::string label, Band& value) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	use_counter.use(it->first);
	
	if(!it->second.compare("EUC_VIS"))
	{
		value = Band::EUC_VIS;
		return true;
	}
	if(!it->second.compare("EUC_Y"))
	{
		value = Band::EUC_Y;
		return true;
	}
	if(!it->second.compare("EUC_J"))
	{
		value = Band::EUC_J;
		return true;
	}
	if(!it->second.compare("EUC_H"))
	{
		value = Band::EUC_H;
		return true;
	}
	if(!it->second.compare("SDSS_U"))
	{
		value = Band::SDSS_U;
		return true;
	}
	if(!it->second.compare("SDSS_G"))
	{
		value = Band::SDSS_G;
		return true;
	}
	if(!it->second.compare("SDSS_R"))
	{
		value = Band::SDSS_R;
		return true;
	}
	if(!it->second.compare("SDSS_I"))
	{
		value = Band::SDSS_I;
		return true;
	}
	if(!it->second.compare("SDSS_Z"))
	{
		value = Band::SDSS_Z;
		return true;
	}
	if(!it->second.compare("J"))
	{
		value = Band::J;
		return true;
	}
	if(!it->second.compare("H"))
	{
		value = Band::H;
		return true;
	}
	if(!it->second.compare("Ks"))
	{
		value = Band::Ks;
		return true;
	}
	if(!it->second.compare("IRAC1"))
	{
		value = Band::IRAC1;
		return true;
	}
	if(!it->second.compare("IRAC2"))
	{
		value = Band::IRAC2;
		return true;
	}
	if(!it->second.compare("F435W"))
	{
		value = Band::F435W;
		return true;
	}
	if(!it->second.compare("F606W"))
	{
		value = Band::F606W;
		return true;
	}
	if(!it->second.compare("F775W"))
	{
		value = Band::F775W;
		return true;
	}
	if(!it->second.compare("F850LP"))
	{
		value = Band::F850LP;
		return true;
	}
	if(!it->second.compare("F814W"))
	{
		value = Band::F814W;
		return true;
	}
	if(!it->second.compare("F110W"))
	{
		value = Band::F110W;
		return true;
	}
	if(!it->second.compare("F160W"))
	{
		value = Band::F160W;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be EUC_VIS,EUC_Y,EUC_J,EUC_H,SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,IRAC1,IRAC2,F435W,F606W,F775W,F850LP,F814W,F110W or F160W!" << std::endl;
	return false;
}



/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.
 *
 * If there is an entry in the parameter file this function will always
 * return it in string format - no type checking.
 */
bool InputParams::get(std::string label,std::string& value) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	use_counter.use(it->first);
	
	value = it->second;
	
	return true;
}

/**
 * Add a new parameter to the parameter list.
 */
void InputParams::put(std::string label, std::string value, std::string comment)
{
	params[label] = value;
	
	if(comment.empty())
		comments.erase(label);
	else
		comments[label] = "# " + comment;

	return;
}

// Check to see if parameter exists.
bool InputParams::exist(std::string label) const
{
	const_iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	return true;
}

std::string to_string(const Band &band){
  
  switch (band) {
    case Band::EUC_VIS:
      return "EUC_VIS";
    case Band::EUC_Y:
      return "EUC_Y";
    case Band::EUC_J:
      return "EUC_J";
    case Band::EUC_H:
      return "EUC_H";
    case Band::SDSS_U:
      return "SDSS_U";
    case Band::SDSS_G:
      return "SDSS_G";
    case Band::SDSS_R:
      return "SDSS_R";
    case Band::SDSS_I:
      return "SDSS_I";
    case Band::SDSS_Z:
      return "SDSS_Z";
    case Band::J:
      return "J";
    case Band::H:
      return "H";
    case Band::Ks:
      return "Ks";
    case Band::IRAC1:
      return "IRAC1";
    case Band::IRAC2:
      return "IRAC2";
    case Band::F435W:
      return "F435W";
    case Band::F606W:
      return "F606W";
    case Band::F775W:
      return "F775W";
    case Band::F850LP:
      return "F850LP";
    case Band::F814W:
      return "F814W";
    case Band::F110W:
      return "F110W";
    case Band::F160W:
      return "F160W";
    default:
      return "UnknownBand";
      break;
  }
}

std::ostream &operator<<(std::ostream &os, Band const &band) {
  
  switch (band) {
    case Band::EUC_VIS:
      return os << "EUC_VIS";
    case Band::EUC_Y:
      return os << "EUC_Y";
    case Band::EUC_J:
      return os << "EUC_J";
    case Band::EUC_H:
      return os << "EUC_H";
    case Band::SDSS_U:
      return os << "SDSS_U";
    case Band::SDSS_G:
      return os << "SDSS_G";
    case Band::SDSS_R:
      return os << "SDSS_R";
    case Band::SDSS_I:
      return os << "SDSS_I";
    case Band::SDSS_Z:
      return os << "SDSS_Z";
    case Band::J:
      return os << "J";
    case Band::H:
      return os << "H";
    case Band::Ks:
      return os << "Ks";
    case Band::IRAC1:
      return os << "IRAC1";
    case Band::IRAC2:
      return os << "IRAC2";
    case Band::F435W:
      return os << "F435W";
    case Band::F606W:
      return os << "F606W";
    case Band::F775W:
      return os << "F775W";
    case Band::F850LP:
      return os << "F850LP";
    case Band::F814W:
      return os << "F814W";
    case Band::F110W:
      return os << "F110W";
    case Band::F160W:
      return os << "F160W";
    default:
      return os << "UnknownBand";
      break;
  }

}

