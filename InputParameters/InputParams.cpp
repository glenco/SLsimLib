/*
 * InputParams.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: bmetcalf
 */
#include "slsimlib.h"

#ifdef ENABLE_FITS
#include <CCfits/CCfits>
#endif

namespace
{
	template<typename T>
	std::string to_str(const T& v)
	{
		std::stringstream sstr;
		sstr << v;
		return sstr.str();
	}
}

/** \brief The constructor reads in all the lines from the parameter file and stores the labels and values of each parameter.
*  If a parameter is provided that should not exist an error message is printed
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

	while(!file_in.eof())
	{
		std::string myline;
		getline(file_in, myline);
		std::cout << myline << std::endl;

		// remove all tabs from the string
		std::size_t np = myline.find('\t');
		while(np != std::string::npos)
		{
			myline.replace(np, 1, " ");
			np = myline.find('\t');
		}

		// strip of comments
		np = myline.find_first_of('#');
		if(np == std::string::npos)
			np = myline.size(); // Case where no comment is found
		if(np > 0)
		{
			std::string comment = myline.substr(np);
			myline = myline.substr(0, np);  // strip out comments

			// parse line
			np = myline.find_first_of(' ');
			if(np > 0)
			{
				std::string label = myline.substr(0, np);
				
				// check if label is known
				if(known_labels.find(label) == known_labels.end())
					std::cerr << "WARNING: unknown parameter " << label << " in file " << paramfile << "!" << std::endl;
				
				myline = myline.substr(np);  // should now just have the value plus white space
				np = myline.find_first_not_of(' ');
				if(np > 0)
					myline = myline.substr(np);  // strip off spaces after value if they exist
				np = myline.find_first_of(' ');
				myline = myline.substr(0,np);
				if(myline.size() <= 0)
				{
					std::cout << "ERROR: Paramters " << label << " does not have a valid value in parameter file " << paramfile_name << std::endl;
					exit(1);
				}

				params.insert(std::make_pair(label, myline));
				if(!comment.empty())
					comments.insert(std::make_pair(label, comment));
			}
		}

		myline.clear();
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
	print_used();
}

/// Print all parameters and values to stdout.
void InputParams::print()
{
	std::cout << "number of lines read: " << params.size() << std::endl;
	for(iterator it = params.begin(); it != params.end(); ++it)
		std::cout << it->first << "\t\t" << it->second << "\t" << comments[it->first] << std::endl;
}

/// Print parameters and values that have been accessed within the code to stdout.
void InputParams::print_used()
{
	std::cout << "###### Used Parameters #######" << std::endl;
	std::size_t n = 0;
	for(iterator it = params.begin(); it != params.end(); ++it)
	{
		if(use_number[it->first] > 0)
		{
			std::cout << it->first << "\t\t" << it->second << "\t" << comments[it->first] << std::endl;
			++n;
		}
	}
	std::cout << std::endl << n << " parameters where USED out of a total of " << params.size() << " paramaters read from the parameter file." << std::endl;
}

/// Print parameters and values that where read in but not accessed within the code to stdout.
void InputParams::print_unused()
{
	std::cout << "###### Unused Parameters #######" << std::endl;
	std::cout << "number of lines read: " << params.size() << std::endl;
	std::size_t n = 0;
	for(iterator it = params.begin(); it != params.end(); ++it)
	{
		if(use_number[it->first] == 0)
		{
			std::cout << it->first << "\t\t" << it->second << "\t" << comments[it->first] << std::endl;
			++n;
		}
	}
	std::cout << std::endl << n << " parameters where UNUSED out of a total of " << params.size() << " paramaters read from the parameter file." << std::endl;
}

/// Print all parameters to a file in the format needed for an input parameter file. Unused parameters can be stripped with the optional second argument.
void InputParams::PrintToFile(std::string filename, bool strip_unused)
{
	paramfile_name = filename;
	std::ofstream file_out(paramfile_name.c_str());

	std::cout << "Creating parameter file: " << paramfile_name;
	file_out << "# This parameter file was created by GLAMER."<< std::endl;
	file_out << "# It can be used as an input parameter file."<< std::endl;
	file_out << "# number of parameters: " << params.size() << std::endl << std::endl;
	
	for(iterator it = params.begin(); it != params.end(); ++it)
	{
		if(strip_unused && use_number[it->first] == 0)
			continue;
		
		iterator comment = comments.find(it->first);
		if(comment == comments.end())
		   file_out << it->first << "\t\t" << it->second << std::endl;
		else
			file_out << it->first << "\t\t" << it->second << "\t" << comment->second << std::endl;
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
#ifdef ENABLE_FITS
	std::cout << "Reading MOKA FITS parameters...\n" << std::endl;
	
	std::string MOKA_input_file;
	if(!get("MOKA_input_file", MOKA_input_file))
		throw new std::runtime_error("Parameter MOKA_input_file must be set for MOKA_input_params to work!");
	
	try
	{
		std::auto_ptr<CCfits::FITS> ff(new CCfits::FITS(MOKA_input_file, CCfits::Read));
		
		CCfits::PHDU* h0 = &ff->pHDU();
		
		double sidel; // box side length in arc seconds
		double zlens; // redshift of lens
		double zsource; // redshift of source
		double omega_m; // omega matter
		double omega_l; // omega_lambda
		double hubble; // hubble constant H/100
		
		h0->readKey("SIDEL", sidel);
		h0->readKey("ZLENS", zlens);
		h0->readKey("ZSOURCE", zsource);
		h0->readKey("OMEGA", omega_m);
		h0->readKey("LAMBDA", omega_l);
		h0->readKey("H", hubble);
		
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
	catch(CCfits::FITS::CantOpen)
	{
		std::cout << "can not open " << MOKA_input_file << std::endl;
		exit(1);
	}
#else
	std::cout << "Please enable the preprocessor flag ENABLE_FITS !" << std::endl;
	exit(1);
#endif
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * bool entries in the parameter file must be 0,1,true or false.
 */
bool InputParams::get(std::string label, bool& value)
{
	iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	if(!it->second.compare("0") || !it->second.compare("false"))
	{
		value = false;
		++use_number[it->first];
		return true;
	}
	
	if(!it->second.compare("1") || !it->second.compare("true"))
	{
		value = true;
		++use_number[it->first];
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
bool InputParams::get(std::string label, MassFuncType& value)
{
	iterator it = params.find(label);
	if(it == params.end())
		return false;

	if(!it->second.compare("0") || !it->second.compare("PS"))
	{
		value = PS;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("1") || !it->second.compare("ST"))
	{
		value = ST;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("2") || !it->second.compare("PowerLaw"))
	{
		value = PL;
		++use_number[it->first];
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

bool InputParams::get(std::string label, LensHaloType& value)
{
	iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	if(!it->second.compare("0") || !it->second.compare("nolens"))
	{
		value = null_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("1") || !it->second.compare("NFW"))
	{
		value = nfw_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("2") || !it->second.compare("PseudoNFW"))
	{
		value = pnfw_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("3") || !it->second.compare("PowerLaw"))
	{
		value = pl_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("4") || !it->second.compare("NSIE"))
	{
		value = nsie_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("5") || !it->second.compare("AnaLens"))
	{
		value = ana_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("6") || !it->second.compare("UniLens")){
		
		value = uni_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("7") || !it->second.compare("MOKALens"))
	{
		value = moka_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("8") || !it->second.compare("DummyLens"))
	{
		value = dummy_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("9") || !it->second.compare("Hernquist"))
	{
		value = hern_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("10") || !it->second.compare("Jaffe"))
	{
		value = jaffe_lens;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("11") || !it->second.compare("MultiDark"))
	{
		value = multi_dark_lens;
		++use_number[it->first];
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

bool InputParams::get(std::string label, GalaxyLensHaloType& value)
{
	iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	if(!it->second.compare("0") || !it->second.compare("none"))
	{
		value = null_gal;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("1") || !it->second.compare("NSIE"))
	{
		value = nsie_gal;
		++use_number[it->first];
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be 0 or none, 1 or NSIE" << std::endl;
	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MainLensType entries in the parameter file must be 0 through 2 or NFW, PowerLaw, or PointMass.
 */

bool InputParams::get(std::string label, ClumpInternal& value)
{
	iterator it = params.find(label);
	if(it == params.end())
		return false;

	if(!it->second.compare("0") || !it->second.compare("NFW"))
	{
		value = nfw;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("1") || !it->second.compare("PowerLaw"))
	{
		value = powerlaw;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("2") || !it->second.compare("PointMass"))
	{
		value = pointmass;
		++use_number[it->first];
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
 * MainLensType entries in the parameter file must be One,Mono,BrokenPowerLaw,Salpeter,SinglePowerLaw,Kroupa or Chabrier
 */

bool InputParams::get(std::string label, IMFtype& value)
{
	iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	if(!it->second.compare("One"))
	{
		value = One;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("Mono"))
	{
		value = Mono;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("BrokenPowerLaw"))
	{
		value = BrokenPowerLaw;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("Salpeter"))
	{
		value = Salpeter;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("SinglePowerLaw"))
	{
		value = SinglePowerLaw;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("Kroupa"))
	{
		value = Kroupa;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("Chabrier"))
	{
		value = Chabrier;
		++use_number[it->first];
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
 * MainLensType entries in the parameter file must be SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,i1, or i2
 */

bool InputParams::get(std::string label, Band& value)
{
	iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	if(!it->second.compare("SDSS_U"))
	{
		value = SDSS_U;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("SDSS_G"))
	{
		value = SDSS_G;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("SDSS_R"))
	{
		value = SDSS_R;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("SDSS_I"))
	{
		value = SDSS_I;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("SDSS_Z"))
	{
		value = SDSS_Z;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("J"))
	{
		value = J;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("H"))
	{
		value = H;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("Ks"))
	{
		value = Ks;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("i1"))
	{
		value = i1;
		++use_number[it->first];
		return true;
	}
	if(!it->second.compare("i2"))
	{
		value = i2;
		++use_number[it->first];
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,i1, or i2!" << std::endl;
	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.
 *
 * If there is an entry in the parameter file this function will always
 * return it in string format - no type checking.
 */
bool InputParams::get(std::string label,std::string& value)
{
	iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	value = it->second;
	++use_number[it->first];
	
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
bool InputParams::exist(std::string label)
{
	iterator it = params.find(label);
	if(it == params.end())
		return false;
	
	return true;
}
