/*
 * InputParams.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: bmetcalf
 */
#include "slsimlib.h"

/// The constructor reads in all the lines from the parameter file and stores the labels and values of each parameter.
InputParams::InputParams(std::string paramfile) {

	std::cout << "Reading parameters from file: " << paramfile << std::endl;
	std::ifstream file_in(paramfile.c_str());
	if(!file_in){
		std::cout << "Can't open parameter file " << paramfile << std::endl;
		exit(1);
	}

	paramfile_name = paramfile;

	std::string myline,comment;
	int np,i=0;

	while(!file_in.eof()){

		getline(file_in,myline);
		std::cout << myline << std::endl;

		// remove all tabs from the string
		np = myline.find('\t');
		while(np >= 0){
			myline.replace(np,1," ");
			np = myline.find('\t');
		}

		// strip of comments
		np = myline.find_first_of('#');
		if(np == std::string::npos) np = myline.size(); // Case where no comment is found
		if(np > 0){
			comment = myline.substr(np);
			myline = myline.substr(0,np);  // strip out comments

			//pars line
			np = myline.find_first_of(' ');
			if(np > 0){
				labels.push_back(myline.substr(0,np));
				comments.push_back(comment);

				myline = myline.substr(np);  // should now just have the value plus white space
				np = myline.find_first_not_of(' ');
				if(np > 0) myline = myline.substr(np);  // strip off spaces after value if they exist
				np = myline.find_first_of(' ');
				myline = myline.substr(0,np);
				if(myline.size() <= 0){
					std::cout << "ERRROR: Paramters " << labels[i] << " does not have a valid value in parameter file " << paramfile_name << std::endl;
					exit(1);
				}

				char_values.push_back(myline);
				//std::cout << labels[i] << "     " << char_values[i] << std::endl;
				++i;
			}
		}

		myline.clear();
	}

	use_number.resize(labels.size(),0);
	std::cout << "number of lines read: " << labels.size() << std::endl;
	//print();
}

InputParams::~InputParams() {
	labels.clear();
	char_values.clear();
	comments.clear();
}

/// Print all parameters and values to stdout.
void InputParams::print(){
	std::cout << "number of lines read: " << labels.size() << std::endl;
	for(int i=0;i<labels.size();++i){
		std::cout << labels[i] << "     " << char_values[i] << "    " << comments[i] << std::endl;
	}
}

/// Print parameters and values that have been accessed within the code to stdout.
void InputParams::print_used(){
	std::cout << "###### Used Parameters #######" << std::endl;
	std::cout << "number of lines read: " << labels.size() << std::endl;
	int n=0;
	for(int i=0;i<labels.size();++i){
		if(use_number[i] > 0){
			std::cout << labels[i] << "     " << char_values[i] << "    " << comments[i] << std::endl;
			++n;
		}
	}
	std::cout << std::endl << n << " Parameters where used out of a total of " << labels.size() << " paramaters read from the parameter file. " << std::endl;
}

/// Print parameters and values that where read in but not accessed within the code to stdout.
void InputParams::print_unused(){
	std::cout << "###### Unused Parameters #######" << std::endl;
	int n=0;
	for(int i=0;i<labels.size();++i){
		if(use_number[i] == 0){
			std::cout << labels[i] << "     " << char_values[i] << "    " << comments[i] << std::endl;
			++n;
		}
	}
	std::cout << std::endl << n << " Parameters where UNUSED out of a total of " << labels.size() << " paramaters read from the parameter file. " << std::endl;
}

/// Print all parameters to a file in the format needed for an input parameter file
void InputParams::PrintToFile(std::string filename){

	paramfile_name = filename;
	std::ofstream file_out(paramfile_name.c_str());

	std::cout << "Creating parameter file: " << paramfile_name;
	file_out << "# This parameter file was created by GLAMER."<< std::endl;
	file_out << "# It can be used as an input parameter file."<< std::endl;
	file_out << "# number of parameters: " << labels.size() << std::endl << std::endl;
	for(int i=0;i<labels.size();++i){
		file_out << labels[i] << "               " << char_values[i] << "         " << comments[i] << std::endl;
	}
}

/** \brief Print used parameters to a file in the format needed for an input parameter file.
* Only the parameters that were accessed at least once are printed.
 */
void InputParams::PrintUsedToFile(std::string filename){

	paramfile_name = filename;
	std::ofstream file_out(paramfile_name.c_str());

	std::cout << "Creating parameter file: " << paramfile_name;
	file_out << "# This parameter file was created by GLAMER."<< std::endl;
	file_out << "# It can be used as an input parameter file."<< std::endl;
	file_out << "# number of parameters: " << labels.size() << std::endl << std::endl;
	for(int i=0;i<labels.size();++i){
		if(use_number[i] > 0) file_out << labels[i] << "               " << char_values[i] << "         " << comments[i] << std::endl;
	}
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * bool entries in the parameter file must be 0,1,true or false.
 */
bool InputParams::get(std::string label,bool& value){
	unsigned int i;
	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	if(!char_values[i].compare("0") || !char_values[i].compare("false")){
		value = false;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("1") || !char_values[i].compare("true")){
		value = true;
		use_number[i]++;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be 1 or 0 representing true or false!"<< std::endl;
	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MassFuncType entries in the parameter file must be 0 through 2 or PS (Press & Schechter), ST (Sheth & Torman) or PowLaw (Power-law).
 */
bool InputParams::get(std::string label,MassFuncType& value){
	unsigned int i;
	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	if(!char_values[i].compare("0") || !char_values[i].compare("PS")){
		value = PS;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("1") || !char_values[i].compare("ST")){
		value = ST;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("2") || !char_values[i].compare("PowerLaw")){
		value = PL;
		use_number[i]++;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be 0, 1 or 2 or PS, ST or PowerLaw!"<< std::endl;
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

bool InputParams::get(std::string label,LensHaloType& value){
	unsigned int i;
	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	if(!char_values[i].compare("0") || !char_values[i].compare("nolens")){
		value = null_lens;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("1") || !char_values[i].compare("NFW")){
		value = nfw_lens;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("2") || !char_values[i].compare("PseudoNFW")){
		value = pnfw_lens;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("3") || !char_values[i].compare("PowerLaw")){
		value = pl_lens;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("4") || !char_values[i].compare("NSIE")){
		value = nsie_lens;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("5") || !char_values[i].compare("AnaLens")){
		value = ana_lens;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("6") || !char_values[i].compare("UniLens")){
		value = uni_lens;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("7") || !char_values[i].compare("MOKALens")){
		value = moka_lens;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("8") || !char_values[i].compare("DummyLens")){
		value = dummy_lens;
		use_number[i]++;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be 0 or nolens, 1 or NFW, 2 or PseudoNFW, " << std::endl;
	std::cout << "3 or PowerLaw, 4 or NSIE, 5 or AnaLens, 6 or UniLens, 7 or MOKALens, 8 or DummyLens"<< std::endl;
	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * GalaxyLensType entries in the parameter file must be 0 or none, 1 or NSIE
 */

bool InputParams::get(std::string label,GalaxyLensHaloType& value){
	unsigned int i;
	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	if(!char_values[i].compare("0") || !char_values[i].compare("none")){
		value = null_gal;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("1") || !char_values[i].compare("NSIE")){
		value = nsie_gal;
		use_number[i]++;
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

bool InputParams::get(std::string label,ClumpInternal& value){
	unsigned int i;
	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	if(!char_values[i].compare("0") || !char_values[i].compare("NFW")){
		value = nfw;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("1") || !char_values[i].compare("PowerLaw")){
		value = powerlaw;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("2") || !char_values[i].compare("PointMass")){
		value = pointmass;
		use_number[i]++;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be 0, 1 or 2 or NFW, PowerLaw, or PointMass!"<< std::endl;
	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MainLensType entries in the parameter file must be One,Mono,BrokenPowerLaw,Salpeter,SinglePowerLaw,Kroupa or Chabrier
 */

bool InputParams::get(std::string label,IMFtype& value){
	unsigned int i;
	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	if(!char_values[i].compare("One")){
			value = One;
			use_number[i]++;
			return true;
	}
	if(!char_values[i].compare("Mono")){
		value = Mono;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("BrokenPowerLaw")){
		value = BrokenPowerLaw;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("Salpeter")){
		value = Salpeter;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("SinglePowerLaw")){
		value = SinglePowerLaw;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("Kroupa")){
		value = Kroupa;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("Chabrier")){
		value = Chabrier;
		use_number[i]++;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be One, Mono, BrokenPowerLaw, Salpeter, SinglePowerLaw, Kroupa or Chabrier!"<< std::endl;
	return false;
}


/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * MainLensType entries in the parameter file must be SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,i1, or i2
 */

bool InputParams::get(std::string label,Band& value){
	unsigned int i;
	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	if(!char_values[i].compare("SDSS_U")){
		value = SDSS_U;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("SDSS_G")){
		value = SDSS_G;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("SDSS_R")){
		value = SDSS_R;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("SDSS_I")){
		value = SDSS_I;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("SDSS_Z")){
		value = SDSS_Z;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("J")){
		value = J;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("H")){
		value = H;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("Ks")){
		value = Ks;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("i1")){
		value = i1;
		use_number[i]++;
		return true;
	}
	if(!char_values[i].compare("i2")){
		value = i2;
		use_number[i]++;
		return true;
	}

	std::cout << label << " in parameter file " << paramfile_name << " needs to be SDSS_U,SDSS_G,SDSS_R,SDSS_I,SDSS_Z,J,H,Ks,i1, or i2!"<< std::endl;
	return false;
}

/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.
 *
 * If there is an entry in the parameter file this function will always
 * return it in string format - no type checking.
 */
bool InputParams::get(std::string label,std::string& value){
	unsigned int i;
	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	value = char_values[i];
	use_number[i]++;
	return true;
}
/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * The entry in the parameter file must have a numerical value.  If it
 * does not, an exception message is printed and false is returned.
 */
bool InputParams::get(std::string label,double& value){
	unsigned int i;
	bool numeric;

	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	value = string_to_double(char_values[i],numeric);
	if(!numeric){
		std::cout << "Expected a numerical value for " << labels[i] << " in parameter file " << paramfile_name << std::endl;
		return false;
	}
	use_number[i]++;
	return true;
}
/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 * The entry in the parameter file must have a numerical value.  If it
 * does not, an exception message is printed and false is returned.
 */
bool InputParams::get(std::string label,float& value){
	unsigned int i;
	bool numeric;

	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	value = string_to_double(char_values[i],numeric);
	if(!numeric){
		std::cout << "Expected a numerical value for " << labels[i] << " in parameter file " << paramfile_name << std::endl;
		return false;
	}
	use_number[i]++;
	return true;
}
/** \brief Assigns to "value" the value of the parameter called "label".
 * If this parameter label does not appear in the parameter file false
 * is returned.  If the parameter in the file does not "match" the type
 * of value false will also be returned and a warning printed to stdout.
 *
 ** The entry in the parameter file must have a numerical value.  If it
 * doesn't an exception message is printed and false is returned.  If
 * the entry is a decimal number it is truncated to an integer and true
 * is returned.  This is a possible source of error if the wrong type is
 * used.
 */
bool InputParams::get(std::string label,int& value){
	unsigned int i;
	bool numeric;

	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	value = string_to_int(char_values[i],numeric);
	if(!numeric){
		std::cout << "Expected a numerical value for " << labels[i] << " in parameter file " << paramfile_name << std::endl;
		return false;
	}
	use_number[i]++;
	return true;
}
/**
 * Add a new parameter to the parameter list.
 */
void InputParams::put(std::string label,std::string value,std::string comment){

	labels.push_back(label);
	use_number.push_back(0);
	char_values.push_back(value);
	if(comment == "") comments.push_back("");
	else comments.push_back("#" + comment);

	return;
}
/**
 * Add a new parameter to the parameter list.
 */
void InputParams::put(std::string label,int value,std::string comment){

	labels.push_back(label);
	use_number.push_back(0);

	// convert to string
	std::ostringstream strs;
	strs << value;
	char_values.push_back(strs.str());
	if(comment == "") comments.push_back("");
	else comments.push_back("#" + comment);

	return;
}
/**
 * Add a new parameter to the parameter list.
 */
void InputParams::put(std::string label,float value,std::string comment){

	labels.push_back(label);
	use_number.push_back(0);

	// convert to string
	std::ostringstream strs;
	strs << value;
	char_values.push_back(strs.str());
	if(comment == "") comments.push_back("");
	else comments.push_back("#" + comment);

	return;
}
/**
 * Add a new parameter to the parameter list.
 */
void InputParams::put(std::string label,double value,std::string comment){

	labels.push_back(label);
	use_number.push_back(0);

	// convert to string
	std::ostringstream strs;
	strs << value;
	char_values.push_back(strs.str());
	if(comment == "") comments.push_back("");
	else comments.push_back("#" + comment);

	return;
}

// Check to see if parameter exists.
bool InputParams::exist(std::string label){
	unsigned int i;
	for(i=0;i<labels.size();++i)
		if(labels[i] == label) break;
	if(i==labels.size()) return false;

	return true;
}

double InputParams::string_to_double( const std::string& s,bool& numeric ){
	std::istringstream i(s);

	double x;
	if (!(i >> x)){
		numeric = false;
		return 0;
	}
	numeric = true;
	return x;
}
int InputParams::string_to_int( const std::string& s,bool& numeric  ){
	std::istringstream i(s);

	int x;
	if (!(i >> x)){
		numeric = false;
		return 0;
	}
	numeric = true;
	return x;
}
