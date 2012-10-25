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
class InputParams {
public:
	InputParams(std::string paramfile);
	virtual ~InputParams();

	bool get(std::string label,bool& value);
	bool get(std::string label,std::string& value);
	bool get(std::string label,int& value);
	bool get(std::string label,double& value);
	bool get(std::string label,IntProfType& value);

	bool exist(std::string label);
	void print();
	/// Returns total number of parameters.
	unsigned int Nparams(){return labels.size();}
	void PrintToFile(std::string filename);

private:
	std::string paramfile_name;

	std::vector<std::string> labels;
	std::vector<std::string> char_values;

	double string_to_double( const std::string& s,bool& numeric);
	int string_to_int( const std::string& s,bool& numeric);
};

#endif /* INPUTPARAMS_H_ */
