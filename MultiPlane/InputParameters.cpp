/*
 * InputParameters.cpp
 *
 *  Created on: Aug 29, 2012
 *      Author: bmetcalf
 */

#include <slsimlib.h>
#include <InputParameters.h>

/// Reads in the parameter file and stores the parameters internally
InputParameters::InputParameters(std::string parameter_filename) {
	// TODO Finish writing this
	  std::cout << "analytic lens: reading from " << parameter_filename << std::endl;

	  std::ifstream file_in(parameter_filename.c_str());
	  if(!file_in){
	    std::cout << "Can't open file " << parameter_filename << std::endl;
	    exit(1);
	  }
}

InputParameters::~InputParameters() {
}
