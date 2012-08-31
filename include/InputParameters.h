/*
 * InputParameters.h
 *
 *  Created on: Aug 29, 2012
 *      Author: bmetcalf
 */

#ifndef INPUTPARAMETERS_H_
#define INPUTPARAMETERS_H_

//TODO BEN/MARGARETA Complete this class and convert source and lens constructors to take to reduce redundancy.
/**
 * This class is for reading in all the parameters from the parameter file and then for
 * passing their values to the other constructors, etc.
 */
class InputParameters {
public:
	InputParameters();
	InputParameters(std::string parameter_filename);
	virtual ~InputParameters();

	int NumberOfParameters;

};

#endif /* INPUTPARAMETERS_H_ */
