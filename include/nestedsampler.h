/*
 * nestedsampler.h
 *
 *  Created on: Mar 7, 2012
 *      Author: bmetcalf
 */

class Sample{
	double *params;
	double lnL;  // log likelihood
	double lnW;  // log of weights some of which is the evidence
};

class NestedSampler{
public:
	NestedSampler(double (*lnlikelihood)(double *),double (*lnprior)(double *),void (*stepper)(double *),int Nparams);
	~NestedSampler();

	void explore(int Nsamples);
	void PrintResults();

	// output
	double getEvidence();
	double getInformation();
	double getSample(std::vector<Sample> sample);
	double getBestFit(double *param);

private:

	double Z;
	double H;

	int Nparameters;  // Number of parameters, dimension of parameter space
	double (*lnlikelihood)(double *params);  //
	double (*lnprior)(double *params);
	void (*stepper)(double *params);

};
