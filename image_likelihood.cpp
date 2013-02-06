#include "slsimlib.h"
#include <iostream>
#include <math.h>
#include <cstdlib>
using namespace std;

double chi_square (PixelMap data, PixelMap model, float background_subtract, float background_noise, float norm, PixelMap* mask)
{
	int data_Npixels = data.getNpixels();
	int model_Npixels = model.getNpixels();
	if(data_Npixels != model_Npixels){
		std::cout << "Size of data and model do not agree!" << std::endl;
		exit(1);
	}
    int relevant_npixels = 0;
	if (mask != NULL)
	{
		int mask_Npixels = mask->getNpixels();
		if(data_Npixels != mask_Npixels){
			std::cout << "Size of data and mask do not agree!" << std::endl;
			exit(1);
		}
	}
	bool good_pix[data_Npixels*data_Npixels];
    double chi = 0.;
    double diff;
    double sigma;
	for (int i = 0; i < data_Npixels*data_Npixels; i++)
	{
		if (mask == NULL || mask->getValue(i) != 0)
		{
			diff = norm*(data.getValue(i)-model.getValue(i)-background_subtract);
			sigma = sqrt(norm*(data.getValue(i)+background_noise));
			chi += diff*diff/sigma/sigma;
			relevant_npixels ++;
		}
	}
	return chi/relevant_npixels;
}
