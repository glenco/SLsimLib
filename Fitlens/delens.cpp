#include "../include/delens.h"

#include <random>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <cmath>

#include "../include/Tree.h"
/**
* \brief Attempts to create a map of the source given a map of the images and
* a lens model.
*
* TODO: This needs comments.  Also it should take a random seed so that it is reproducable.
*/
PixelMap delens(Lens& lens, const PixelMap& lensed, std::size_t Nrays, double res_factor, double mag_min)
{
	// number of points to shoot per batch, to not overflow memory
	constexpr std::size_t Nbatch = 4096;
	
	// the resulting image should match the original, with optional resolution change
	PixelMap delensed = ((res_factor == 1) ? PixelMap(lensed) : PixelMap(lensed, res_factor));
	
	// clear the delensed image
	delensed.Clean();
	
	// inverse pixel area of lensed image
	double a_l = lensed.getResolution()*lensed.getResolution();
	
	// pixel area of delensed image
	double a_d = delensed.getResolution()*delensed.getResolution();
	
	// maximum inverse magnification
	double inv_mag_max = (mag_min == 0) ? 0 : std::abs(1./mag_min);
	
	// random number generator
	std::mt19937 rng((std::random_device())());
	
	// uniform random number generator for positioning within a pixel
	std::uniform_real_distribution<double> unif(-0.5*lensed.getResolution(), 0.5*lensed.getResolution());
	
	// inverse magnifications for each lensed pixel
	std::vector<double> inv_mag(lensed.size());
	
	// calculate magnification for each pixel
	{
		// points at center each pixel
		Point* points = NewPointArray(Nbatch);
		
		// points in source plane
		Point* spoints = LinkToSourcePoints(points, Nbatch);
		
		// number of pixels
		std::size_t Npixels = lensed.size();
		
		// current pixel index
		std::size_t i = 0;
		
		// do in batches
		while(Npixels > 0)
		{
			// number of points in batch
			std::size_t Npoints = std::min(Nbatch, Npixels);
			
			// set the positions of each pixel
			for(std::size_t j = 0; j < Npoints; ++j)
				Utilities::PositionFromIndex(i+j, points[j].x, lensed.getNpixels(), lensed.getRange(), lensed.getCenter());
			
			// shoot rays to get magnifications
			lens.rayshooterInternal(Npoints, points, false);
			
			// store absolute magnification for each pixel
			for(std::size_t j = 0; j < Npoints; ++j, ++i)
				inv_mag[i] = std::abs(points[j].invmag);
			
			// decrese number of pixels
			Npixels -= Npoints;
		}
		
		// free arrays
		FreePointArray(points);
		FreePointArray(spoints);
	}
	
	// delens image
	{
		// normalization for Monte Carlo integrals
		std::vector<double> norm(delensed.size(), 0);
		
		// create an array of points for each ray
		Point* points = NewPointArray(Nbatch);
		
		// create source points
		Point* spoints = LinkToSourcePoints(points, Nbatch);
		
		// current pixel index
		std::size_t i = 0;
		
		// current ray index
		std::size_t j = 0;
		
		// run in batches
		while(true)
		{
			// current point index
			std::size_t Npoints = 0;
			
			// set positions for rays of each pixel
			for(; i < lensed.size(); ++i)
			{
				// number of rays in the pixel
				std::size_t Npix = std::max(1., std::ceil(inv_mag[i]*Nrays)) + 0.5;
				
				// rays for pixel
				for(; j < Npix && Npoints < Nbatch; ++j, ++Npoints)
				{
					// center of the pixel
					Utilities::PositionFromIndex(i, points[Npoints].x, lensed.getNpixels(), lensed.getRange(), lensed.getCenter());
					
					// randomize position within pixel
					points[Npoints].x[0] += unif(rng);
					points[Npoints].x[1] += unif(rng);
				}
				
				// check if batch is full
				if(Npoints == Nbatch)
					break;
				
				// reset ray index if done with pixel
				if(j == Npix)
					j = 0;
			}
			
			// shoot rays
			lens.rayshooterInternal(Npoints, points, false);
			
			// add the delensed rays to output image
			for(std::size_t k = 0; k < Npoints; ++k)
			{
				// check maximum magnitude
				if(inv_mag_max && std::abs(points[k].invmag) > inv_mag_max)
					continue;
				
				// get index for lensed pixel
				long l = Utilities::IndexFromPosition(points[k].x, lensed.getNpixels(), lensed.getRange(), lensed.getCenter());
				
				// get index for delensed pixel
				long d = Utilities::IndexFromPosition(spoints[k].x, delensed.getNpixels(), delensed.getRange(), delensed.getCenter());
				
				// check if the delensed pixel is in the image
				if(d == -1)
					continue;
				
				// number of rays in the pixel
				std::size_t Npix = std::max(1., std::ceil(inv_mag[l]*Nrays)) + 0.5;
				
				// importance weight
				double w = std::abs(points[k].invmag)/Npix*a_l;
				
				// add to sum
				delensed.AddValue(d, w*lensed[l]/a_l);
				
				// add to normalization
				norm[d] += w;
			}
			
			// stop if all pixels and rays are done
			if(i == lensed.size() && j == 0)
				break;
		}
		
		// delete point arrays
		FreePointArray(points);
		FreePointArray(spoints);
		
		// normalize the Monte Carlo integrals
		for(std::size_t i = 0; i < delensed.size(); ++i)
			if(norm[i] > 0)
				delensed.AssignValue(i, a_d*delensed[i]/norm[i]);
	}
	
	return delensed;
}
