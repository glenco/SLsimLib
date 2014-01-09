#include "../include/delens.h"

#include <random>
#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <cmath>

#include "../include/Tree.h"

PixelMap delens(Lens& lens, const PixelMap& lensed, std::size_t Nrays, double res_factor)
{
	// number of points to shoot per batch, to not overflow memory
	constexpr std::size_t Nbatch = 4096;
	
	// the resulting image should match the original
	PixelMap delensed = ((res_factor == 1) ? PixelMap(lensed) : PixelMap(lensed, res_factor));
	
	// clear the image
	delensed.Clean();
	
	// random number generator
	std::mt19937 rng((std::random_device())());
	
	// uniform random number generator for positioning within a pixel
	std::uniform_real_distribution<double> unif(-0.5*lensed.getResolution(), 0.5*lensed.getResolution());
	
	// inverse magnifications for each pixel
	std::vector<double> inv_mag(lensed.size());
	
	// calculate magnification for each pixel
	{
		// number of pixels
		std::size_t Npixels = lensed.size();
		
		// current pixel index
		std::size_t i = 0;
		
		// do in batches
		while(Npixels > 0)
		{
			// number of points in batch
			std::size_t Npoints = std::min(Nbatch, Npixels);
			
			// points of each pixel
			Point* points = NewPointArray(Npoints, true);
			
			// set the positions of each pixel
			for(std::size_t j = 0; j < Npoints; ++j)
				Utilities::PositionFromIndex(i+j, points[j].x, lensed.getNpixels(), lensed.getRange(), lensed.getCenter());
			
			// points in source plane
			Point* spoints = LinkToSourcePoints(points, Npoints);
			
			// shoot rays to get magnifications
			lens.rayshooterInternal(Npoints, points, false);
			
			// store magnification and sign for each pixel
			for(std::size_t j = 0; j < Npoints; ++j, ++i)
				inv_mag[i] = std::abs(points[j].invmag);
			
			// free arrays
			FreePointArray(points);
			FreePointArray(spoints);
			
			// decrese number of pixels
			Npixels -= Npoints;
		}
	}
	
	// delens image
	{
		// current pixel index
		std::size_t i = 0;
		
		// current ray index
		std::size_t j = 0;
		
		// run in batches
		while(true)
		{
			// create an array of points for each ray
			Point* points = NewPointArray(Nbatch, true);
			
			// current point index
			std::size_t Npoints = 0;
			
			// set positions for rays of each pixel
			for(; i < lensed.size(); ++i)
			{
				// center of the pixel
				double pos[2];
				Utilities::PositionFromIndex(i, pos, lensed.getNpixels(), lensed.getRange(), lensed.getCenter());
				
				// number of rays in the pixel
				std::size_t Npix = std::ceil(inv_mag[i]*Nrays);
				
				// create random ray positions within pixel
				for(; j < Npix && Npoints < Nbatch; ++j, ++Npoints)
				{
					points[Npoints].x[0] = pos[0] + unif(rng);
					points[Npoints].x[1] = pos[1] + unif(rng);
				}
				
				// check if batch is full
				if(Npoints == Nbatch)
					break;
				
				// reset ray index if done with pixel
				if(j == Npix)
					j = 0;
			}
			
			// create source points
			Point* spoints = LinkToSourcePoints(points, Npoints);
			
			// shoot rays
			lens.rayshooterInternal(Npoints, points, false);
			
			// create an image of the delensed rays
			for(std::size_t k = 0; k < Npoints; ++k)
			{
				// get index for lensed pixel
				long l = Utilities::IndexFromPosition(points[k].x, lensed.getNpixels(), lensed.getRange(), lensed.getCenter());
				
				// get index for delensed pixel
				long d = Utilities::IndexFromPosition(spoints[k].x, delensed.getNpixels(), delensed.getRange(), delensed.getCenter());
				
				// number of rays for pixel
				std::size_t Npix = std::ceil(inv_mag[l]*Nrays);
				
				// if the pixel is in the image, add one count
				if(d != -1)
					delensed.AddValue(d, std::abs(points[k].invmag)*lensed[l]/Npix);
			}
			
			// delete point arrays
			FreePointArray(points);
			if(spoints)
				FreePointArray(spoints);
			
			// stop if all pixels and rays are done
			if(i == lensed.size() && j == 0)
				break;
		}
	}
	
	return delensed;
}
