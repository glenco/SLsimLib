#include "slsimlib.h"
#include "image_likelihood.h"

#include <algorithm>

ImageLikelihood::ImageLikelihood(PixelData data)
: dta(data),
dof(0), dim(0),
images_size(100), images(0),
grid_points(64)
{
	// create ImageInfo array
	images = new ImageInfo[images_size];
	
	// set grid to match data
	grid_range = dta.getRange();
	std::copy(dta.getCenter(), dta.getCenter() + 2, grid_center);
	
	// calculate dof
	redof();
}

ImageLikelihood::~ImageLikelihood()
{
	delete[] images;
}

unsigned long ImageLikelihood::dimension() const
{
	return dim;
}

void ImageLikelihood::dimension(unsigned long n)
{
	dim = n; redof();
}

PixelData ImageLikelihood::data() const
{
	return dta;
}

void ImageLikelihood::data(PixelData data)
{
	// copy data
	swap(dta, data);
	
	// change grid to match data
	grid_range = dta.getRange();
	std::copy(dta.getCenter(), dta.getCenter() + 2, grid_center);
	
	// recalculate degrees of freedom
	redof();
}

std::size_t ImageLikelihood::imagesSize() const
{
	return images_size;
}

void ImageLikelihood::imagesSize(std::size_t size)
{
	images_size = size;
	
	delete[] images;
	images = new ImageInfo[images_size];
}

std::size_t ImageLikelihood::gridPoints() const
{
	return grid_points;
}

void ImageLikelihood::gridPoints(std::size_t n)
{
	grid_points = n;
}

void ImageLikelihood::redof()
{
	dof = std::pow(double(dta.getNpixels()), 2) - dim;
}
