#include "../include/image_likelihood.h"

#include <algorithm>

#include "grid_maintenance.h"
#include "map_images.h"

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

double ImageLikelihood::operator()(Model& model)
{
	// create a grid for the model
	Grid grid(model.lens, (int)grid_points, grid_center, grid_range);
	
	// number of images generated
	int image_count;
	
	// do the mapping
	map_images(
		model.lens, // Lens* lens
		model.source, // Source* source
		&grid, // Grid* grid
		&image_count, // int* Nimages
		images, // ImageInfo* imageinfo
		(int)images_size, // int Nimagesmax
		model.source->getRadius(), // double xmax
		0.1*model.source->getRadius(), // double xmin
		0, // double initial_size
		EachImage, // ExitCriterion criterion
		true, // bool kappa_off
		false, // bool FindCenter
		true // bool divide_images
	);
	
	// build the fit image
	// TODO: loop for multi source
	PixelMap image(dta.getCenter(), dta.getNpixels(), dta.getResolution());
	image.AddImages(images, image_count, false);
	
	// return multinormal log-likelihood
	return -0.5*dta.chi_square(image);
}

void ImageLikelihood::redof()
{
	dof = std::pow(dta.getNpixels(), 2) - dim;
}
