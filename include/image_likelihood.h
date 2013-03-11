#ifndef IMAGE_LIKELIHOOD_H_
#define IMAGE_LIKELIHOOD_H_

#include "parameters.h"
#include "image_processing.h"

#include <algorithm>
#include <cstddef>
#include <cmath>

class Source;
class Lens;

/**
 * \brief Calculate the likelihood for a source and lens based on image data.
 * 
 * This class takes a source and lens, generates a PixelMap by simulating the
 * lensing, and compares it to a given data PixelMap. It uses the chi_square
 * function to calculate a likelihood.
 * 
 * To render the image PixelMap from source and lens, an internal grid is
 * created with an initial number of points (default: 64) and a center and
 * range which are taken from the data PixelMap. The initial number of grid
 * points can be adjusted with the `gridPoints()` method. Changing the data or
 * number of initial grid points will delete the current grid and create a new
 * instance with the given number of initial points, center and range.
 * 
 * Each source can have a maximum number of lensed images (default: 100).
 * This maximum number can be changed using the `imagesSize()` method.
 */
template<typename SourceType, typename LensType>
class ImageLikelihood
{
public:
	/**
	 * Indicate that this is a log-likelihood.
	 */
	static const bool is_logarithmic = true;
	
	/**
	 * Parameters for the source and lens.
	 */
	typedef SourceLensParameters<SourceType, LensType> parameter_type;
	
	/**
	 * \brief Construct ImageLikelihood for the given source and lens.
	 * 
	 * Use the provided source and lens to calculate the likelihood. The
	 * ImageLikelihood takes ownership of both objects (due to Model) and will
	 * delete them on destruction.
	 * 
	 * Before the likelihood can be calculated, the data PixelMap has to be set
	 * using the `data()` method.
	 * 
	 * \param source The source object.
	 * \param lens The lens object.
	 * \param data The PixelData to be compared to simulated images.
	 */
	ImageLikelihood(SourceType* source, LensType* lens, PixelData data)
	: source(source), lens(lens),
	  dta(data),
	  dof(0), dim(0),
	  images_size(100), images(0),
	  grid(0), grid_points(64)
	{
		// create ImageInfo array
		images = new ImageInfo[images_size];
		
		// create cosmology
		cosmo = new COSMOLOGY();
		
		// bind source, lens and cosmology together using model
		model = new Model(lens, source, cosmo);
		
		// set grid to match data
		grid_range = dta.getRange();
		std::copy(dta.getCenter(), dta.getCenter() + 2, grid_center);
		
		// create grid
		regrid();
		
		// calculate dof
		redof();
	}
	
	/**
	 * \brief Destructor.
	 * 
	 * Destruct the ImageLikelihood. The Model destructor will also destroy
	 * the source and lens.
	 */
	~ImageLikelihood()
	{
		delete[] images;
		
		delete grid;
		
		delete model; // deletes lens, source and cosmo
	}
	
	/** Get the number of free parameters. */
	unsigned long dimension() const { return dim; }
	
	/** Set the number of free parameters. */
	void dimension(unsigned long n) { dim = n; redof(); }
	
	/** Get the data PixelMap. */
	PixelData data() const { return dta; }
	
	/** Set the data PixelMap. */
	void data(PixelData data)
	{
		// copy data
		swap(dta, data);
		
		// change grid to match data
		grid_range = dta.getRange();
		std::copy(dta.getCenter(), dta.getCenter() + 2, grid_center);
		
		// recreate grid
		regrid();
		
		// recalculate degrees of freedom
		redof();
	}
	
	/** Get the maximum number of lensed images. */
	std::size_t imagesSize() const { return images_size; }
	
	/** Set the maximum number of lensed images. */
	void imagesSize(std::size_t size)
	{
		images_size = size;
		
		delete[] images;
		images = new ImageInfo[images_size];
	}
	
	/** Get the number of initial grid points. */
	std::size_t gridPoints() const { return grid_points; }
	/** Set the number of initial grid points. */
	void gridPoints(std::size_t n) { grid_points = n; regrid(); }
	
	/**
	 * Calculate the likelihood for the given source and lens parameters.
	 */
	double operator()(const parameter_type& params)
	{
		using std::log;
		
		// load parameters into source and lens
		params.source >> (*source);
		params.lens >> (*lens);
		
		// flush the surface brightness
		grid->RefreshSurfaceBrightnesses(source);
		
		// reinitialize the grid
		grid->ReInitializeGrid(lens);
		
		// number of images generated
		int image_count;
		
		// do the mapping
		map_images(
			lens, // Lens* lens
			source, // Source* source
			grid, // Grid* grid
			&image_count, // int* Nimages
			images, // ImageInfo* imageinfo
			(int)images_size, // int Nimagesmax
			source->getRadius(), // double xmax
			0.1*source->getRadius(), // double xmin
			0, // double initial_size
			EachImage, // ExitCriterion criterion
			true, // bool kappa_off
			false, // bool FindCenter
			true // bool divide_images
		);
		
		// build the fit image
		// TODO: loop for multi source
		PixelMap image(dta.getCenter(),dta.getNpixels(),dta.getResolution());
		image.AddImages(images, image_count, false);
		
		// calculate chi^2 for image and data
		double chi2 = dta.chi_square(image);
		
		// return chi^2 log-likelihood
		return (-chi2/2);//+ log(chi2)*(0.5*dof - 1);
	}
	
private:
	inline void redof()
	{
		// TODO: calculate degrees of freedom from dta and dim
		dof = 0;
	}
	
	inline void regrid()
	{
		delete grid;
		grid = new Grid(lens, (int)grid_points, grid_center, grid_range);
	}
	
	SourceType* source;
	LensType* lens;
	
	PixelData dta;
	
	unsigned long dof;
	unsigned long dim;
	
	std::size_t images_size;
	ImageInfo* images;
	
	COSMOLOGY* cosmo;
	Model* model;
	
	std::size_t grid_points;
	double grid_center[2];
	double grid_range;
	Grid* grid;
};

#endif
