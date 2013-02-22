#ifndef IMAGE_LIKELIHOOD_H_
#define IMAGE_LIKELIHOOD_H_

#include "slsimlib.h"

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
template<typename Source, typename Lens>
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
	typedef SourceLensParameters<Source, Lens> parameter_type;
	
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
	 */
	ImageLikelihood(Source* source, Lens* lens)
	: source(source), lens(lens),
	  dof(0),
	  off(0), ns(0), norm(1),
	  images_size(100), images(0),
	  grid_points(64), grid(0)
	{
		// create ImageInfo array
		images = new ImageInfo[images_size];
		
		// create cosmology
		cosmo = new COSMOLOGY();
		
		// bind source, lens and cosmology together using model
		model = new Model(lens, source, cosmo);
		
		// these are then set in data()
		std::fill(grid_center, grid_center + 2, 0);
		grid_range = 0;
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
	
	/** Get the data PixelMap. */
	PixelMap data() const { return dta; }
	
	/** Set the data PixelMap. */
	void data(PixelMap data)
	{
		// copy pixel map
		swap(dta, data);
		
		// check if mask is compatible or create empty
		if(dta.size() != msk.base_size())
			msk = PixelMask(dta.size());
		
		// change grid to match data
		grid_range = dta.getRange();
		std::copy(dta.getCenter(), dta.getCenter() + 2, grid_center);
		
		// recreate grid
		regrid();
		
		// recalculate degrees of freedom
		redof();
	}
	
	/** Get the constant background offset. */
	double offset() const { return off; }
	/** Set the constant background offset. */
	void offset(double offset) { off = offset; }
	
	/** Get the constant background noise. */
	double noise() const { return ns; }
	/** Set the constant background noise. */
	void noise(double noise) { ns = noise; }
	
	/** Get the pixel count normalization. */
	double normalization() const { return norm; }
	/** Set the pixel count normalization. */
	void normalization(double n) { norm = n; }
	
	/** Get the mask. */
	PixelMask mask() const { return msk; }
	
	/** Set the mask. */
	bool mask(PixelMask mask)
	{
		// make sure size of data and mask PixelMap agree
		if(dta.valid() && mask.base_size() != dta.size())
			return false;
		
		// swap current mask and given one
		swap(msk, mask);
		
		// recalculate degrees of freedom
		redof();
		
		// success
		return true;
	}
	
	/** Set the mask using a PixelMap. */
	bool mask(PixelMap mask_map)
	{
		return mask(PixelMask(mask_map));
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
		
		// make sure there is data
		if(!dta.valid())
			return -1;
		
		// load parameters into source and lens
		set(*source, params.source);
		set(*lens, params.lens);
		
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
		PixelMap image(dta.getNpixels(), dta.getRange(), dta.getCenter());
		image.AddImages(images, image_count, false);
		
		// calculate chi^2 for image and data
		double chi2 = chi_square(dta, image, off, ns, norm, msk);
		
		// return chi^2 log-likelihood
		return (-chi2/2) + log(chi2)*(0.5*dof - 1);
	}
	
private:
	inline void redof()
	{
		// calculate degrees of freedom
		// TODO: take number of parameters into account
		dof = static_cast<unsigned long>(msk.valid() ? dta.size() : msk.size());
	}
	
	inline void regrid()
	{
		delete grid;
		grid = new Grid(lens, (int)grid_points, grid_center, grid_range);
	}
	
	Source* source;
	Lens* lens;
	
	PixelMap dta;
	unsigned long int dof;
	
	double off;
	double ns;
	double norm;
	
	PixelMask msk;
	
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
