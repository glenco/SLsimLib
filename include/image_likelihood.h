#ifndef IMAGE_LIKELIHOOD_H_
#define IMAGE_LIKELIHOOD_H_

#include "slsimlib.h"

template<typename Source, typename Lens>
class ImageLikelihood
{
public:
	typedef SourceLensParameters<Source, Lens> parameter_type;
	
	ImageLikelihood(Source* source, Lens* lens)
	: source(source), lens(lens),
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
		
		// these are actually set in data()
		pixels = 0;
		rng = 0;
		std::fill(ctr, ctr + 2, 0);
	}
	
	~ImageLikelihood()
	{
		delete[] images;
		
		delete model; // deletes lens, source and cosmo
		
		delete grid;
	}
	
	PixelMap data() const { return dta; }
	
	void data(PixelMap data)
	{
		// copy pixmap
		using std::swap;
		swap(dta, data);
		
		// get properties
		pixels = dta.getNpixels();
		rng = dta.getRange();
		std::copy(dta.getCenter(), dta.getCenter() + 2, ctr);
		
		// recreate grid
		regrid();
	}
	
	double offset() const { return off; }
	void offset(double offset) { off = offset; }

	double noise() const { return ns; }
	void noise(double noise) { ns = noise; }
	
	double normalization() const { return norm; }
	void normalization(double n) { norm = n; }
	
	PixelMap mask() const { return msk; }
	void mask(PixelMap mask) { using std::swap; swap(msk, mask); }
	
	double range() const { return rng; }
	void range(double r) { rng = r; regrid(); }
	
	const double* center() const { return ctr; }
	void center(const double* c) { std::copy(c, c + 2, ctr); regrid(); }
	
	void dimensions(const double* center, double range)
	{
		std::copy(center, center + 2, ctr);
		rng = range;
		regrid();
	}
	
	std::size_t imagesSize() const
	{
		return images_size;
	}
	
	void imageSize(std::size_t size)
	{
		images_size = size;
		
		delete[] images;
		images = new ImageInfo[images_size];
	}
	
	std::size_t gridPoints() const
	{
		return grid_points;
	}
	
	void gridPoints(std::size_t n)
	{
		grid_points = n;
		regrid();
	}
	
	double operator()(const parameter_type& params)
	{
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
			images_size, // int Nimagesmax
			source->getRadius(), // double xmax
			0.1*source->getRadius(), // double xmin
			0, // double initial_size
			EachImage, // ExitCriterion criterion
			true, // bool kappa_off
			false, // bool FindCenter
			true // bool divide_images
		);
		
		// build the fit image
		PixelMap image(pixels, rng, ctr);
		image.AddImages(images, image_count, false);
		
		return chi_square(dta, image, off, ns, norm, msk);
	}
	
private:
	inline void regrid()
	{
		delete grid;
		grid = new Grid(lens, grid_points, ctr, rng);
	}
	
	Source* source;
	Lens* lens;
	
	PixelMap dta;
	
	double off;
	double ns;
	double norm;
	
	PixelMap msk;
	
	std::size_t pixels;
	double rng;
	double ctr[2];
	
	std::size_t images_size;
	ImageInfo* images;
	
	COSMOLOGY* cosmo;
	Model* model;
	
	std::size_t grid_points;
	Grid* grid;
};

#endif
