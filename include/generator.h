#ifndef GENERATOR_H_
#define GENERATOR_H_

#include "parameters.h"

/**
 * \brief Utility class for creating generators.
 * 
 * This utility class can be subclassed by generators. It automatically adds
 * the necessary static information from the template arguments.
 */
template<typename SourceType, typename LensType>
class Generator
{
public:
	/**
	 * \brief The underlying source type.
	 * 
	 * Stores the source type the generator is made for. This is not the
	 * parameter type the generator is working on, see `source_parameters`
	 * below.
	 */
	typedef SourceType source_type;
	
	/**
	 * \brief The source parameter type.
	 * 
	 * This is the type of source parameters the generator is working on.
	 */
	typedef Parameters<SourceType> source_params;
	
	/**
	 * \brief The underlying lens type.
	 * 
	 * Stores the lens type the generator is made for. This is not the
	 * parameter type the generator is working on, see `lens_parameters`
	 * below.
	 */
	typedef LensType lens_type;
	
	/**
	 * \brief The lens parameter type.
	 * 
	 * This is the type of lens parameters the generator is working on.
	 */
	typedef Parameters<LensType> lens_params;
	
	/**
	 * \brief Type of combined source and lens parameters.
	 * 
	 * This is the full parameter type for both source and lens.
	 */
	typedef SourceLensParameters<SourceType, LensType> source_lens_params;
	
	/**
	 * The value type required by monaco.
	 */
	typedef source_lens_params value_type;
};

#endif
