#include "../include/InputParams.h"

// array of all parameters
namespace {
	const char* parameter_list[] = {
		// General
		"outputfile",
		"deflection_off",
		"read_redshift_planes",
		"redshift_planes_fil",
		"z_lens",
		"z_source",
		
		// Cosmology
		"Omega_matter",
		"Omega_lambda",
		"Omega_baryon",
		"Omega_neutrino",
		"hubble",
		"sigm_8",

		// Main halos type
		"main_halo_on",
		"main_DM_halo_type",
		"main_galaxy_halo_type",

		// Field halos
		"field_off",
		"field_Nplanes",

		// Field halos type
		"field_internal_profile",
		"field_internal_profile_galaxy",
		"field_prof_internal_slope_pnfw",
		"field_prof_internal_slope_pl",

		// Field halos from a mass function
		"field_mass_func_type",
		"field_fov",
		"field_buffer",
		"field_min_mass",
		
		// Field halos from an input file
		"field_input_simulation_file",
		
		// AnaNSIE lens halo model
		"main_sigma",
		"main_core",
		"main_axis_ratio",
		"main_pos_angle",
		"main_NDistortionModes",
		"main_perturb_beta",
		"main_perturb_kappa",
		"main_perturb_gamma",
		"main_perturb_monopole",
		"main_perturb_quadrapole",
		"main_perturb_hexopole",
		"main_perturb_octopole",
		"main_sub_Ndensity",
		"main_sub_beta",
		"main_sub_alpha",
		"main_sub_Rmax",
		"main_sub_mass_max",
		"main_sub_mass_min",
		"main_sub_type",
		
		// MOKA lens halo model
		"MOKA_input_file",
		"MOKA_input_params",
		"MOKA_analyze",
		"background_field",

		// NFW lens halo model
		"mass_nfw",
		"Rmax_nfw",
		"zlens_nfw",
		"concentration_nfw",

		// Pseudo NFW lens halo model
		"mass_pnfw",
		"Rmax_pnfw",
		"zlens_pnfw",
		"concentration_pnfw",
		"slope_pnfw",

		// Power law lens halo model
		"mass_pl",
		"Rmax_pl",
		"zlens_pl",
		"slope_pl",

		// Simple NSIE lens halo model
		"mass_nsie",

		// Hernquist lens halo model
		"mass_hernquist",
		"Rmax_hernquist",
		"zlens_hernquist",
		"rscale_hernquist",
		
		// Jaffe lens halo model
		"mass_jaffe",
		"Rmax_jaffe",
		"zlens_jaffe",
		"rscale_jaffe",

		// Uniform lens halo model
		"zlens_uniform",
		"kappa_uniform",
		"gamma_uniform_1",
		"gamma_uniform_2",
		"zsource_reference",

		// Stars
		"main_stars_N",
		"main_stars_fraction",
		"main_stars_mass",
		"main_stars_mass_function",
		"main_stars_min_mass",
		"main_stars_max_mass",
		"main_stars_bending_point",
		"main_stars_lo_mass_slope",
		"main_stars_hi_mass_slope",

		// MultiDark lenses
		"MultiDark_input_file",
		
		// Type of source SB model
		"SourceSBType",
		
		// Gaussian source model
		"gauss_r2",
		
		// BLR source model
		"BHmass",
		"gamma",
		"inclin",
		"opening_ang",
		"r_in",
		"r_out",
		"nuo",
		"source_sigma"
	};
}

// initialize the static set of all known labels
const std::set<std::string> InputParams::known_labels(parameter_list, parameter_list + sizeof(parameter_list)/sizeof(parameter_list[0]));
