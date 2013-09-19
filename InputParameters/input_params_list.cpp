#include "../include/InputParams.h"

// array of all parameters
namespace {
	const char* parameter_list[] = {
		// General
		"outputfile",
		"deflection_off",
    "lensing_off",
		"read_redshift_planes",
		"redshift_planes_file",
		"z_lens",
		"z_source",
		
		// Cosmology
		"Omega_matter",
		"Omega_lambda",
		"Omega_baryon",
		"Omega_neutrino",
		"hubble",
		"sigma_8",

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
		"field_mass_func_alpha",
		"field_fov",
		"field_buffer",
		"field_min_mass",
		
		// Field halos from an input file
		"field_input_simulation_file",
		
		// Main halos
		"main_mass",
		"main_zlens",
		"main_Rmax",
		"main_concentration",
		"main_slope",
		"main_sigma",
		"main_core",
		"main_axis_ratio",
		"main_pos_angle",
		"main_rscale",

		// AnaNSIE perturbations
		"main_NDistortionModes",
		"main_perturb_beta",
		"main_perturb_kappa",
		"main_perturb_gamma",
		"main_perturb_monopole",
		"main_perturb_quadrapole",
		"main_perturb_hexopole",
		"main_perturb_octopole",

		// AnaNSIE substructures
		"main_sub_Ndensity",
		"main_sub_beta",
		"main_sub_alpha",
		"main_sub_Rmax",
		"main_sub_mass_max",
		"main_sub_mass_min",
		"main_sub_type",
		
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

		// MOKA lens halo model
		"MOKA_input_file",
		"MOKA_input_params",
		"MOKA_analyze",
		"MOKA_background_field",

		// Uniform lens halo model
		"kappa_uniform",
		"gamma_uniform_1",
		"gamma_uniform_2",

		// Reference redshift for halo quantities that depend on source z
		"zsource_reference",

		// MultiDark lenses
		"MultiDark_input_file",
		
		// Type of source SB model
		"SourceSBType",
		
		// Gaussian source model
		"gauss_r2",
		
		// BLR source model
		"source_BHmass",
		"source_gamma",
		"source_inclin",
		"source_opening_ang",
		"source_r_in",
		"source_r_out",
		"source_nuo",
		"souce_fK",

		// SourceMultiAnaGalaxy & Sky
		"source_input_galaxy_file",
		"source_band",
		"source_mag_limit",
		"source_sb_limit",
		"input_sersic_file"
	};
}

// initialize the static set of all known labels
const std::set<std::string> InputParams::known_labels(parameter_list, parameter_list + sizeof(parameter_list)/sizeof(parameter_list[0]));
