#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "lens_halos.h"
#include "particle_halo.h"
#include "particle_halo2.h"
#include "slsimlib.h"

using HaloType = LensHaloNFW;

constexpr double ARCSEC_TO_RAD = PI / 180.0 / 60.0 / 60.0;
constexpr double RELATIVE_TOLERANCE = 1.0e-6;

struct HaloSpec {
	double mass;
	double x_arcsec;
	double y_arcsec;
};

struct ForceResult {
	double alpha[2] = {0.0, 0.0};
	KappaType kappa = 0.0;
	KappaType gamma[3] = {0.0, 0.0, 0.0};
	KappaType phi = 0.0;
};

struct ComparisonResult {
	ForceResult exact;
	ForceResult lens_h;
	ForceResult lens_halos;
};

template <class Halo>
ForceResult evaluate(Halo &halo, double x_phys, double y_phys) {
	ForceResult result;
	double position[2] = {x_phys, y_phys};
	halo.force_halo(
		result.alpha,
		&result.kappa,
		result.gamma,
		&result.phi,
		position
	);
	return result;
}

void add_result(ForceResult &sum, const ForceResult &result) {
	sum.alpha[0] += result.alpha[0];
	sum.alpha[1] += result.alpha[1];
	sum.kappa += result.kappa;
	sum.gamma[0] += result.gamma[0];
	sum.gamma[1] += result.gamma[1];
	sum.gamma[2] += result.gamma[2];
	sum.phi += result.phi;
}

double vector_norm(double x, double y) {
	return std::sqrt(x * x + y * y);
}

double relative_vector_error(
	const double actual[2],
	const double expected[2]
) {
	return vector_norm(
		actual[0] - expected[0],
		actual[1] - expected[1]
	) / std::max(vector_norm(expected[0], expected[1]), 1.0e-30);
}

double relative_scalar_error(double actual, double expected) {
	return std::abs(actual - expected) / std::max(std::abs(expected), 1.0);
}

void expect_force_near(
	const ForceResult &actual,
	const ForceResult &expected,
	double tolerance
) {
	EXPECT_LE(relative_vector_error(actual.alpha, expected.alpha), tolerance)
		<< "alpha differs from the direct halo sum";
	EXPECT_LE(relative_vector_error(actual.gamma, expected.gamma), tolerance)
		<< "gamma differs from the direct halo sum";
	EXPECT_LE(relative_scalar_error(actual.kappa, expected.kappa), tolerance)
		<< "kappa differs from the direct halo sum";
}

std::vector<HaloSpec> make_halo_specs(double scale) {
	std::vector<HaloSpec> specs;
	specs.reserve(49);

	int index = 0;
	for(int iy = -3; iy <= 3; ++iy) {
		for(int ix = -3; ix <= 3; ++ix) {
			HaloSpec spec;
			spec.mass = 3.0e11 * (1.0 + ((index * 11) % 30));
			spec.x_arcsec = (ix * 35.0 + 3.0 * iy + 11.0) * scale;
			spec.y_arcsec = (iy * 30.0 - 2.0 * ix - 9.0) * scale;
			specs.push_back(spec);
			++index;
		}
	}

	return specs;
}

std::vector<HaloType> build_halos(
	const std::vector<HaloSpec> &specs,
	std::size_t count,
	double lens_redshift,
	COSMOLOGY &cosmology
) {
	std::vector<HaloType> halos;
	halos.reserve(count);

	for(std::size_t i = 0; i < count; ++i) {
		const HaloSpec &spec = specs[i];
		HALOCalculator calculator(&cosmology, spec.mass, lens_redshift);

		halos.emplace_back(
			spec.mass,
			calculator.getRvir(),
			lens_redshift,
			calculator.getConcentration(),
			1.0,
			0.0,
			cosmology
		);
		halos.back().setTheta(
			spec.x_arcsec * ARCSEC_TO_RAD,
			spec.y_arcsec * ARCSEC_TO_RAD
		);
	}

	return halos;
}

ForceResult exact_sum(
	std::vector<HaloType> &halos,
	const std::vector<HaloSpec> &specs,
	std::size_t count,
	double lens_distance
) {
	ForceResult sum;

	for(std::size_t i = 0; i < count; ++i) {
		const double x = -specs[i].x_arcsec * ARCSEC_TO_RAD * lens_distance;
		const double y = -specs[i].y_arcsec * ARCSEC_TO_RAD * lens_distance;
		add_result(sum, evaluate(halos[i], x, y));
	}

	return sum;
}

ComparisonResult compare_configuration(
	double scale,
	std::size_t count,
	double theta_open,
	double lens_redshift,
	COSMOLOGY &cosmology
) {
	const std::vector<HaloSpec> specs = make_halo_specs(scale);
	std::vector<HaloType> exact_halos =
		build_halos(specs, count, lens_redshift, cosmology);
	std::vector<HaloType> halos_for_h(exact_halos);
	std::vector<HaloType> halos_for_halos(exact_halos);

	const double inverse_area = 0.0;
	const double hole_radius = 5.0 / 60.0 * PI / 180.0;
	const int bucket_size = 5;

	LensHaloH<HaloType> lens_h(
		halos_for_h,
		lens_redshift,
		inverse_area,
		hole_radius,
		cosmology,
		bucket_size,
		theta_open
	);
	LensHaloHalos<HaloType> lens_halos(
		halos_for_halos,
		lens_redshift,
		cosmology,
		inverse_area,
		theta_open
	);

	ComparisonResult result;
	result.exact = exact_sum(
		exact_halos,
		specs,
		count,
		cosmology.angDist(lens_redshift)
	);
	result.lens_h = evaluate(lens_h, 0.0, 0.0);
	result.lens_halos = evaluate(lens_halos, 0.0, 0.0);
	return result;
}

COSMOLOGY &test_cosmology() {
	static COSMOLOGY cosmology(CosmoParamSet::Planck18);
	return cosmology;
}

TEST(LensHaloHalosTest, SingleHaloWrappersMatchBareHalo) {
	COSMOLOGY &cosmology = test_cosmology();
	const double lens_redshift = 0.5;
	const double mass = 1.0e13;
	HALOCalculator calculator(&cosmology, mass, lens_redshift);

	std::vector<HaloType> halos;
	halos.emplace_back(
		mass,
		calculator.getRvir(),
		lens_redshift,
		calculator.getConcentration(),
		1.0,
		0.0,
		cosmology
	);
	halos.front().setTheta(0.0, 0.0);

	HaloType bare = halos.front();
	std::vector<HaloType> halos_for_h(halos);
	std::vector<HaloType> halos_for_halos(halos);
	LensHaloH<HaloType> lens_h(
		halos_for_h,
		lens_redshift,
		0.0,
		5.0 / 60.0 * PI / 180.0,
		cosmology,
		5,
		0.05
	);
	LensHaloHalos<HaloType> lens_halos(
		halos_for_halos,
		lens_redshift,
		cosmology,
		0.0,
		0.05
	);

	const double radius =
		2.0 * ARCSEC_TO_RAD * cosmology.angDist(lens_redshift);
	const ForceResult expected = evaluate(bare, radius, 0.0);
	expect_force_near(evaluate(lens_h, radius, 0.0), expected, 1.0e-12);
	expect_force_near(evaluate(lens_halos, radius, 0.0), expected, 1.0e-12);
}

TEST(LensHaloHalosTest, LensHaloHMatchesDirectSum) {
	COSMOLOGY &cosmology = test_cosmology();
	const double lens_redshift = 0.5;
	const std::size_t counts[] = {1, 2, 3, 4, 5, 6, 8, 12, 20, 30, 49};

	for(std::size_t count : counts) {
		SCOPED_TRACE(testing::Message() << "halo count: " << count);
		const ComparisonResult result = compare_configuration(
			1.0, count, 0.0, lens_redshift, cosmology
		);
		expect_force_near(result.lens_h, result.exact, RELATIVE_TOLERANCE);
	}
}

TEST(LensHaloHalosTest, MatchesDirectSumForSeparatedHalos) {
	ComparisonResult result = compare_configuration(
		1.25, 49, 0.0, 0.5, test_cosmology()
	);
	expect_force_near(result.lens_halos, result.exact, RELATIVE_TOLERANCE);
}

TEST(LensHaloHalosTest, MatchesDirectSumForCompactHalos) {
	ComparisonResult result = compare_configuration(
		1.0, 49, 0.05, 0.5, test_cosmology()
	);
	expect_force_near(result.lens_halos, result.exact, RELATIVE_TOLERANCE);
}

TEST(LensHaloHalosTest, MatchesDirectSumAcrossTreeThreshold) {
	COSMOLOGY &cosmology = test_cosmology();
	const std::size_t counts[] = {1, 2, 3, 4, 5, 6, 8, 12, 20, 30, 49};

	for(std::size_t count : counts) {
		SCOPED_TRACE(testing::Message() << "halo count: " << count);
		const ComparisonResult result = compare_configuration(
			1.0, count, 0.0, 0.5, cosmology
		);
		expect_force_near(
			result.lens_halos,
			result.exact,
			RELATIVE_TOLERANCE
		);
	}
}
