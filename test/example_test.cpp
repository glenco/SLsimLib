// must be included
#include <gtest/gtest.h>

#include "quadrature.h"

namespace {

// A two-point closed rule is sufficient to demonstrate the quadrature API
// with a constant integrand.
struct TwoPointRule {
	static const int n = 2;
	static const double p[n];
	static const double w[n];
	static const double e[n];
};

const double TwoPointRule::p[TwoPointRule::n] = {-1.0, 1.0};
const double TwoPointRule::w[TwoPointRule::n] = {1.0, 1.0};
const double TwoPointRule::e[TwoPointRule::n] = {0.0, 0.0};


TEST(QuadratureTest, IntegratesAConstantOverARectangle) {
	const ISOP::quadrature_result result = ISOP::quadrature<TwoPointRule>(
		[](double, double) { return 3.0; },
		-1.0, 2.0,
		0.0, 4.0,
		1.0e-12, 1.0e-12
	);

	EXPECT_DOUBLE_EQ(result.value, 36.0);
	EXPECT_DOUBLE_EQ(result.error, 0.0);
}

}  // namespace
