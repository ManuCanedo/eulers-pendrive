#include <iostream>
#include <array>

#include "euler_math.h"

int main()
{
	const auto ode = [](double t, double x) {
		const double t2 = t * t;
		return t2 * x + t2 * std::sin(t2 * t);
	};
	const double ini_t = 0;
	const double ini_val = 1;
	const double fin_t = 2;

	const auto sym_sol = [](double t) {
		const double t3 = t * t * t;
		return -0.3 * std::cos(t3) - 0.1 * std::sin(t3) +
		       1.3 * std::exp(t3 / 3);
	};
	const double sol2 = sym_sol(2);
	printf("Real value: %0.5f\n\n", sol2);

	const std::array<double, 6> hs = {
		0.2, 0.1, 0.01, 1.0 / 1000, 1.0 / 10000, 1.0 / 100000
	};
	for (auto &h : hs) {
		const double eu = euler_ode_solve(ini_t, ini_val, fin_t, ode, h);
		printf("Euler's approximation for h=%0.5f:\n\tresult: %0.6f\n\terror: %0.20f\n",
		       h, eu, std::abs(eu - sol2));

		const double rk4 = runge_kutta_ode_solve(ini_t, ini_val, fin_t, ode, h);
		printf("Runge-Kutta's approximation for h=%0.5f:\n\tresult: %0.6f\n\terror: %0.20f\n\n",
		       h, rk4, std::abs(rk4 - sol2));
	}

	return 0;
}