#include <cmath>
#include <vector>
#include <functional>

double runge_kutta_ode_solve(double ini_t, double ini_val, double fin_t,
			     const std::function<double(double, double)> &ode,
			     double h)
{
	std::vector<double> iters((fin_t - ini_t) / h + 1);

	for (std::size_t k = 0; k < iters.size(); ++k) {
		if (!k) {
			iters[k] = ini_val;
			continue;
		}
		const double t_k = ini_t + h * k;

		const double f_k1 = ode(t_k, iters[k - 1]);
		const double f_k2 =
			ode(t_k + h / 2, iters[k - 1] + f_k1 * h / 2);
		const double f_k3 =
			ode(t_k + h / 2, iters[k - 1] + f_k2 * h / 2);
		const double f_k4 = ode(t_k + h, iters[k - 1] + f_k3 * h);

		iters[k] = iters[k - 1] +
			   h / 6 * (f_k1 + 2 * f_k2 + 2 * f_k3 + f_k4);
	}
	return iters.back();
}

double euler_ode_solve(double ini_t, double ini_val, double fin_t,
		       const std::function<double(double, double)> &ode,
		       double h)
{
	std::vector<double> iters((fin_t - ini_t) / h + 1);

	for (std::size_t k = 0; k < iters.size(); ++k) {
		if (!k) {
			iters[k] = ini_val;
			continue;
		}
		const double t_k = ini_t + h * k;
		const double f_k = ode(t_k, iters[k - 1]);
		iters[k] = iters[k - 1] + h * f_k;
	}
	return iters.back();
}

/*
*   Calculate the value an Ordinary Differential Equation 'ode' will return in a time 'fin_t'
*   based on an initial time 'ini_t' and value 'ini_val'.
*/
double ode_solve(double ini_t, double ini_val, double fin_t,
		 const std::function<double(double, double)> &ode,
		 double h = (1.0 / 1000000))
{
	/*
    *   Using the Runge-Kutta's method for high precission.
    *   If not, using the simpler Euler's method
    */
	if (h < 1.0 / 100000)
		return runge_kutta_ode_solve(ini_t, ini_val, fin_t, ode, h);
	return euler_ode_solve(ini_t, ini_val, fin_t, ode, h);
}