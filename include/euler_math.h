
#include <cmath>
#include <vector>
#include <functional>

double runge_kutta_ode_solve(double ini_t, double ini_val, double fin_t,
			     const std::function<double(double, double)> &ode,
			     double h)
{
	double val = ini_val;
	for (double t = 1; t < fin_t + h; t += h) {
		const double f_k1 = ode(t, val);
		const double f_k2 = ode(t + h / 2, val + f_k1 * h / 2);
		const double f_k3 = ode(t + h / 2, val + f_k2 * h / 2);
		const double f_k4 = ode(t + h, val + f_k3 * h);
		val += h * (f_k1 + 2 * f_k2 + 2 * f_k3 + f_k4) / 6;
	}
	return val;
}

double euler_ode_solve(double ini_t, double ini_val, double fin_t,
		       const std::function<double(double, double)> &ode,
		       double h)
{
	double val = ini_val;
	for (double t = 1; t < fin_t + h; t += h)
		val += h * ode(t, val);

	return val;
}

/*
*   Calculate the value an Ordinary Differential Equation 'ode' will return in a time 'fin_t'
*   based on an initial time 'ini_t' and value 'ini_val'.
*/
double ode_solve(double ini_t, double ini_val, double fin_t,
		 const std::function<double(double, double)> &ode,
		 double h = 0.1)
{
	/*
    *   Using the Runge-Kutta's method for high precission.
    *   If not, using the simpler Euler's method
    */
	if (h > 0.1)
		return euler_ode_solve(ini_t, ini_val, fin_t, ode, h);
	return runge_kutta_ode_solve(ini_t, ini_val, fin_t, ode, h);
}

