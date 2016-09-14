// Schrodinger equation:
// Ut = K*Uxx + P*U

#include <cmath>
#include <complex>

typedef std::complex<double> complexd;

const complexd imag_unit(0,1);
const double Plank_const = 1.0;
const double PI = std::acos(-1);

double gaussian (double X, double m, double s)
{
	return exp(-(X-m)*(X-m) / (2*s*s))/(s * sqrt(2*PI));	
}

enum Scheme {Explicit, Implicit};
enum Solver_mode {Approx_solution, Exact_solution, Compare};

namespace parameters
{
	using namespace std;

	// Solver mode (approximate/exact/compare):
	Solver_mode solver_mode = Approx_solution;

	// Scheme type (explicit/implicit):
	Scheme scheme = Implicit;

	// Measure times:
	//vector<double> times {5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100}; // C++11 initialization
	vector<double> times; // Initialization is below

	// Laplas operator coefficient (kinetic energy)
	complexd K (double X, double T)
	{
		return -1;
	}

	// Potential energy operator coefficient
	complexd P (double X, double T)
	{
		/*const double small = 0.0001;
		const double big = 1000.0;
		if (X <= 1 && X >= -1)
			return small;
		else
			return big;*/
		return X*X;
	}

	// Time step:
	const double dT = 0.001;
	// Length step:
	const double dX = 0.001;
	// Left border:
	const double X_l = -3;
	// Right border:
	const double X_r = 3;
	// Start time:
	const double T_start = 0;

	// Initial condition:
	// U(X,0) = initial(X)
	complexd initial(double X)
	{
		/*double sigma = 0.1;
		double mu = 0;
		return gaussian(X,mu,sigma);*/
		return cos(PI*X/6);
	}

	// Border conditions:
	// a1*U(X_l,T) + b1*Ux(X_l,T) = border_l(T)
	// a2*U(X_r,T) + b2*Ux(X_r,T) = border_r(T)

	const complexd a1 = 1;
	const complexd b1 = 0;
	const complexd a2 = 1;
	const complexd b2 = 0;

	complexd border_l(double T)
	{
		return 0;
	}
	/*{
		return initial(X_l);
	}*/
	complexd border_r(double T)
	{
		return 0;
	}
	/*{
		return initial(X_r);
	}*/

	// Exact solution for test:
	complexd exact_solution (double X, double T)
	{
		return 0;
	}

	// Check functions:
	bool check_borders ()
	{
		return (a1 != 0.0 || b1 != 0.0) && (a2 != 0.0 || b2 != 0.0);
	}
	bool check_length_borders ()
	{
		return (X_r-X_l >= dX);
	}
	bool check_time_borders ()
	{
		return (T_start >= 0);
	}
	bool check_left_border_cond ()
	{
		return ((a1 != 0.0) || (b1 != 0.0));
	}
	bool check_right_border_cond ()
	{
		return ((a2 != 0.0) || (b2 != 0.0));
	}
	void init_time()
	{
		int time_vec_size = 10000;
		times.resize(time_vec_size);
		double T = 0, T_step = dT;
		for (int i = 0; i < times.size(); ++i)
		{
			times[i] = T;
			T += T_step;
		}
	}
}

namespace P = parameters;
