// Heat equation:
// Ut = A*Uxx

#include <cmath>

const double PI = std::acos(-1);
enum Scheme {Explicit, Implicit};
enum Mode {Approx_solution, Exact_solution, Compare};

namespace parameters
{
	using namespace std;

	// Solver mode (approximate/exact/compare):
	Mode solver_mode = Approx_solution;

	// Scheme type (explicit/implicit):
	Scheme scheme = Implicit;

	// Measure times:
	vector<double> m_times {5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100}; // C++11 initialization
	//vector<double> m_times {1,2,3,4,5,6,8}; // C++11 initialization

	// Equation coefficient (must be positive):
	const double A = 0.0029;

	// Time step:
	const double dT = 0.05;
	// Length step:
	const double dX = 0.01;
	// Left border:
	const double X_l = 0;
	// Right border:
	const double X_r = 1;
	// Start time:
	const double T_start = 0;

	// Initial condition:
	// U(X,0) = initial(X)
	double initial(double X)
	{
		return cos(PI*X) + sin(PI*X);
	}

	// Border conditions:
	// a1*U(X_l,T) + b1*Ux(X_l,T) = border_l(T)
	// a2*U(X_r,T) + b2*Ux(X_r,T) = border_r(T)

	//const double a1 = 1.0/2.0;
	//const double b1 = 1.0/2.0/PI;
	//const double a2 = -1.0/2.0;
	//const double b2 = -1.0/2.0/PI;

	const double a1 = 0;
	const double b1 = 1.0/PI;
	const double a2 = 0;
	const double b2 = -1.0/PI;

	//const double a1 = 1;
	//const double b1 = 0;
	//const double a2 = -1;
	//const double b2 = 0;

	double border_l(double T)
	{
		return exp(-A*PI*PI*T);
	}
	/*{
		return initial(X_l);
	}*/
	double border_r(double T)
	{
		return exp(-A*PI*PI*T);
	}
	/*{
		return initial(X_r);
	}*/

	// Exact solution for test:
	double exact_solution (double X, double T)
	{
		return (cos(PI*X) + sin(PI*X))*exp(-A*PI*PI*T);
	}

	// Check functions:
	bool check_courant ()
	{
		return (A*dT/dX/dX <= 1.0/2.0);
	}
	bool check_A ()
	{
		return (A > 0);
	}
	bool check_borders ()
	{
		return (a1 != 0 || b1 != 0) && (a2 != 0 || b2 != 0);
	}
	bool check_length_borders ()
	{
		return ((X_r-X_l >= dX) && (X_l >= 0) && (X_r > 0));
	}
	bool check_time_borders ()
	{
		return (T_start >= 0);
	}
	bool check_left_border_cond ()
	{
		return ((a1 != 0) || (b1 != 0));
	}
	bool check_right_border_cond ()
	{
		return ((a2 != 0) || (b2 != 0));
	}
}

namespace P = parameters;
