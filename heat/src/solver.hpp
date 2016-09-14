#include <cstring>
#include <cstdio>
#include <exception>

using namespace std;

// Solver exception class:

class Solver_exception: public std::exception
{
	mutable char* errstr; 

	public:

	Solver_exception(const char* str = "")
	{
		errstr = const_cast <char*> (str);
	}
	~Solver_exception() throw()
	{
		delete [] errstr;
	}
	virtual const char* what() const throw()
	{
		char* tmp = errstr;
		char* prefix = const_cast <char*> ("Solver class error: ");
		try
		{
			errstr = new char [strlen(prefix)+strlen(errstr)+2];
		}
		catch (std::exception& smth)
		{
			return "Couldn't generate an error message (there is no memory)\n";
		}
		sprintf(errstr, "%s%s.\n", prefix, tmp);
		return errstr;
	}
};

// Solver for triadiagonal matrix using shuttle method:
vector<double> shuttle_method (vector<double>& a, vector<double>& c, vector<double>& b, vector<double>& f)
// n - dimension
// a - lower diagonal
// b - upper diagonal
// c - central diagonal
// f - right-hand side
{
	int n = f.size();
	vector<double> x(n);
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i]/c[i-1];
		c[i] = c[i] - m*b[i-1];
		f[i] = f[i] - m*f[i-1];
	}

	x[n-1] = f[n-1]/c[n-1];

	for (int i = n - 2; i >= 0; i--)
		x[i] = (f[i] - b[i]*x[i+1]) / c[i];

	return x;
}

// Function which searches for the nearest value to value stated as 2nd argument:
double nearest_value (const vector<double>& values, const double value)
{
	double nearest = values[0];
	double dst = abs(value - values[0]);
	for (int i = 1; i < values.size(); ++i)
	{
		if (abs(value - values[i]) < dst)
		{
			dst = abs(value - values[i]);
			nearest = values[i];
		}
		if (abs(value - values[i]) == dst && values[i] > nearest)
			nearest = values[i];
	}
	return nearest;
}

// Max element from vector:
double max (vector<double>& v)
{
	double max = v[0];
	for (int i = 1; i < v.size(); ++i)
		if (v[i] > max)
			max = v[i];
	return max;
}

// Measure of proximity of solution:
double proximity_measure (const vector<double>& X, const double T, const vector<double>& U)
// max(di), di=|exact-approx|
/*{
	double max_dst = 0;
	for (int i = 0; i < X.size(); ++i)
	{
		double U_exact = P::exact_solution(X[i], T);
		if (abs((U_exact-U[i])/U_exact) > max_dst)
			max_dst = abs((U_exact-U[i])/U_exact);
	}
	return max_dst;
}*/
// max(ei), ei=|exact-approx|/|exact|
{
	double max_dst = 0;
	for (int i = 0; i < X.size(); ++i)
	{
		double U_exact = P::exact_solution(X[i], T);
		if (abs((U_exact-U[i])/U_exact) > max_dst)
			max_dst = abs((U_exact-U[i])/U_exact);
	}
	return max_dst;
}

//---------------------------------------------------------------------------------------------------------

void Solver::check_params () const
{
	if (!P::check_A())
		throw Solver_exception("incorrect equation coefficient (A)");
	if (!P::check_borders())
		throw Solver_exception("incorrect border conditions");
	if (!P::check_length_borders())
		throw Solver_exception("incorrect length borders");
	if (!P::check_time_borders())
		throw Solver_exception("incorrect time borders");
	if (!P::check_left_border_cond())
		throw Solver_exception("incorrect left border conditions");
	if (!P::check_right_border_cond())
		throw Solver_exception("incorrect right border conditions");
	if (s == Explicit && !P::check_courant())
		cerr << "WARNING!\nCourant condition is not observed (A*dT/dX/dX > 1/2)\n";
}

void Solver::init ()
{
	s = P::scheme;
	T = 0;
	T_fin = max(P::m_times);
	check_params();
	for (double Xi = P::X_l; Xi <= P::X_r; Xi += P::dX)
	{
		X.push_back(Xi);
		U.push_back(P::initial(Xi));
	}
}

vector<double> Solver::evolution_step ()
{
	// Explicit scheme:
	if (s == Explicit)
	{
		vector<double> tmp_U(U.size());

		for (int i = 1; i < U.size() - 1; ++i)
			tmp_U[i] = (U[i+1] - 2*U[i] + U[i-1]) * P::A * P::dT / P::dX / P::dX + U[i];
		tmp_U[0] = (P::border_l(T) - P::b1/P::dX*tmp_U[1]) / (P::a1 - P::b1/P::dX);
		tmp_U[U.size()-1] = (P::border_r(T) + P::b2/P::dX*tmp_U[U.size()-2]) / (P::a2 + P::b2/P::dX);
		return tmp_U;
	}

	// Implicit scheme:
	if (s == Implicit)
	{
		// Shuttle method:
		vector<double> upper(U.size()), centr(U.size()), lower(U.size()), right(U.size());
		for (int i = 0; i < U.size(); ++i)
		{
			if (i == 0)
			{
				lower[i] = 0;
				centr[i] = P::a1 - P::b1/P::dX;
				upper[i] = P::b1/P::dX;
				right[i] = P::border_l(T);
				continue;
			}
			if (i == U.size()-1)
			{
				lower[i] = -P::b2/P::dX;
				centr[i] = P::a2 + P::b2/P::dX;
				upper[i] = 0;
				right[i] = P::border_r(T);
				continue;
			}
			double c = P::A * P::dT / P::dX / P::dX;
			lower[i] = -c;
			centr[i] = 1 + 2*c;
			upper[i] = -c;
			right[i] = U[i];
		}

		return shuttle_method(lower, centr, upper, right);
	}
}

void Solver::solve (Mode mode, ostream& out)
{
	export_header(out);

	export_values(out);
	for (T = P::T_start + P::dT; T <= T_fin; T += P::dT)
	{
		if (mode == Approx_solution)
		// Standard evolution
		{
			U = evolution_step();
			if (abs(nearest_value(P::m_times, T) - T) <= P::dT/2)
				export_values(out);
		}
		if (mode == Compare)
		// Compre approximate solution with exact one
		{
			U = evolution_step();
			if (abs(nearest_value(P::m_times, T) - T) <= P::dT/2)
				out << proximity_measure(X,T,U) << endl;
		}
		if (mode == Exact_solution)
		// Evolution of exact solution
		{
			for (int i = 0; i < U.size(); ++i)
				U[i] = P::exact_solution(X[i],T);
			if (abs(nearest_value(P::m_times, T) - T) <= P::dT/2)
				export_values(out);
		}
	}
}

void Solver::export_header (ostream& out) const
{
	/*
	out << "Length sampling (" << X.size() << " elements):\n";
	for (int i = 0; i < X.size(); ++i)
	{
		out << X[i];
		if (i != X.size()-1)
			out << ' ';
		else
			out << endl;
	}
	out << "Measure times (" << P::m_times.size() << " elements):\n";
	for (int i = 0; i < P::m_times.size(); ++i)
	{
		if (P::m_times[i] - P::T_start >= 0)
		{
			double m_time = round((P::m_times[i] - P::T_start) / P::dT) * P::dT;
			out << m_time;
			if (i != P::m_times.size()-1)
				out << ' ';
			else
				out << endl;
		}
	}
	out << endl;
	*/
	
	for (int i = 0; i < X.size(); ++i)
	{
		out << X[i];
		if (i != X.size()-1)
			out << ' ';
		else
			out << endl;
	}
}

void Solver::export_values (ostream& out) const
{
	out << T << endl;
	
	for (int i = 0; i < U.size(); ++i)
	{
		out << U[i];
		if (i != U.size()-1)
			out << ' ';
		else
			out << endl;
	}	
}
