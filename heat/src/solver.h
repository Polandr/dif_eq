#include <iostream>
#include <fstream>
#include <vector>

#include "parameters.hpp"

class Solver
{
	Scheme s;
	double T;
	double T_fin;
	std::vector<double> U;
	std::vector<double> X;
	std::vector<double> evolution_step ();

public:

	void check_params () const;
	void init ();
	void solve (Mode mode, std::ostream& out);
	void export_header (std::ostream& out) const;
	void export_values (std::ostream& out) const;
};

#include "solver.hpp"
