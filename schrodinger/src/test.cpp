#include <iostream>
#include "solver.h"

// Schrodinger equation:
// 
// i*h*Ut = K(x,t)*Uxx + P(x,t,U)
// U(t=0) = initial(X)
// U(X=X_l) = border_l(T)
// U(X=X_r) = border_r(T)

using namespace std;

int main(int argc, char *argv[])
{
	/*cout << "Measure times:\n";
	for (int i = 0; i < P::m_times.size(); ++i)
	{
		cout << P::m_times[i];
		if (i != P::m_times.size()-1)
			cout << ' ';
		else
			cout << endl;
	}
	cout << endl;*/

	/*int time_size = 10000;
	P::m_times.resize(time_size);
	double T = 0;
	for (int i = 0; i < P::m_times.size(); ++i)
	{
		P::m_times[i] = T;
		T += P::dT;
	}*/
	P::init_time();

	Solver solver;
	solver.init();
	if (argc == 1)
		solver.solve(P::solver_mode, cout);
	else
	{
		ofstream file(argv[1]);
		solver.solve(P::solver_mode, file);
	}
	return 0;
}
