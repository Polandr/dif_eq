#include <iostream>
#include "solver.h"

// Equation:
// 
// Ut = A*Uxx
// U (t=0) = cos(PI*x) + sin(PI*x)
// 1/2*U + 1/(2*PI)*Ux (x=0) = exp(-A*PI*PI*t)
// -1/2*U - 1/(2*PI)*Ux (x=1) = exp(-A*PI*PI*t)
// 0 <= x <= 1
// t > 0
//
// Exact solution:
// U = (cos(PI*x) + sin(PI*x))*exp(-A*PI*PI*t)

using namespace std;

int main(int argc, char *argv[])
{
	/*cout << "A*dT/dX/dX = ";
	cout << P::A*P::dT/P::dX/P::dX << endl;
	cout << "Measure times:\n";
	for (int i = 0; i < P::m_times.size(); ++i)
	{
		cout << P::m_times[i];
		if (i != P::m_times.size()-1)
			cout << ' ';
		else
			cout << endl;
	}
	cout << endl;*/

	Solver solver;
	solver.init();
	if (argc == 1)
		solver.solve(P::solver_mode,cout);
	else
	{
		ofstream file(argv[1]);
		solver.solve(P::solver_mode,file);
	}
	return 0;
}
