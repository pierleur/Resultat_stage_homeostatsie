#include <iostream>
#include <chrono>

#include "ibex.h"
#include <math.h>

#include "ibex_Ctc3BCid_Dir.h"
#include "Solve_prob.h"
#include <list>


using namespace std;
using namespace ibex;



int main()
{	
	int nb_para=2;
	
	std::chrono::high_resolution_clock::time_point a,b;
	unsigned int time;
	
	System sys("FeinbergShinar2_cl.txt");
	DefaultSolver solver(sys,1e-07);
	
	CtcHC4 hc44cid(sys.ctrs,1e-3,true);
	
	IntervalVector result_Y2(sys.box);
	
	Ctc3BCid exe2(hc44cid);
	
	Solve_prob solv(exe2,result_Y2,nb_para,sys);
	
	a = std::chrono::high_resolution_clock::now();
	solv.contract(1e-2);
	b = std::chrono::high_resolution_clock::now();
	
	
	time = chrono::duration_cast<std::chrono::microseconds>(b - a).count();
	cout << "Resultat :\n" << result_Y2 << endl;
	cout << "Time s : " << time/1000000. << " volume : "<<result_Y2.subvector(0,sys.box.size()-nb_para-1).volume()<<endl<<endl;
	cout << "Perimetre : " << result_Y2.subvector(0,sys.box.size()-nb_para-1).perimeter() << endl<<endl;

	return 0;
}
