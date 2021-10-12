#ifndef __SOLVE_PROB_H__
#define __SOLVE_PROB_H__

#include <iostream>
#include <ibex.h>
#include "ibex_Ctc3BCid_Dir.h"
#include <fstream>
#include <list>
#include <cmath>

using namespace std;
using namespace ibex;

struct Matrice {
	vector<double> value;
	int x,y;
};

class Solve_prob
{
public:
	Solve_prob(Ctc& ctc, IntervalVector& boxE, int nb_para, System& system);
	void contract(double prec);
	void contract(int nb_iter);
protected:
	void contraction(int dim);
	void mise_a_jour_boxI();
	int contact_with_E(IntervalVector& current_box);
	void HullSolution(IntervalVector& box,DefaultSolver& solver);
	
	void genereBoxI();
	void cutWithBoxI();
	
	IntervalVector& m_boxE;
	IntervalVector  m_boxI;
	
	int m_nb_para;
	
	Ctc3BCidDir m_ctc;
	Ctc3BCid    m_ctc_3B;
	list<IntervalVector> m_list_segments;

	System& m_system;
	NormalizedSystem m_Nsystem;
	int m_test;
};

#endif
