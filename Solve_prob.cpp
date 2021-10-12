#include "Solve_prob.h"

Solve_prob::Solve_prob(Ctc& ctc, IntervalVector& boxE, int nb_para, System& system) :
	m_ctc(ctc,&boxE,nb_para),m_boxE(boxE),m_nb_para(nb_para),
	m_boxI(boxE.size(),Interval(1,0)),m_ctc_3B(ctc),m_system(system),m_Nsystem(system,1e-7)
{
	m_test=0;
}

void Solve_prob::HullSolution(IntervalVector& box,DefaultSolver& solver)
{
	std::pair<IntervalVector, double> result;
	LoupFinderDefault find(m_Nsystem);
	IntervalVector tempo(box.size());
	IntervalVector vide(box.size(),Interval(1,0));
	for (int i=0;i<solver.get_data().nb_solution();i++)
	{
		try{
			tempo = solver.get_data().solution(i);
			result = find.find(tempo,vide,1e9);
			box |= result.first;
		}catch(...){}
	}
	for (int i=0;i<solver.get_data().nb_unknown();i++)
	{
		try{
			tempo = solver.get_data().unknown(i);
			result = find.find(tempo,vide,1e9);
			box |= result.first;
		}catch(...){}
	}
}

void Solve_prob::contract(double prec)
{
	m_ctc.contract(m_boxE);
	m_list_segments.push_back(m_boxE);
	
	int n=0;
	while (m_boxE.perimeter()>m_boxI.perimeter()*(1+prec)){
		contraction(n%m_boxE.size());	
		n++;	
	}
	cout <<endl <<"boxI :\n"<< m_boxI<<endl;
}

void Solve_prob::contract(int nb_iter)
{
	m_ctc.contract(m_boxE);
	m_list_segments.push_back(m_boxE);
	
	int n=0;
	while (n<nb_iter){
		contraction(n%m_boxE.size());	
		n++;	
	}
	cout <<endl <<"boxI :\n"<< m_boxI<<endl;
}

int Solve_prob::contact_with_E(IntervalVector& current_box)
{
	for (int i=0;i<m_boxE.size()-m_nb_para;i++)
		if ( (current_box[i].lb()==m_boxE[i].lb()) || (current_box[i].ub()==m_boxE[i].ub()))
			return i+1;
	return 0;
}

void Solve_prob::genereBoxI()
{
	auto current = m_list_segments.begin();
	int taille = m_list_segments.size();
	
	int id_touch;

	if (m_nb_para>0)
	for (;current != m_list_segments.end();++current)
	{
		if (id_touch=contact_with_E(*current))
		{
			id_touch--;
			IntervalVector tempo = *current;
			
			for (int k=m_boxE.size()-m_nb_para;k<m_boxE.size();k++)
			{
				double middle = ( tempo[k].lb()+tempo[k].ub() )/2.;
				tempo[k] = Interval(middle,middle);
			}
			
			DefaultSolver solver(m_system,tempo.extr_diam_index(false)/100.);
			solver.solve(tempo);
			HullSolution(m_boxI,solver);
		}
	}
	
	for (int i=m_boxE.size()-m_nb_para;i<m_boxE.size();i++)
		m_boxI[i] = m_boxE[i];
}

void Solve_prob::cutWithBoxI()
{
	IntervalVector *result;
	for (auto id = m_list_segments.begin(); id != m_list_segments.end(); ++id)
	{
		int n=id->diff(m_boxI,result);
		
		if (n==0)
		{
			auto memo = id;
			--memo;
			m_list_segments.erase(id);
			id = memo;
		}
		else
		{
			if (n==1)
				*id = result[0];
		}
			
			/*
		for (int i=1; i<n; i++) {
			m_list_segments.push_front(result[i]);
		}*/
		
		delete[] result;
	}
}

void Solve_prob::contraction(int dim)
{
	auto current = m_list_segments.begin();
	int taille = m_list_segments.size();
		
	for (int i = 0;i<taille;i++)
	{
		if (contact_with_E(*current))
		{
			if ((*current)[dim].is_bisectable())
			{
				pair<IntervalVector,IntervalVector> subboxes = current->bisect(dim);
				m_ctc.contract(subboxes.first );
				m_ctc.contract(subboxes.second);
			    	*current = subboxes.first;
			    	m_list_segments.push_back(subboxes.second);
		    	}
		}
		++current;
	}
	
	genereBoxI();
	cutWithBoxI();
	
	current = m_list_segments.begin();
	IntervalVector newE(*current);
	++current;
	for (;current != m_list_segments.end();++current)
		newE |= *current;
	m_boxE = newE;
	
	cout << "Taille : " << m_list_segments.size() << " perimetre : " << m_boxE.perimeter() << " BoxI: " << m_boxI.perimeter() << endl;
}
