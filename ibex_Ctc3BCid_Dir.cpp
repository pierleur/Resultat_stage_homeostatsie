#include "ibex_Ctc3BCid_Dir.h"

using namespace std;
namespace ibex {


Ctc3BCidDir::Ctc3BCidDir(Ctc& ctc, IntervalVector* boxE, int nb_para, int s3b, int scid,
			int vhandled, double var_min_width)
			:Ctc3BCid(ctc,s3b,scid,vhandled,var_min_width),m_boxE(boxE),m_nb_para(nb_para)
			{}

void Ctc3BCidDir::contract(IntervalVector& box) {
	ContractContext context(box);
	contract(box,context);
}


void Ctc3BCidDir::contract(IntervalVector& box, ContractContext& context) {
	int var;                                           // [gch] variable to be carCIDed
	m_initBox = IntervalVector(box);

	start_var=nb_var-1;
	//  patch pour l'optim  A RETIRER ??

	this->context = &context;
	

	for (int k=0; k<vhandled; k++) {                   // [gch] k counts the number of varCIDed variables [gch]

		var=(start_var+k)%nb_var;
		
		if ( borne(var) )//borneDown(var) || borneUp(var) )
		{
			var3BCID(box, var);
			
			if (box.is_empty()) {
				context.output_flags.add(FIXPOINT);
				this->context = NULL;
				return;
			}
		}
	}

	context.prop.update(BoxEvent(box,BoxEvent::CONTRACT));
	this->context = NULL;
	//	start_var=(start_var+vhandled)%nb_var;             //  en contradiction avec le patch pour l'optim
}

bool Ctc3BCidDir::borne(int var)
{	
	return m_initBox[var].ub() == (*m_boxE)[var].ub() || m_initBox[var].lb() == (*m_boxE)[var].lb();
}

}
