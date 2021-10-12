#ifndef __IBEX_CTC_3B_CID_DIR_H__
#define __IBEX_CTC_3B_CID_DIR_H__

#include <ibex.h>

namespace ibex {

class Ctc3BCidDir : public Ctc3BCid
{

public:
	
	Ctc3BCidDir(Ctc& ctc, IntervalVector* boxE,int nb_para, int s3b=default_s3b, int scid=default_scid,
			int vhandled=-1, double var_min_width=default_var_min_width);
	
	void contract(IntervalVector& box) override;
	
	virtual void contract(IntervalVector& box, ContractContext& context);
	
protected:
	//Boite enveloppe
	IntervalVector* m_boxE;
	
	bool borne(int var);
private:
	IntervalVector  m_initBox;
	int m_nb_para;
};

}

#endif
