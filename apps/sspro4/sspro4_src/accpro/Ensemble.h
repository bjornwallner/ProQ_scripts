#include <math.h>
#include "Model.h"
#include "Sequence.h"

#ifndef _ENSEMBLE_H
#define _ENSEMBLE_H

#define METHOD 1   //1: average 2: vote 
//vote result is worse than average, thus remove voting code. 

class Ensemble
{

	private:
		int		m_modelNum; 
		int 		m_classNum; 
		Model**		m_pModelList; 	
		int		m_totalAminoAcids; 
		float *		m_errors; 
	public:
		Ensemble(istream& is); 
		~Ensemble(); 

		//called by PredictSet
		void predict(Sequence* pSeq, bool detailed = false);

		//called by PredictSeq, class index start from 0 to m_classNum - 1
		void predict(Sequence* pSeq, int classIndex);
		void resetStatistics(); 
		void printStatistics(); 
};

#endif 
