#ifndef _ENSEMBLE_H
#define _ENSEMBLE_H

#include "Model.h"
#include "Sequence.h"

class Ensemble
{

private:
	int m_modelNum;
	int m_classNum; 

	Model ** m_pModelList; 

	int m_totalAminoAcids;
	int m_helixNum;
	int m_sheetNum;
	int m_coilNum;

	int m_totalErrors;
	int m_helixErrors;
	int m_sheetErrors;
	int m_coilErrors; 

public:
	Ensemble(istream& is);

	~Ensemble();

	void predict(Sequence* seq, bool detailed = false);
	float getHelixRecall();
	float getSheetRecall();
	float getCoilRecall();	
	float getRecall();
	void resetStatistics(); 
	void printStatistics(); 

};

#endif
