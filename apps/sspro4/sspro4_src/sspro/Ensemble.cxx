#include "Ensemble.h"
/**
	read a list of model from files.
	format of model definition file:
	first line: model_number class_number 
	following lines: the full path to the model file
*/ 
Ensemble::Ensemble(istream& is) : m_modelNum(0), m_pModelList(NULL), m_totalAminoAcids(0),
	m_helixNum(0), m_sheetNum(0), m_coilNum(0), m_totalErrors(0), m_helixErrors(0), 
	m_sheetErrors(0), m_coilErrors(0)
{
	is >> m_modelNum >> m_classNum; 	
	m_pModelList = new Model*[m_modelNum]; 

	char filename[200]; 
	for ( int i = 0; i < m_modelNum; ++i)
	{
		is >> filename;
		ifstream fin(filename);
		if (!fin.is_open())
		{
			cout << "can't open the model file, program died.\n";
			exit(1); 
		}			
		m_pModelList[i] = new Model(fin); 
	}
	
}

Ensemble::~Ensemble()
{
}

void Ensemble::predict(Sequence* seq, bool detailed)
{
	int length = seq->length; 
	float * results = new float[3 * length]; 
	memset(results, 0, 3 * length * sizeof(float)); 
	
	//take the average output
	for (int i = 0; i < m_modelNum; ++i)
	{
		m_pModelList[i]->predict(seq); 
		float * pOutput = m_pModelList[i]->getOutput();
		for (int j = 0; j <= 3*length - 1; j++)
		{
			results[j] += pOutput[j] / m_modelNum; 
		}
	} 	

	//make prediction	
	int t; 
	for (t=1; t<=seq->length; t++) 
	{
  		if (seq->y[t]<0) 
    		{
			seq->y_pred[t]=-1;
    			continue;
    		}

		Float pred=0.0;
		int label=-1;

		seq->alpha[t] = results[3*(t-1)];
		seq->beta[t]  = results[3*(t-1)+1];
		seq->gamma[t] = results[3*(t-1)+2];

		for (int c=0; c<3; c++) 
		{
			if (results[3*(t-1)+c] > pred) 
			{
				pred = results[3*(t-1)+c];
				label=c;
			}
		}
		seq->y_pred[t]=label;
	}

	if (detailed)
	{
		//print out true secondary structure
		for (t=1; t<=seq->length; t++) 
		{
			if (seq->y[t]==0) cout << "H";
			if (seq->y[t]==1) cout << "E";
			if (seq->y[t]==2) cout << "C";
		}
		cout << endl; 
	}

	for (t=1; t<=seq->length; t++) 
	{
		m_totalAminoAcids++; 
		switch (seq->y[t])
		{
			case 0: m_helixNum++; break;
			case 1: m_sheetNum++; break;
			case 2: m_coilNum++; break;
			default: break;
		}

		if (seq->y[t]!=seq->y_pred[t])
		{
			m_totalErrors++; 
			if (seq->y[t]==0) m_helixErrors++;
			if (seq->y[t]==1) m_sheetErrors++;
			if (seq->y[t]==2) m_coilErrors++;
		}
		//output the prediction details
		if (detailed)
		{
			if (seq->y_pred[t]==0) cout << "H";
			if (seq->y_pred[t]==1) cout << "E";
			if (seq->y_pred[t]==2) cout << "C";
		}
	}
	if (detailed)
	{
		cout << endl; 
	}
	delete [] results;
}

float Ensemble::getHelixRecall()
{
	return (float)(m_helixNum - m_helixErrors) / m_helixNum; 
}
float Ensemble::getSheetRecall()
{
	return (float)(m_sheetNum - m_sheetErrors) / m_sheetNum; 
}
float Ensemble::getCoilRecall()
{
	return (float)(m_coilNum - m_coilErrors) / m_coilNum; 
}

float Ensemble::getRecall()
{
	return (float)(m_totalAminoAcids - m_totalErrors) / m_totalAminoAcids; 
}
void Ensemble::resetStatistics()
{
	m_totalAminoAcids = m_helixNum = m_sheetNum = m_coilNum = 0; 
	m_totalErrors = m_helixErrors = m_sheetErrors = m_coilErrors  = 0; 
}

void Ensemble::printStatistics()
{
	cout << "Overall error rate: " << 1 - getRecall() << " (" << m_totalErrors << "/" << m_totalAminoAcids << ")\n";
	cout << "helix error rate: " << 1 - getHelixRecall() << " (" << m_helixErrors << "/" << m_helixNum << ")\n";
	cout << "beta sheet error rate: " << 1 - getSheetRecall() << " (" << m_sheetErrors << "/" << m_sheetNum << ")\n";
	cout << "coil error rate: " << 1 - getCoilRecall() << " (" << m_coilErrors << "/" << m_coilNum << ")\n";
} 
