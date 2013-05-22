#include <math.h>
#include "Model.h"
#include "Sequence.h"
#include "Ensemble.h"



Ensemble::Ensemble(istream& is): m_totalAminoAcids(0)
{
	is >> m_modelNum >> m_classNum;
	char modelName[500];
	m_pModelList = new Model*[m_modelNum]; 
	m_errors = new float[m_classNum]; 
	for (int j = 0; j < m_classNum; j++)
	{
		m_errors[j] = 0; 
	}
	for ( int i = 0; i < m_modelNum; i++)
	{
		is >> modelName; 
		ifstream mstream(modelName);
		if (!mstream.is_open())
		{
			cout << "can't read model: " << modelName << endl; 
			exit(1); 
		}
		m_pModelList[i] = new Model(mstream); 
		
	}
	
}

Ensemble::~Ensemble()
{
	delete [] m_errors; 
	for (int i = 0; i < m_modelNum; i++)
	{
		delete m_pModelList[i]; 
	}
	delete [] m_pModelList; 
}

void Ensemble::resetStatistics()
{
	m_totalAminoAcids = 0; 
	for (int k = 0; k < m_modelNum; k++)
	{
		m_pModelList[k]->resetNErrors(); 
	}
	for (int i = 0; i < m_classNum; i++)
	{
		m_errors[i] = 0; 
	}
}

void Ensemble::printStatistics()
{
	cout << "Prediction error rate for each class: " << endl; 
	
	for (int i = 0; i < m_classNum; i++)
	{
		cout << "Threshold " << i  << "=" << i*0.05 << " : "; 
		cout << m_errors[i] / m_totalAminoAcids << endl; 	
	}
}

void Ensemble::predict(Sequence* pSeq, int classIndex)
{
    int length = pSeq->length;

    //new buffer to store average result
    Float *pResult = new Float[length * m_classNum];

    //new buffer to store result of each model
    //#define PFLOAT Float*
    Float ** results = new Float*[m_modelNum];
    //Float ** results = new PFLOAT[m_modelNum];
    for(int m=0; m < m_modelNum; ++m)
    {
         results[m] = new Float[length * m_classNum];
    	memset(results[m], 0, length*m_classNum*sizeof(Float));
    }

    memset(pResult, 0, length*m_classNum*sizeof(Float));

    for (int i = 0; i < m_modelNum; i++)
    {
       m_pModelList[i]->predict(pSeq);
       Float* pModelResult = m_pModelList[i]->out();
       for (int j = 0; j < length * m_classNum; ++j)
       {
	       pResult[j] += pModelResult[j];
               results[i][j]= pModelResult[j];

	       if (i == m_modelNum -1) //average the result
	       {
		       pResult[j] /= m_modelNum;
	       }
       }
    }
	
    //calculate error number
    for (int t =1; t <= length; ++t)
    {
            if (pSeq->y[t] < 0)
            {
                continue;
            } 
            ++m_totalAminoAcids;
	    int * predict = new int[m_classNum];
	    for (int c =0; c < m_classNum; c++)
	    {
                    //predict using average 
                    if (pResult[m_classNum*(t-1)+c]  > 0.5)
		    {
			    if (c == classIndex)
				{
					cout << 'e'; 
				}
			    predict[c] = 1;
		    }
		    else
		    {
			    if (c == classIndex)
				{
					cout << 'b'; 
				}
			    predict[c] = 0;
		    }

	    }
	    delete [] predict;
    } 
	cout << endl; 

    delete[] pResult;
    for(int m=0; m < m_modelNum; ++m)
    {
         delete  [] results[m];
    }
    delete [] results;
}

void Ensemble::predict(Sequence* pSeq, bool detailed)
{
    int length = pSeq->length;

    //new buffer to store average result
    Float *pResult = new Float[length * m_classNum];

    //new buffer to store result of each model
    Float ** results = new Float*[m_modelNum];
    for(int m=0; m < m_modelNum; ++m)
    {
         results[m] = new Float[length * m_classNum];
    	memset(results[m], 0, length*m_classNum*sizeof(Float));
    }

    memset(pResult, 0, length*m_classNum*sizeof(Float));

    for (int i = 0; i < m_modelNum; i++)
    {
       m_pModelList[i]->predict(pSeq);
       Float* pModelResult = m_pModelList[i]->out();
       for (int j = 0; j < length * m_classNum; ++j)
       {
	       pResult[j] += pModelResult[j];
               results[i][j]= pModelResult[j];

	       if (i == m_modelNum -1) //average the result
	       {
		       pResult[j] /= m_modelNum;
	       }
       }
    }

    //calculate error number
    float step = 1.0/m_classNum;
    for (int t =1; t <= length; ++t)
    {
            if (pSeq->y[t] < 0)
            {
                continue;
            } 
            ++m_totalAminoAcids;
	    int * predict = new int[m_classNum];
	    int * target = new int[m_classNum];
	    for (int c =0; c < m_classNum; c++)
	    {
                    //predict using average 
                    if (pResult[m_classNum*(t-1)+c]  > 0.5)
		    {
			    predict[c] = 1;
		    }
		    else
		    {
			    predict[c] = 0;
		    }

		    if ( pSeq->racc[t] <= step * c)
		    {
			    target[c] = 0;
		    }
		    else
		    {
			    target[c]= 1;
		    }
		    if (target[c] != predict[c])
		    {
			    ++m_errors[c];
		    }
	    }
	    delete [] predict;
	    delete [] target;
    } 

    delete[] pResult;
    for(int m=0; m < m_modelNum; ++m)
    {
         delete  [] results[m];
    }
    delete [] results;

}

