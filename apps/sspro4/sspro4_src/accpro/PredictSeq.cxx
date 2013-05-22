#include <math.h>
#include "Sequence.h"
#include "Ensemble.h"

/*
	Use ensemble to predict the solvent acc for one sequence. 
	Model Definition File Format: line 1: model_number class_number, line 2..., the path to each model. 
	Data Set Format: 1: 9-line, 2: two line(name, sequence in compact form), in addition to that, the dataset must have a title line (e.g. 1 20 3)
	12/18/2003, Jianlin Cheng
*/


int
main(int argc, char** argv)
{
	if (argc!=6) 
	{
		cerr << "Usage: " << argv[0] << " model_definition dataset_file alignment_directory dataset_format threshold_index\n";
		//dataset_foramt: 1: 9-line format. 
		exit(1);
	}

	char model[256];
	char prot[256];
	char alig[256];
	char format[256]; 
	char classIndex[10]; 

	strcpy(model,argv[1]);
//	cout << "model definition: " << model << endl; 
	strcpy(prot,argv[2]);
//	cout << "data set: " << prot << endl; 
	strcpy(alig,argv[3]);
//	cout << "alignment directory: " << alig << endl; 
	strcpy(format, argv[4]);
//	cout << "dataset foramt: " << format << endl;
	strcpy(classIndex, argv[5]);
//	cout << "threshold index: " << classIndex << endl;

	int iFormat = atoi(format); 
	if (iFormat != 1 && iFormat != 2)
	{
		cout << "dataset must be 9-line or 2-line format.\n";
		exit(1); 
	}

	//acc thredhold = index * 5%, support 20 different thresholds. 
	//e.g.: index(threshold): 0(0), 1(5%), 2(10%), 3(15%), ... 19(95%) 
	int index = atoi(classIndex); 
	if (index < 0)
	{
		cout << "threshold index must be equal or bigger than 0. \n";
		exit(1); 
	}
  

	//read models
	ifstream mstream(model);
	if (!mstream.is_open())
	{
		cout << "can't open model definition file.\n";
		exit(1); 	
	}
	Ensemble ensemble(mstream);
//	cout << "models are read from file. \n"; 

 // 	cout << "Reading test dataset\n";
  	ifstream tstream(prot);
	if (!tstream.is_open())
	{
		cout << "can't open dataset file.\n";
		exit(1); 
	}

  	DataSet dataset(tstream, alig, iFormat);
  	dataset.set_belief(0);


	//make prediction
	for (int i = 0; i < 1; i++)
	{
		ensemble.predict(dataset.sequence[i], index);
	}
  	return 0;
}
