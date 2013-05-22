
#include "Ensemble.h"
#include "Sequence.h"
#include "fstream.h"

/*
	Predict the secondary structure for one sequence. Data set file only contains one sequence. 
*/
int
main(int argc, char** argv)
{
	if (argc!=5) 
	{
		cerr << "Usage: " << argv[0] << " model_definition dataset_file alignment_directory dataset_format\n";
		//dataset_foramt: 0: 3-line, 1: 9-line, 2: 12 line (9 line + predicted contacts of 8, 12, 16 A). 
		exit(1);
	}

	char model[256];
	char prot[256];
	char alig[256];
	char format[256]; 

	strcpy(model,argv[1]);
//	cout << "model definition: " << model << endl; 
	strcpy(prot,argv[2]);
//	cout << "data set: " << prot << endl; 
	strcpy(alig,argv[3]);
//	cout << "alignment directory: " << alig << endl; 
	strcpy(format, argv[4]);
//	cout << "dataset foramt: " << format << endl;

	
	int iFormat = atoi(format); 

	//prediction program accept five different formats
	/*
		0: 3-line (name, sequence, and true_secondary_structure. all of them are strings not seperated by space.
		1: 9-line
		2: 9 lines +  3 lines of contacts
		3: 9-lines + 3 lines of contacts + 1 line of predicted secondary structure
		4: line 1: name, line 2: sequence length line 3: sequence,  line 4-5: predicted contacts, line6: predicted secondary strucutre
		Notice: format 4 is only used by prediction program, not by training program and data set prediction.
	*/
  	if (iFormat != 0 && iFormat != 1 && iFormat != 2 && iFormat !=3 && iFormat != 4)
  	{
		cout << "data set format is not acceptable.\n"; 
		exit(1); 
  	}

//	cout << "read model file...." << endl; 
	//read models
	ifstream mstream(model);
	if (!mstream.is_open())
	{
		cout << "can't open model definition file.\n";
		exit(1); 	
	}

	Ensemble ensemble(mstream);
	
//	cout << "finish reading model file...." << endl; 
	
	//read test data set
	ifstream dstream(prot); 	
	if (!dstream.is_open())
	{
		cout << "can't open dataset file.\n";
		exit(1); 
	}

//	cout << "reading dataset...." << endl; 
	DataSet dataset(dstream, alig, iFormat);
	dataset.set_belief(0); 	

	//make prediction for one sequence only
	for (int i = 0; i < 1; i++)
	{
		ensemble.predict(dataset.sequence[i], 1);
	}
//	cout << endl; 
//	ensemble.printStatistics(); 
	return 0; 
}
