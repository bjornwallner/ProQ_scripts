
#ifndef Model_h
#define Model_h 1

#include <stdlib.h>
#include <math.h>
#include "Sequence.h"
#include "NN.h"
#include "NNt.h"


class Model {
 private:
  int NU;
  int NY;
  int NH;

  int context;
  int Moore; 

  int NF;
  int NB;
  int NH2;

  int CoF;
  int CoB;
  int Step;

  Float threshold;

  int modular;


  NN* NetYOne;

//  NN* NetOut;
  NNt* NetOut;
  NNt* NetF;
  NNt* NetB;

  Float** FF;
  Float** BB;
  Float** FFbp;
  Float** BBbp;

  Float* P_F;
  Float* P_B;

  Float* Y;

  int** Conf;

int classes;

  Float temp_error;
  int temp_aas;
  
  Float squared_error; 
  Float squared_errorF; 
  Float squared_errorB; 
  int* nerrors;
  int* nerrors_0;
  int* nerrors_1;
  int* nall_0;
  int* nall_1;

  Float epsilon;

  void alloc();


 public:

  Model(int NU, int NY, int NH,  int context ,int Moore, int NF, int NB, int NH2, int CoF, int CoB, int Step, Float threshold=0.0, int classes=1);
  Model(istream& is);
  void read(istream& is);
  void write(ostream& os);

  void randomize(int seed);


  void F1_F(Sequence* seq, int t);
  void B1_B(Sequence* seq, int t);
  void propagate(Sequence* seq);
  void propagate(Sequence* seq, int T, int W);

  void forward(Sequence* seq, int t);

  void Feed(Sequence* seq);
  void predict(Sequence* seq);
  void predict(Sequence* seq, int W);
  Float* out() {return Y;}
  int** getConf() {return Conf;}

  int* getNErrors() { return nerrors;};

  int* getNErrors_0() { return nerrors_0;};
  int* getNErrors_1() { return nerrors_1;};

  int* getNAll_0() { return nall_0;};
  int* getNAll_1() { return nall_1;};



   void resetNErrors() { 

	   for (int c=0;c<classes;c++) {
		   nerrors[c]=0;
		   nerrors_0[c]=0;
		   nerrors_1[c]=0;
		   nall_0[c]=0;
		   nall_1[c]=0;
	   }


	for (int p=0;p<NY;p++)
	  for (int y=0;y<NY;y++)
		Conf[p][y]=0;
	};

  Float get_tot_error() { 
	return -temp_error;
	};
  Float get_squared_error() { 
	return -squared_error;
	};
  Float get_squared_errorF() { 
	return squared_errorF;
	};
  Float get_squared_errorB() { 
	return squared_errorB;
	};
  void reset_squared_error() { 
	temp_error=0.0;
	squared_error=0.0;
	squared_errorF=0.0;
	squared_errorB=0.0;
	};

  void setEpsilon(Float eps) { epsilon=eps; };


};


#endif // Model_h
