
#ifndef Model_h
#define Model_h 1

#include <stdlib.h>
#include <math.h>
#include "Sequence.h"
#include "NN.h"
#include "NNt.h"

#define INPUT_SIZE 1024 

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
  Float m_input[INPUT_SIZE];
  Float m_3dinput[INPUT_SIZE]; 

  int modular;


  NN* NetYOne;

  NN* NetOut;
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

  Float temp_error;
  int temp_aas;
  
  Float squared_error; 
  Float squared_errorF; 
  Float squared_errorB; 
  int nerrors;
  int nerrors_alpha;
  int nerrors_beta;
  int nerrors_gamma;

  Float epsilon;

  void alloc();

  int m_use3d; 
  int m_useacc; 
  int m_binNum; 
  bool m_full3d; 

  private:
   Float calcDistance(Float x1, Float y1, Float z1, Float x2, Float y2, Float z2);

 public:
   //use3d: indicate where or not the model will use 3d information. default value is 0 (not use 3d information). 
  Model(int NU, int NY, int NH,  int context ,int Moore, int NF, int NB, int NH2, int CoF, int CoB, int Step, Float threshold=1000.0, int use3d = 0, int useacc = 0, int binNum = 0);
  Model(istream& is);
  void read(istream& is);
  void write(ostream& os);

  void randomize(int seed);


  void F1_F(Sequence* seq, int t);

  Float * prepareInput(Sequence* seq, int t); 
  Float * prepareInputWith3d(Sequence * seq, int t); 

  void B1_B(Sequence* seq, int t);
  void propagate(Sequence* seq);
  void propagate(Sequence* seq, int T, int W);

  void forward(Sequence* seq, int t);


  void Feed(Sequence* seq);
  void predict(Sequence* seq);
  void predict(Sequence* seq, int W);
  Float* out() {return Y;}
  int** getConf() {return Conf;}

  int getNErrors() { return nerrors;};

  int getNErrors_alpha() { return nerrors_alpha;};
  int getNErrors_beta() { return nerrors_beta;};
  int getNErrors_gamma() { return nerrors_gamma;};

   void resetNErrors() { 
	nerrors=0;
	nerrors_alpha=0;
	nerrors_beta=0;
	nerrors_gamma=0;
	for (int p=0;p<NY;p++)
	  for (int y=0;y<NY;y++)
		Conf[p][y]=0;
	};

  Float get_tot_error() { 
	return temp_error;
	};
  Float get_squared_error() { 
	return squared_error;
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

  //Lenght of Y is hard coded: 8192
  Float * getOutput() { return Y; }; 


};


#endif // Model_h
