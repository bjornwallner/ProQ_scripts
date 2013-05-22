

#ifndef NN_h
#define NN_h 1
#include "Layer.h"
#include "General.h"
#include <math.h>

// NN ver. 2.01 (8/9/2000)
//
// One-hidden layer Feedforward neural net.
// Input: categorical data (one-hot), real valued or mixed. 
// ouput: multinomial probability P(class|attributes) with normalized
// exponentials (softmax).
// Cost function: -log likelihood for the multinomial (cross entropy)
// Gradient: plain backpropagation (no momentum)

class NN 
{
private:
  int NI;
  int NIr;
  int NH;
  int NO;
  int* NK;
  int* NK2;

  Float* backprop;

  Layer_soft* upper;
  Layer_tanh* lower;


public:
  // Constructor. Parameters:
  // Number of input attributes, number of hidden units, number of output units;
  // t_NK contains the cardinalities of the attribute spans.

  NN(int t_NI, int t_NH, int t_NO, int* t_NK) :
    NI(t_NI), NH(t_NH), NO(t_NO)
    {
      NK=new int[NI];
      for (int i=0; i<NI; i++) NK[i]=t_NK[i];
      upper= new Layer_soft(NO,NH);
      upper->set_as_output();
      lower= new Layer_tanh(NH,NK,NI);
      lower->set_as_input();
      NIr=0;
    };

  // Constructor for a net with mixed inputs.
  // NI = number of input attributes (categorical inputs)
  // NIr = number of inputs (real valued)
  // ..
  // outp = output or non-output network (for backprop signal)

  NN(int t_NI,int t_NIr, int t_NH, int t_NO, int* t_NK, int outp=1) :
    NI(t_NI), NIr(t_NIr), NH(t_NH), NO(t_NO)
    {
      int i;
      NK=new int[NI];
      for (i=0; i<NI; i++)
	NK[i]=t_NK[i];
      NK2=new int[NIr];
      for (i=0; i<NIr; i++)
	NK2[i]=1;
      upper= new Layer_soft(NO,NH);
      if (outp) 
        upper->set_as_output();
      lower= new Layer_tanh(NH,NK,NI,NIr);
//      lower->set_as_input();
      backprop=new Float[NIr];
    };

  // Create/read a net from file
  NN(istream& is, int outp=1);
  void read(istream& is, int outp=1);

  // Forward pass
  void forward(int* I);
  void forward(Float* I);
  void forward(int* I1, Float* I2);
  void forward(Float* I1, Float* I2);


  Float f_cost(Float* t) {
	return upper->f_cost(t);
  }

  // Update gradients 
  void gradient(int* I, Float* t);
  void gradient(Float* I, Float* t);
  void gradient(int* I1, Float* I2, Float* t);
  void gradient(Float* I1, Float* I2, Float* t);

  // Update weights
  void resetGradient();
  virtual void initWeights(int seed);
  inline Float* out() { return upper->out(); };
  void write(ostream& os);

  inline int get_NI() { return NI; };
  inline int get_NIr() { return NIr; };
  inline int get_NO() { return NO; };
  inline int get_NH() { return NH; };

};



#endif // NN_h
