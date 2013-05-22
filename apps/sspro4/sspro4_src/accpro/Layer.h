
#ifndef Layer_h
#define Layer_h 1
#include "General.h"
#include <math.h>



// Layer ver 2.01
// 9/8/2000
// 
// ANN Layers
// Linear, tanh and softmax outputs
// Categorical (one-hot), real-valued and mixed inputs.




class Layer
{
protected:
int NY;
int NU;
int NUr;
int* NK;

Float* Y;
Float* A;
Float* U;	 //NU*NK
Float* delta;    //NY
Float* backprop; //NU*NK

Float*** W;
Float*** dW;

Float* B;        //NY
Float* dB;       //NY

int level; 	//1=in,2=out


void alloc(int NY, int NU, int* NK);


public:

void softmax();
void squash();


// Constructor
// Categorical inputs

Layer(int t_NY, int* t_NK, int t_NU) :
	NY(t_NY), NU(t_NU)
{
NK=new int[NU];
for (int i=0; i<NU; i++)
	NK[i]=t_NK[i];
alloc(NY,NU,NK);
level=0;
}

// Constructor
// Real-valued inputs

Layer(int t_NY, int t_NU) :
	NY(t_NY), NU(t_NU)
{
NK=new int[NU];
for (int i=0; i<NU; i++)
	NK[i]=1;
NUr=0;
alloc(NY,NU,NK);
level=0;
}

// Constructor
// Mixed inputs (NU categorical attributes, NUr real-valued)

Layer(int t_NY, int* t_NK, int t_NU, int t_NUr) :
	NY(t_NY), NU(t_NU), NUr(t_NUr)
{
int i;
NK=new int[NU+NUr];
for (i=0; i<NU; i++)
	NK[i]=t_NK[i];
for (i=NU; i<NU+NUr; i++)
	NK[i]=1;
alloc(NY,NU+NUr,NK);
level=0;
}

Layer(istream& is);



void set_as_input() { 
  if (level==2) {level=3;}
  else {level=1;}
}

void set_as_output() { 
  if (level==1) {level=3;}
  else {level=2;}
}


void read(istream& is);
void write(ostream& os);

virtual void forward(int* I);
virtual void forward(Float* I);
virtual void forward(int* I1, Float* I2);
virtual void forward(Float* I1, Float* I2);

virtual Float f1(Float a);
virtual Float f_cost(Float* t);
Float log_cost(Float* t);
Float sq_cost(Float* t);


void gradient(int* I);
void gradient(Float* I);
void gradient(int* I1, Float* I2);
void gradient(Float* I1, Float* I2);
void gradient();


void resetGradient();
virtual void initWeights(int seed);

inline Float* back_out() { return backprop; }
inline Float* Aout() { return A; }
inline Float* out() { return Y; }


inline int get_NY() { return NY; }
inline int get_NU() { return NU; }
inline int* get_NK() { return NK; }

inline Float*** get_dW() { return dW; }

void set_dW(Float*** newdW);


};


class Layer_tanh : public Layer
{

public:


Layer_tanh(int t_NY, int* t_NK, int t_NU) :
Layer(t_NY, t_NK, t_NU)
{
}

Layer_tanh(int t_NY, int t_NU) :
Layer(t_NY, t_NU)
{
}

Layer_tanh(int t_NY, int* t_NK, int t_NU, int t_NUr) :
Layer(t_NY, t_NK, t_NU, t_NUr)
{
}

Layer_tanh(istream& is) :
Layer(is)
{
}

void forward(int* I)
{
Layer::forward(I);
squash();
}

void forward(Float* I)
{
Layer::forward(I);
squash();
}

void forward(int* I1,Float* I2)
{
Layer::forward(I1,I2);
squash();
}

void forward(Float* I1,Float* I2)
{
Layer::forward(I1,I2);
squash();
}





Float f1(Float a);

Float f_cost(Float* t);


void initWeights(int seed)
{
Layer::initWeights(seed);
}

};



class Layer_soft : public Layer
{

public :

Layer_soft(int t_NY, int* t_NK, int t_NU) :
Layer(t_NY, t_NK, t_NU)
{
}

Layer_soft(int t_NY, int t_NU) :
Layer(t_NY, t_NU)
{
}

Layer_soft(int t_NY, int* t_NK, int t_NU, int t_NUr) :
Layer(t_NY, t_NK, t_NU, t_NUr)
{
}

Layer_soft(istream& is) :
Layer(is)
{
}

void forward(int* I)
{
Layer::forward(I);
softmax();
}

void forward(Float* I)
{
Layer::forward(I);
softmax();
}

void forward(int* I1,Float* I2)
{
Layer::forward(I1,I2);
softmax();
}

void forward(Float* I1,Float* I2)
{
Layer::forward(I1,I2);
softmax();
}




Float f1(Float a);

Float f_cost(Float* t);

void initWeights(int seed)
{
Layer::initWeights(seed);
}

};




/*

class Selector : public Layer
{
int* edges;

public:


Selector(int t_NY, int* t_NK, int t_NU) :
Layer(t_NY, t_NK, t_NU)
{
edges=new int[NY];
}

Selector(int t_NY, int t_NU) :
Layer(t_NY, t_NU)
{
edges=new int[NY];
}

Selector(int t_NY, int* t_NK, int t_NU, int t_NUr) :
Layer(t_NY, t_NK, t_NU, t_NUr)
{
edges=new int[NY];
}

Selector(istream& is) :
Layer(is)
{
edges=new int[NY];
}


void forward(int* I)
{}
void forward(Float* I);
void forward(int* I1, Float* I2)
{}
void forward(Float* I1, Float* I2)
{}

Float f1(Float a);

Float f_cost(Float* t);


void setEdges();

void initWeights(int seed)
{
Layer::initWeights(seed);
setEdges();
}

};



*/

#endif
