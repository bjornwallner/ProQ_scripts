
// NNt ver. 2.01 (8/9/2000)

#include "NNt.h"


NNt::NNt(istream& is, int outp)
{
  is >> NO >> NH >> NI >> NIr;
  upper=new Layer_tanh(is);
  if (outp)
    upper->set_as_output();
  lower=new Layer_tanh(is);
//  lower->set_as_input();

  int i;
  NK=new int[NI];
  for (i=0; i<NI; i++) 
	NK[i]=lower->get_NK()[i];
  NK2=new int[NIr];
  for (i=0; i<NIr; i++) 
	NK2[i]=1;
  backprop=new Float[NIr];
}


void
NNt::read(istream& is, int outp)
{
  is >> NO >> NH >> NI >> NIr;
  upper->read(is);
  if (outp)
  	upper->set_as_output();
  lower->read(is);
//  lower->set_as_input();

  int i;
  for (i=0; i<NI; i++) 
	NK[i]=lower->get_NK()[i];

  for (i=0; i<NIr; i++) 
	NK2[i]=1;
}



void
NNt::forward(int* I)
{
  lower->forward(I);
  upper->forward(lower->out());
}

void
NNt::forward(Float* I)
{
  lower->forward(I);
  upper->forward(lower->out());
}

void
NNt::forward(int* I1,Float* I2)
{
  lower->forward(I1,I2);
  upper->forward(lower->out());
}
void
NNt::forward(Float* I1,Float* I2)
{
  lower->forward(I1,I2);
  upper->forward(lower->out());
}

void
NNt::gradient(int* I, Float* t)
{
  upper->gradient();
  lower->gradient(I);
}

void
NNt::gradient(Float* I, Float* t)
{
  upper->gradient();
  lower->gradient(I);
}

void
NNt::gradient(int* I1,Float* I2, Float* t)
{
  upper->gradient();
  lower->gradient(I1,I2);
}
void
NNt::gradient(Float* I1,Float* I2, Float* t)
{
  upper->gradient();
  lower->gradient(I1,I2);
}


void
NNt::resetGradient()
{
  lower->resetGradient();
  upper->resetGradient();
}

void
NNt::initWeights(int seed)
{
  lower->initWeights(seed);
  upper->initWeights(seed);
}


void
NNt::write(ostream& os)
{
  os << NO << " " << NH<< " " << NI<< " " << NIr <<"\n";
  upper->write(os);
  lower->write(os);
}



