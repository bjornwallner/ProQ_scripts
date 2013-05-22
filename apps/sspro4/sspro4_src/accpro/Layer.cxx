
// Layer ver 2.01
// 9/8/2000

#include "Layer.h"
#include <stdlib.h>

void
Layer::softmax() {
  int y;

  int overflow=0;
  Float max=A[0];
  int amax=0;
  Float norm=0;
  for (y=0; y<NY; y++) {
    if (A[y]>85) {
      overflow=1;
      }
    else {
      norm += (Float)exp(A[y]);
      }
    if (A[y]>max) {
      max = A[y];
      amax = y;
      }
    }

  if (overflow) {
    memset(Y, 0, NY*sizeof(Float));
    Y[amax]=1.0;
    return;
    }
  else {
    for (y=0; y<NY; y++) {
      Y[y] = (Float)exp(A[y])/norm;
      }
    }
}



void
Layer::squash() {
for (int y=0; y<NY; y++)
	Y[y]=(Float)tanh(A[y]);
}






void
Layer::alloc(int NY, int nu, int* NK)
{
int y,u;
int NUtot=0;

for (u=0; u<nu; u++)
	NUtot += NK[u];

Y=new Float[NY];
A=new Float[NY];
U=new Float[NUtot];

delta=new Float[NY];
backprop=new Float[NUtot];

W=new Float**[NY];
dW=new Float**[NY];
B=new Float[NY];
dB=new Float[NY];
for (y=0; y<NY; y++) {
	W[y]=new Float*[nu];
	dW[y]=new Float*[nu];
	for (u=0; u<nu; u++) {
		W[y][u]=new Float[NK[u]];
		dW[y][u]=new Float[NK[u]];
		}
	}

}




Layer::Layer(istream& is)
{
int y,u,k;

is >> NY;
is >> NU;
is >> NUr;
NK=new int[NU+NUr];

for (u=0; u<NU+NUr; u++) is >> NK[u];

alloc(NY,NU+NUr,NK);

for (y=0; y<NY; y++) {
  for (u=0; u<NU+NUr; u++) {
    for (k=0; k<NK[u]; k++) {
      is >> W[y][u][k];
      }
    }
  is >> B[y];
  }
level=0;
}



void
Layer::read(istream& is)
{
int y,u,k;

is >> NY;
is >> NU;
is >> NUr;

for (u=0; u<NU+NUr; u++) is >> NK[u];

for (y=0; y<NY; y++) {
  for (u=0; u<NU+NUr; u++) {
    for (k=0; k<NK[u]; k++) {
      is >> W[y][u][k];
      }
    }
  is >> B[y];
  }
level=0;
}





void Layer::write(ostream& os)
{
int y,u,k;

os << NY << "\n";
os << NU << "\n";
os << NUr << "\n";

for (u=0; u<NU+NUr; u++)
	os << NK[u] << " ";
os << "\n";

for (y=0; y<NY; y++) {
  for (u=0; u<NU+NUr; u++) {
    for (k=0; k<NK[u]; k++) {
      os << W[y][u][k] << " ";
      }
    }
  os << B[y] << "\n";
  }
}




void
Layer::forward(int* I)
{
int y,u;

  for (y=0; y<NY; y++) {
    Float a=B[y];
    for (u=0; u<NU; u++) {
      if (I[u]>=0)
        a += W[y][u][I[u]];
      }
    Y[y]=A[y]=a;
    }
}




void
Layer::forward(Float* I)
{
int y,u,k,i;

  for (y=0; y<NY; y++) {
    i=0;
    Float a=B[y];
    for (u=0; u<NU; u++) {
      for (k=0; k<NK[u]; k++) {
        U[i]=I[i];
        a += W[y][u][k]*I[i++];
        }
      }
    Y[y]=A[y]=a;
    }
}


void
Layer::forward(int* I1, Float* I2)
{
int y,u,k,i;

  for (y=0; y<NY; y++) {
    i=0;
    Float a=B[y];
    for (u=0; u<NU; u++) {
      if (I1[u]>=0)
        a += W[y][u][I1[u]];
      }
    for (u=NU; u<NU+NUr; u++) {
      for (k=0; k<NK[u]; k++) {
        U[i]=I2[i];
        a += W[y][u][k]*I2[i++];
        }
      }
    Y[y]=A[y]=a;
    }
}

void
Layer::forward(Float* I1, Float* I2)
{
int y,u,k,i1,i2;

  for (y=0; y<NY; y++) {
    i1=0;
    i2=0;
    Float a=B[y];
    for (u=0; u<NU; u++) {
      for (k=0; k<NK[u]; k++)
        a += W[y][u][k]*I1[i1++];
      }
    for (u=NU; u<NU+NUr; u++) {
      for (k=0; k<NK[u]; k++) {
        U[i2]=I2[i2];
        a += W[y][u][k]*I2[i2++];
        }
      }
    Y[y]=A[y]=a;
    }
}

Float
Layer::f1(Float a)
{
return 1.0;
}







Float
Layer::f_cost(Float* t)
{
Float sum=0.0;

for (int y=0; y<NY; y++)
	sum += (t[y]-Y[y])*(t[y]-Y[y]);
return sum;
}






Float
Layer::log_cost(Float* t)
{
Float sum=0.0;

for (int y=0; y<NY; y++) {
   if ((t[y]) && (Y[y]))
	sum += t[y]*(Float)log(Y[y]);
  }
return sum;
}




Float
Layer::sq_cost(Float* t)
{
Float sum=0.0;

for (int y=0; y<NY; y++) {
	sum += (t[y]-Y[y])*(t[y]-Y[y]);
  }
return sum;
}




void
Layer::gradient(int* I)
{
int y,u;

for (y=0; y<NY; y++) {
  for (u=0; u<NU; u++) {
    if (I[u]>=0) {
      dW[y][u][I[u]] += delta[y];
      }
    }
  dB[y] += delta[y];
  }
}




void
Layer::gradient(Float* I)
{
int y,u,k;
int i;

for (y=0; y<NY; y++) {
  i=0;
  for (u=0; u<NU; u++) {
    for (k=0; k<NK[u]; k++) {
      	dW[y][u][k] += delta[y]*I[i++];
	}
      }
  dB[y] += delta[y];
  }
}


void
Layer::gradient(int* I1,Float* I2)
{
int y,u,k;
int i;

for (y=0; y<NY; y++) {
  for (u=0; u<NU; u++) {
    if (I1[u]>=0) {
      dW[y][u][I1[u]] += delta[y];
      }
    }
  i=0;
  for (u=NU; u<NU+NUr; u++) {
    for (k=0; k<NK[u]; k++) {
      	dW[y][u][k] += delta[y]*I2[i++];
	}
      }
  dB[y] += delta[y];
  }
}


void
Layer::gradient(Float* I1,Float* I2)
{
int y,u,k;
int i1,i2;

for (y=0; y<NY; y++) {
  i1=0;
  for (u=0; u<NU; u++) {
    for (k=0; k<NK[u]; k++) {
      dW[y][u][k] += delta[y]*I1[i1++];
      }
    }
  i2=0;
  for (u=NU; u<NU+NUr; u++) {
    for (k=0; k<NK[u]; k++) {
      	dW[y][u][k] += delta[y]*I2[i2++];
	}
      }
  dB[y] += delta[y];
  }
}



void
Layer::gradient()
{
gradient(U);
}



void
Layer::resetGradient()
{
int y,u;
for (y=0; y<NY; y++) {
  for (u=0; u<NU+NUr; u++) {
    memset(dW[y][u], 0, NK[u]*sizeof(Float));
    }
  dB[y] = 0.0;
  }
}



void
Layer::initWeights(int seed)
{
int y,u,k;
Float D=(Float)(NU+NUr);
if (D>20) D=20.0;

//srand48(seed);
srand(seed);
for (y=0; y<NY; y++) {
  for (u=0; u<NU+NUr; u++) {
    for (k=0; k<NK[u]; k++) {
      W[y][u][k] = (Float)(0.5-(Float)rand()/(RAND_MAX+1))/D;
	if (rand() % 2 == 0)
	{
		W[y][u][k] = (0 - W[y][u][k]); 
	}
      }
    }
  B[y] = (Float)(0.5-(Float)rand()/(RAND_MAX+1))/D;
	if (rand() % 2 == 0)
	{
		B[y] = (0 - B[y]); 
	}
  }
}



void
Layer::set_dW(Float*** newdW) {
for (int y=0; y<NY; y++) {
  for (int u=0; u<NU+NUr; u++) {
    for (int k=0; k<NK[u]; k++) {
	dW[y][u][k]=newdW[y][u][k];
      }
    }
  }
}



// Layer_soft


Float
Layer_soft::f1(Float a)
{
int qua = -1;
Float sum=0;
int overflow=0;

for (int y=0;y<NY;y++) {
	if (A[y]==a) qua=y;
	if (A[y]>80) {
		overflow=1;
		}
	else
		sum += (Float)exp(A[y]);
	}
if (qua == -1)
	cerr << " AMBIGUITY!! ";
Float parz;
if (overflow) {
	parz=0.0;
	}
else
	parz = (Float)exp(a)/sum;

return (parz - parz*parz);
}


Float
Layer_soft::f_cost(Float* t)
{
return Layer::log_cost(t);
}



// Layer_tanh


Float
Layer_tanh::f1(Float a)
{
return 1.0-(Float)(tanh(a)*tanh(a));
}




Float
Layer_tanh::f_cost(Float* t)
{
return Layer::sq_cost(t);
}



