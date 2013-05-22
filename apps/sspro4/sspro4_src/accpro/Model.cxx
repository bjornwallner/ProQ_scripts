



#include "Model.h"




void
Model::alloc() {
int t;

FF = new Float*[MAX_T];
FFbp = new Float*[MAX_T];
BB = new Float*[MAX_T];
BBbp = new Float*[MAX_T];
P_F = new Float[NF];
P_B = new Float[NB];

for (t=0;t<MAX_T;t++) {
	FF[t] = new Float[NF];
	FFbp[t] = new Float[NF];
	BB[t] = new Float[NB];
	BBbp[t] = new Float[NB];
	}

Y=new Float[MAX_T*classes];	// here classes instead of the usual NY

Conf=new int*[NY];
for (int y=0;y<NY;y++)
	Conf[y]=new int[NY];
nerrors=new int[classes];
nerrors_0=new int[classes];
nerrors_1=new int[classes];
nall_0=new int[classes];
nall_1=new int[classes];



for (int f=0;f<NF;f++)
	P_F[f]=0;
for (int b=0;b<NB;b++)
	P_B[b]=0;
}



Model::Model(int the_NU, int the_NY, int the_NH, int the_context, int the_Moore, int the_NF,
	int the_NB, int the_NH2, int the_CoF, int the_CoB, int the_Step, Float the_threshold,
	int the_classes) :
NU(the_NU), NY(the_NY), NH(the_NH), context(the_context), Moore(the_Moore), NF(the_NF), 
NB(the_NB), NH2(the_NH2), CoF(the_CoF), CoB(the_CoB), Step(the_Step), 
threshold(the_threshold), classes(the_classes)
{

int NK[1024];

for (int c=0;c<2*context+1;c++) {
	NK[c]=NU;
	}
if (Moore)
	NetOut = new NNt((2*context+1), (2*CoF+1)*NF+(2*CoB+1)*NB, NH, classes, NK);
else
	NetOut = new NNt(0, (2*CoF+1)*NF+(2*CoB+1)*NB, NH, classes, NK);
NetF = new NNt((2*context+1), NF, NH2, NF, NK, 0);
NetB = new NNt((2*context+1), NB, NH2, NB, NK, 0);

NetOut->resetGradient();
NetF->resetGradient();
NetB->resetGradient();


alloc();
}





Model::Model(istream& is) {
is >> NU >> NY >> NH >> context >> classes;
is >> NF >> NB >> NH2 >> CoF >> CoB >> Step;

NetOut = new NNt(is);
NetF = new NNt(is, 0);
NetB = new NNt(is, 0);

Moore=NetOut->get_NI();
if (Moore) Moore=1;

threshold=0.0;

NetOut->resetGradient();
NetF->resetGradient();
NetB->resetGradient();


alloc();
}





void
Model::read(istream& is) {
is >> NU >> NY >> NH >> context >> classes;
is >> NF >> NB >> NH2 >> CoF >> CoB >> Step;

NetOut->read(is);
NetF->read(is, 0);
NetB->read(is, 0);

Moore=NetOut->get_NI();
if (Moore) Moore=1;

NetOut->resetGradient();
NetF->resetGradient();
NetB->resetGradient();

}




void
Model::write(ostream& os) {
os << NU << " " << NY << " " << NH << " " << context << " " << classes << "\n";
os <<NF<<" "<<NB<<" "<<NH2<<" "<<CoF<<" "<<CoB<<" "<<Step<<"\n";

NetOut->write(os);
NetF->write(os);
NetB->write(os);
}



void
Model::randomize(int seed) {

NetOut->initWeights(seed);
NetF->initWeights(seed++);
NetB->initWeights(seed++);

}



void
Model::F1_F(Sequence *seq, int t) {

seq->generate_profile();
Float I[512];
  for (int c=-context; c<=context; c++) {
    if (t+c <= 0 || t+c > seq->length) {
	for (int i=0;i<20;i++)
		I[20*(context+c)+i] = 0.0;
      } 
    else {
	for (int i=0;i<20;i++)
		I[20*(context+c)+i] = seq->HeP[20*(t+c-1)+i];
      }
    }

NetF->forward(I,FF[t-1]);
for (int f=0; f<NF; f++)
	FF[t][f] = NetF->out()[f];

}


void
Model::B1_B(Sequence* seq, int t) {

seq->generate_profile();
Float I[512];
  for (int c=-context; c<=context; c++) {
    if (t+c <= 0 || t+c > seq->length) {
	for (int i=0;i<20;i++)
		I[20*(context+c)+i] = 0.0;
      } 
    else {
	for (int i=0;i<20;i++)
		I[20*(context+c)+i] = seq->HeP[20*(t+c-1)+i];
      }
    }

NetB->forward(I,BB[t+1]);
for (int b=0; b<NB; b++)
	BB[t][b] = NetB->out()[b];

}


void
Model::propagate(Sequence* seq) {
int T=seq->length;
int t,f,b;

for (f=0; f<NF; f++)
	FF[0][f] = P_F[f];
for (b=0; b<NB; b++)
	BB[T+1][b] = P_B[b];

for (t=1; t<=T; t++)
	F1_F(seq,t);
for (t=T; t>0; t--)
	B1_B(seq,t);
}

void
Model::propagate(Sequence* seq, int T, int W) {
int T0=T-W;
int T1=T+W;
if (T0<1)
	T0=1;
if (T1>seq->length) 
	T1=seq->length;
int t,f,b;

for (f=0; f<NF; f++)
	FF[T0-1][f] = P_F[f];
for (b=0; b<NB; b++)
	BB[T1+1][b] = P_B[b];

for (t=T0; t<=T1; t++)
	F1_F(seq,t);
for (t=T1; t>=T0; t--)
	B1_B(seq,t);
}




void
Model::forward(Sequence* seq, int t) {
int f,b,v;
Float I[512];

// Next cycle sets the input vector
if (Moore) {
seq->generate_profile();
  for (int c=-context; c<=context; c++) {
    if (t+c <= 0 || t+c > seq->length) {
	for (int i=0;i<20;i++)
		I[20*(context+c)+i] = 0.0;
      } 
    else {
	for (int i=0;i<20;i++)
		I[20*(context+c)+i] = seq->HeP[20*(t+c-1)+i];
      }
    }
  }

  Float* X=new Float[(2*CoF+1)*NF+(2*CoB+1)*NB];
// We now set the hidden inputs
for (v=-CoF;v<=CoF;v++) {
  if ((t+(Step*v))<0 || (t+(Step*v))>seq->length) {
    for (f=0;f<NF;f++)
	X[NF*(CoF+v)+f]=0;
    }
  else {
    for (f=0;f<NF;f++)
	X[NF*(CoF+v)+f]=FF[t+(Step*v)][f];
    }
  }

for (v=-CoB;v<=CoB;v++) {
  if ((t+(Step*v))<1 || (t+(Step*v))>seq->length+1) {
    for (b=0;b<NB;b++)
	X[NF*(2*CoF+1) + NB*(CoB+v)+b]=0;
    }
  else {
    for (b=0;b<NB;b++)
	X[NF*(2*CoF+1) + NB*(CoB+v)+b]=BB[t+(Step*v)][b];
    }
  }

  NetOut->forward(I,X);

  for (int y=0;y<classes;y++)
	Y[classes*(t-1)+y] = NetOut->out()[y];
delete[] X;
}



void
Model::Feed(Sequence* seq) {

int t;

propagate(seq);

for (t=1; t<=seq->length; t++) {
	forward(seq,t);
	}

}



void
Model::predict(Sequence* seq) {

int t;
Float stepc=(Float)1/(Float)classes;

Feed(seq);

for (t=1; t<=seq->length; t++) {

//  if (seq->y[t]<0) 
//    {
//    seq->y_pred[t]=-1;
//    continue;
//    }


  
	int* predc=new int[classes]; 
	int* target=new int[classes]; 
	Float pred=0.0;
	int arg=-1;
	int c;

  for (c=0; c<classes; c++) {

	  if (Y[classes*(t-1)+c] > 0.5) {
		predc[c]=1;
	  } else {
		predc[c]=0;
	  }
	  if (seq->racc[t] <= stepc*c) {
		  target[c]=0;
		  nall_0[c] ++;
	  } else {
		  target[c]=1;
		  nall_1[c] ++;
	  }

	  if (target[c] != predc[c]) {
		  nerrors[c] ++;
		  if (target[c] == 0) {
			  nerrors_0[c] ++;
		  } else {
			  nerrors_1[c] ++;
		  }
	  }
  }

  seq->y_pred[t]=arg;

  delete[] target;
  delete[] predc;
}

}








void
Model::predict(Sequence* seq, int W) {

int t;
Float stepc=(Float)1/(Float)classes;

Feed(seq);

for (t=1; t<=seq->length; t++) {

int T0=1;
if ((t-W)>T0) T0=t-W;
int T1=seq->length;
if ((t+W)<T1) T1=t+W;

for (int me=T0-1;me<=T1+1;me++) {
	memset(FF[me],0,NF*sizeof(Float));
	memset(BB[me],0,NB*sizeof(Float));
	}

  if (seq->y[t]<0) 
    {
    seq->y_pred[t]=-1;
    continue;
    }
  propagate(seq, t, W);
  forward(seq, t);


  
	int* predc=new int[classes]; 
	int* target=new int[classes]; 
	Float pred=0.0;
	int arg=-1;
	int c;

  for (c=0; c<classes; c++) {

	  if (Y[classes*(t-1)+c] > 0.5) {
		predc[c]=1;
	  } else {
		predc[c]=0;
	  }
	  if (seq->racc[t] <= stepc*c) {
		  target[c]=0;
		  nall_0[c] ++;
	  } else {
		  target[c]=1;
		  nall_1[c] ++;
	  }

	  if (target[c] != predc[c]) {
		  nerrors[c] ++;
		  if (target[c] == 0) {
			  nerrors_0[c] ++;
		  } else {
			  nerrors_1[c] ++;
		  }
	  }
  }

  seq->y_pred[t]=arg;

  delete[] target;
  delete[] predc;
}

}
