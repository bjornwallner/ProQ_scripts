#include "Model.h"
#include <cmath>
#include <cstdlib>

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

Y=new Float[8196];

Conf=new int*[NY];
for (int y=0;y<NY;y++)
	Conf[y]=new int[NY];

for (int f=0;f<NF;f++)
	P_F[f]=0;
for (int b=0;b<NB;b++)
	P_B[b]=0;
}



Model::Model(int the_NU, int the_NY, int the_NH, int the_context, int the_Moore, int the_NF,
	int the_NB, int the_NH2, int the_CoF, int the_CoB, int the_Step, Float the_threshold, int use3d, int useacc, int binNum) :
NU(the_NU), NY(the_NY), NH(the_NH), context(the_context), Moore(the_Moore), NF(the_NF), NB(the_NB),
NH2(the_NH2), CoF(the_CoF), CoB(the_CoB), Step(the_Step), threshold(the_threshold), m_use3d(use3d), m_useacc(useacc), m_binNum(binNum)
{

m_full3d = true;  //true: forward/backward use 3d info. false: forward/backward not use 3d info. , need to save to model, thus decide how to read weights. 

if (m_use3d != 0 && the_NU <= 20) 
{
   cout << "if 3d information is used, the input must be bigger than 20.\n"; 
   exit(1); 
}

int NK[1024];

for (int c=0;c<2*context+1;c++) {
		if (m_useacc == 0)
		{
			//don't use acc information
			NK[c]=NU;
		}
		else
		{
			//one extra input for solv acc
			NK[c] = NU + 1;  
		}
	}

if (Moore)
	NetOut = new NN((2*context+1), (2*CoF+1)*NF+(2*CoB+1)*NB, NH, NY, NK);
else
	NetOut = new NN(0, (2*CoF+1)*NF+(2*CoB+1)*NB, NH, NY, NK);

if (m_use3d == 0)
{
	NetF = new NNt((2*context+1), NF, NH2, NF, NK, 0);
	NetB = new NNt((2*context+1), NB, NH2, NB, NK, 0);
}
else
{
	for (int c=0;c<2*context+1;c++) {
		if (m_full3d)
		{
			if (m_useacc == 0)
			{
				NK[c] = NU;
			}
			else
			{
				//one extra input for solvent acc
				NK[c] = NU + 1; 
			}
		}
		else
		{
			if (m_useacc = 0)
			{
				NK[c]=20;
			}
			else
			{
				NK[c] = 21; //one input for solvent acc
			}
		}

	}

	NetF = new NNt((2*context+1), NF, NH2, NF, NK, 0);
	NetB = new NNt((2*context+1), NB, NH2, NB, NK, 0);
}

alloc();
}





Model::Model(istream& is) {
//temporary comment for non-full 3d
m_full3d = true;

//m_full3d = false;

m_useacc = 0; 
is >> NU >> NY >> NH >> context;

//read in whether or not 3d information is used in training. 
//not compatible with model trained use old sspro. 
is >> m_use3d ;
is >> m_useacc; 
is >> m_binNum; 

is >> NF >> NB >> NH2 >> CoF >> CoB >> Step;

NetOut = new NN(is);
NetF = new NNt(is, 0);
NetB = new NNt(is, 0);

Moore=NetOut->get_NI();
if (Moore) Moore=1;
//cout << "Moore: " << Moore << endl; 

//fix bug that program can't learn when read from a model
threshold=0;

alloc();
}





void
Model::read(istream& is) {
m_full3d = true;
m_useacc = 0; 
is >> NU >> NY >> NH >> context;

//read in whether or not 3d information is used in training. 
//not compatible with model trained use old sspro. 
is >> m_use3d ;
is >> m_useacc; 
is >> m_binNum; 

is >> NF >> NB >> NH2 >> CoF >> CoB >> Step;

NetOut->read(is);
NetF->read(is, 0);
NetB->read(is, 0);

Moore=NetOut->get_NI();
if (Moore) Moore=1;
//cout << "Moore: " << Moore << endl; 
threshold = 0; 

}




void
Model::write(ostream& os) {

//add: output whether or not 3d inforamtion is used. 
//add: output whether or not acc information is used. 
//if want to use old models, must add this number in model file manually. 
os << NU << " " << NY << " " << NH << " " << context << " " << m_use3d << " " << m_useacc << " " << m_binNum << "\n";

os <<NF<<" "<<NB<<" "<<NH2<<" "<<CoF<<" "<<CoB<<" "<<Step<<"\n";

NetOut->write(os);
NetF->write(os);
NetB->write(os);
}



Float * Model::prepareInput(Sequence* seq, int t)
{

  if (m_full3d && m_use3d != 0)
  {
	return prepareInputWith3d(seq, t); 
  }
  int profileLength = 20;
  int inputLength = 20;
  if (m_useacc != 0)
  {
	//one extra input for solv acc
	inputLength = 21; 
  }

  seq->generate_profile();
  for (int c=-context; c<=context; c++) {
    if (t+c <= 0 || t+c > seq->length) {
	for (int i=0;i<inputLength;i++)
		m_input[inputLength*(context+c)+i] = 0.0;
      } 
    else {
	//copy profile from sequence
	for (int i=0;i<profileLength;i++)
		m_input[inputLength*(context+c)+i] = seq->HeP[profileLength*(t+c)+i];
	//copy solve acc info if it need to be used. 
	
	if (m_useacc != 0)
	{
  		if (seq->accessibility[t] > 25)
			m_input[inputLength*(context+c)+profileLength] = 1;
  		else
			m_input[inputLength*(context+c)+profileLength] = 0;
	}

      }
    }
    return m_input;
}

Float Model::calcDistance(Float x1, Float y1, Float z1, Float x2, Float y2, Float z2)
{
	double total = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2); 
	return sqrt(total); 
}

Float * Model::prepareInputWith3d(Sequence * seq, int t)
{
	if (m_use3d == 0)
	{
		cout << "model doesn't use 3d, can't create 3d input, program exit. \n"; 

		exit(1); 
	}

	if (m_binNum != 1 && m_binNum != 2 && m_binNum != 3)
	{
		cout << " if use 3d information, bin number should be 2 or 3.\n"; 
		exit(1); 
	}
	//for testing if normalization of count helps????
	bool bNormalized = true; 

	//currently lenght is 89: 0-19 for profile of current aa, 20-42 for neighbor in 8A, 43-65 in 12A, 66-88 in 16A

	//lenght of input for each range of neighbors
	int extraLength = 23; 

	//standard length of profile
	const int profileLength = 20; 
	
	// lenght of profile + one solvent acc if solv acc is used. 
	int stepLength = 20; 
	if (m_useacc != 0)
	{
		//one extra input for solvent accessibility for each amino acids
		stepLength = 21; 
	}

	//number of bins
	const int binNum = m_binNum; 

	//distance threshold of each bin
	int threshold[3] = {8, 12, 16}; 
	if (binNum == 2)
	{
		threshold[0] = 12;
		threshold[1] = 16;
	}
	else if (binNum == 3)
	{
		threshold[0] = 8;
		threshold[1] = 12;
		threshold[2] = 16; 
	}
	else if (binNum == 1)
	{
		threshold[0] = 12; 
	}

	//numerical seperation threshold
	const int dist = 5; 
	//specify a range of seperation
	const int minSeperation = 5;
	const int maxSeperation = 10000; 

	int helixNum[3];
	int sheetNum[3];
	int coilNum[3];

	int totalLength = stepLength + extraLength * binNum; 

	if (m_useacc == 0 && totalLength != NU)
	{
		cout << " in preparing 3d input, total length doesn't match with NU.\n"; 
		exit(1); 
	}
	if (m_useacc == 1 && totalLength != NU + 1)
	{
		cout << " in preparing 3d input, total length doesn't match with NU.\n"; 
		exit(1); 
	}

	memset(m_3dinput, 0, sizeof(Float) * INPUT_SIZE);  

	//generate profile of current amino acid. 
	seq->generate_profile();


	for (int c=-context; c<=context; c++) 
	{
		int index = t + c;  //start from 1

		if (index <= 0 || index > seq->length) 
		{
			continue; 
		} 

		//generate profile for amino acide at position: index
		for (int i=0; i<profileLength; i++)
			m_3dinput[totalLength*(context+c)+i] = seq->HeP[profileLength * index + i];

		//set solvent acc if necessary 
		if (m_useacc != 0)
		{
  			if (seq->accessibility[index] > 25)
				m_3dinput[totalLength*(context+c)+profileLength] = 1;
  			else
				m_3dinput[totalLength*(context+c)+profileLength] = 0;
		}

		//generate profile for neighbors according to thresholds
		for (int m = 0; m < binNum; m++)
		{
			helixNum[m] = 0;
			sheetNum[m] = 0;
			coilNum[m] = 0; 
		}

		for (int i = 1; i <= seq->length; i++)
		{
		//	if ( abs(i - index) <= dist)
		//	{
		//		continue; 	
		//	}
			if (abs(i - index) <= minSeperation || abs(i -index) >= maxSeperation)
			{
				continue; 
			}

	//		double distance = calcDistance(seq->getX(index), seq->getY(index), seq->getZ(index), 
	//					seq->getX(i), seq->getY(i),seq->getZ(i) ); 
			
			//reengineering: add isContact(int index1, int index2, float threshold) to Sequence
			//thus to encapsulate the changes. 

				
			for (int j = 0; j < binNum; j++)
			{
	//			if (distance < threshold[j])
				if ( seq->isContacted(index, i, threshold[j]) )
				{
					if (seq->getSS(i) == 0)
					{
						helixNum[j]++;
					}
					else if (seq->getSS(i) == 1)
					{
						sheetNum[j]++; 
					}
					else
					{
						coilNum[j]++; 
					}

					for (int k=0;k < profileLength; k++)
						m_3dinput[totalLength*(context+c)+stepLength + extraLength * j + k ] = 
							seq->HeP[profileLength * i +k];
				}
			}
		
		}
		for (int j = 0; j < binNum; j++)
		{
			float totalNum = helixNum[j] + sheetNum[j] + coilNum[j];
			if (totalNum > 0)
			{	
				for (int k=0;k<profileLength;k++)
					m_3dinput[totalLength*(context+c)+ stepLength + extraLength * j + k ] /= totalNum;
			}
			
			if (!bNormalized) //testing code to see if normalized help
			{
			m_3dinput[totalLength*(context+c)+stepLength + extraLength * j + profileLength ] = helixNum[j] ;
			m_3dinput[totalLength*(context+c)+stepLength + extraLength * j + profileLength + 1] = sheetNum[j] ;
			m_3dinput[totalLength*(context+c)+stepLength + extraLength * j + profileLength + 2] = coilNum[j] ;
			}
			else
			{
				if (totalNum > 0)
				{
					m_3dinput[totalLength*(context+c)+stepLength + extraLength * j + profileLength ] = helixNum[j]/totalNum ;
					m_3dinput[totalLength*(context+c)+stepLength + extraLength * j + profileLength + 1] = sheetNum[j]/totalNum ;
					m_3dinput[totalLength*(context+c)+stepLength + extraLength * j + profileLength + 2] = coilNum[j]/totalNum ;
				}
			}
		}


    	}
	return m_3dinput; 
}

void
Model::F1_F(Sequence *seq, int t) {
/*
seq->generate_profile();
Float I[512];
  for (int c=-context; c<=context; c++) {
    if (t+c <= 0 || t+c > seq->length) {
	for (int i=0;i<20;i++)
		I[20*(context+c)+i] = 0.0;
      } 
    else {
	for (int i=0;i<20;i++)
		I[20*(context+c)+i] = seq->HeP[20*(t+c)+i];
      }
    }
*/
Float * I = prepareInput(seq, t); 
NetF->forward(I,FF[t-1]);
for (int f=0; f<NF; f++)
	FF[t][f] = NetF->out()[f];

}



void
Model::B1_B(Sequence* seq, int t) {

Float * I = prepareInput(seq, t); 

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
//Float I[512];
Float * I = NULL; 
Float temp[INPUT_SIZE]; 

// Next cycle sets the input vector
if (Moore) {
/*
seq->generate_profile();
  for (int c=-context; c<=context; c++) {
    if (t+c <= 0 || t+c > seq->length) {
	for (int i=0;i<20;i++)
		I[20*(context+c)+i] = 0.0;
      } 
    else {
	for (int i=0;i<20;i++)
		I[20*(context+c)+i] = seq->HeP[20*(t+c)+i];
      }
    }
*/ 
	if (m_use3d == 0)
	{
		I = prepareInput(seq, t); 
	}
	else
	{
   		 I = prepareInputWith3d(seq, t); 
	}
  }
  else //hack: to maintain previous behavior. 
  {
	I = temp; 
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

  for (int y=0;y<3;y++)
	Y[3*(t-1)+y] = NetOut->out()[y];
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

Feed(seq);

for (t=1; t<=seq->length; t++) {

  if (seq->y[t]<0) 
    {
    seq->y_pred[t]=-1;
    continue;
    }

  Float pred=0.0;
  int arg=-1;

  seq->alpha[t] = Y[3*(t-1)];
  seq->beta[t]  = Y[3*(t-1)+1];
  seq->gamma[t] = Y[3*(t-1)+2];

  for (int c=0; c<3; c++) 
    {
    if (Y[3*(t-1)+c] > pred) 
      {
      pred = Y[3*(t-1)+c];
      arg=c;
      }
    }
  seq->y_pred[t]=arg;

  }


for (t=1; t<=seq->length; t++) {
  Conf[seq->y_pred[t]][seq->y[t]]++;
  if (seq->y[t]!=seq->y_pred[t]) {
	nerrors++;
	if (seq->y[t]==0) nerrors_alpha++;
	if (seq->y[t]==1) nerrors_beta++;
	if (seq->y[t]==2) nerrors_gamma++;
	}
  }


}


void
Model::predict(Sequence* seq, int W) {

int t,c;

for (t=1; t<=seq->length; t++) {

int T0=1;
if ((t-W)>T0) T0=t-W;
int T1=seq->length;
if ((t+W)<T1) T1=t+W;

for (c=T0-1;c<=T1+1;c++) {
	memset(FF[c],0,NF*sizeof(Float));
	memset(BB[c],0,NB*sizeof(Float));
	}

  if (seq->y[t]<0) 
    {
    seq->y_pred[t]=-1;
    continue;
    }
  propagate(seq, t, W);
  forward(seq, t);

  Float pred=0.0;
  int arg=-1;

  for (c=0; c<3; c++) 
    {
    if (Y[3*(t-1)+c] > pred) 
      {
      pred = Y[3*(t-1)+c];
      arg=c;
      }
    }
  seq->y_pred[t]=arg;

  }


for (t=1; t<=seq->length; t++) {
  if (seq->y[t]!=seq->y_pred[t]) {
  Conf[seq->y_pred[t]][seq->y[t]]++;
	nerrors++;
	if (seq->y[t]==0) nerrors_alpha++;
	if (seq->y[t]==1) nerrors_beta++;
	if (seq->y[t]==2) nerrors_gamma++;
	}
  }


}



