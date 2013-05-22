


#ifndef Sequence_h
#define Sequence_h 1

#include <iostream.h>
#include "General.h"




static int translateU[26] = {
  0,	// A
  -1,	// B
  1,	// C
  2,	// D
  3,	// E
  4,	// F
  5,	// G
  6,	// H
  7,	// I
  -1,	// J
  8,	// K
  9,	// L
  10,	// M
  11,	// N
  -1,	// O
  12,	// P
  13,	// Q
  14,	// R
  15,	// S
  16,	// T
  -1,	// U
  17,	// V
  18,	// W
  0,	// X
  19,	// Y
  -1	// Z
};


static char Utranslate[26] = "ACDEFGHIKLMNPQRSTVWY";

// Computed on SWISS-PROT v.38
static Float frequencies[21] = {
(Float)0.0758211047836989,
(Float)0.016617224149173,
(Float)0.0527726688799907,
(Float)0.063666101502907,
(Float)0.0410241159267021,
(Float)0.0684499551587853,
(Float)0.0224934259530327,
(Float)0.058104312509487,
(Float)0.0594638685702881,
(Float)0.0943875852150685,
(Float)0.0237836702340803,
(Float)0.0444472101922697,
(Float)0.0492088194426418,
(Float)0.0397119366677365,
(Float)0.0516215294902541,
(Float)0.0713318261917733,
(Float)0.0567592995453305,
(Float)0.0658317851926178,
(Float)0.0124027172555561,
(Float)0.0319064538515397,
(Float)0.000194389287066804
};




#define MAX_T 8192

class Sequence {
public:
  char ali[512];
  char alidir[512];
  int length;
  int* u;
  int* y;
  int* y_pred;

  Float* alpha;
  Float* beta;
  Float* gamma;

  int* contacts;
  int* avg_contacts;

  int* ss;
  int* acc;
  Float* racc;

  Float* entropy;
  Float* present;
  int* settato;

  int alignments_loaded;
  int** ALIGNMENTS;
  int* SCORES;
  int* INTERS;

  Float* HeP;
  int HePl;

  Float belief;

  Sequence()
    {
    belief=0.0;
  };
  Sequence(istream& is, char* the_alidir)
    {
	  strcpy(alidir,the_alidir);
	  belief=0.0;
      int t;
      is >> length;
      if (length >= MAX_T-1) {
	FAULT("Sequence too long (" << length << ") MAX_T is " << MAX_T);
      }
      u = new int[length+1];
      y = new int[length+1];
      y_pred = new int[length+1];

      alpha = new Float[length+1];
      beta = new Float[length+1];
      gamma = new Float[length+1];

      entropy = new Float[length+1];
      present = new Float[length+1];

      settato = new int[length+1];

      u[0]=y[0]=0;

	  u[0]=0;

	for (t=1; t<=length; t++) {
		is >> u[t];
		is >> y[t];
		settato[t]=0;
      }
    };

  ~Sequence() {
      delete[] u;
      delete[] y;
      delete[] y_pred;
      delete[] alpha;
      delete[] beta;
      delete[] gamma;

      delete[] entropy;
      delete[] present;
      delete[] settato;
  }

  virtual void set_belief(Float be) {};
  virtual int read_u(istream& is)=0;

  virtual int load_alignments() {return 0;};
  virtual int unload_alignments() {return 0;};
  virtual int Aps(int al,int offs) {return 0;};
  virtual int best_alignments() {return 0;};
  virtual int best_alignments(int AL) {return 0;};
  virtual Float profile_entropy() {return 0;}
  virtual int generate_profile() {return 0;}
  virtual void unload_profile() {}
  virtual void write(ostream& os)
    {
      int t;
      os << length << "\n";
      for (t=1; t<=length; t++) {
	os << u[t] << " ";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	os << y[t] << " ";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	os << y_pred[t] << " ";
      }
      os << "\n";
    };
  virtual void write_u(ostream& os)=0;

};








class SymbolicSequence: public Sequence {
public:
  int off;
  char des[512];
  char des2[512];
  char des3[512];


  SymbolicSequence() : Sequence() {}

  SymbolicSequence(istream& is, char* the_alidir) : Sequence()
    {
	  strcpy(alidir,the_alidir);
	  belief=0.0;
      HePl = 0;
      char tmp[MAX_T];
      int t;
      is >> des;

//	  cerr << des << "\n";

// PSI next 2 lines
//	  is >> des2;
//	  is >> des2;

      is >> tmp;

// PSI next line instead of the following two
//      strcpy(ali,des);
      strcpy(ali,"");
      ext_ali();

      alignments_loaded=0;

      length = strlen(tmp);
      if (length >= MAX_T-1) {
	FAULT("Sequence too long (" << length << ") MAX_T is " << MAX_T);
      }
      u = new int[length+1];
      y = new int[length+1];
      y_pred = new int[length+1];
      alpha = new Float[length+1];
      beta = new Float[length+1];
      gamma = new Float[length+1];

      entropy = new Float[length+1];
      present = new Float[length+1];
      settato = new int[length+1];

      u[0]=y[0]=0;
      for (t=1; t<=length; t++) {
	u[t]=translateU[tmp[t-1]-'A'];
//	   if (u[t]<0 || u[t]>19) {u[t]=-1;}
      }
      is >> tmp;
      if (strlen(tmp)!=length) {
		  cerr << tmp << " " << strlen(tmp) << " " << length << "\n";
		  FAULT("I/O length mismatch");
	  }
      for (t=1; t<=length; t++) {
	settato[t]=0;
	if (tmp[t-1]=='-')
	  y[t]=-1;
	else if (tmp[t-1]=='H')
	  y[t]=0;
	else if (tmp[t-1]=='G')
	  y[t]=1;
	else if (tmp[t-1]=='I')
	  y[t]=2;
	else if (tmp[t-1]=='E')
	  y[t]=3;
	else if (tmp[t-1]=='B')
	  y[t]=4;
	else if (tmp[t-1]=='S')
	  y[t]=5;
	else if (tmp[t-1]=='T')
	  y[t]=6;
	else if (tmp[t-1]=='.')
	  y[t]=7;
      }
    };

  SymbolicSequence(int l) : Sequence()
  {
	  belief=0.0;
      HePl = 0;
      alignments_loaded=0;

      length = l;

      u = new int[length+1];
	  memset(u,0,(length+1)*sizeof(int));
      y = new int[length+1];
      y_pred = new int[length+1];

      alpha = new Float[length+1];
      beta = new Float[length+1];
      gamma = new Float[length+1];

      entropy = new Float[length+1];
      present = new Float[length+1];
      settato = new int[length+1];

      u[0]=y[0]=0;

  }

  ~SymbolicSequence() {
      delete[] u;
      delete[] y;
      delete[] y_pred;
      delete[] alpha;
      delete[] beta;
      delete[] gamma;

      delete[] entropy;
      delete[] present;
      delete[] settato;
  }

void set_belief(Float be) {
	belief=be;
	};

int read_u(istream& is)
    {
      HePl = 0;
      char tmp[MAX_T];
	  strcpy(tmp,"");
      int t;
      is >> tmp;
      length = strlen(tmp);
      if (length >= MAX_T-1) {
		FAULT("Sequence too long (" << length << ") MAX_T is " << MAX_T);
      }
      u[0]=0;
      for (t=1; t<=length; t++) {
		u[t]=translateU[tmp[t-1]-'A'];
      }
	return length;
    };


  int copy_alignment(Sequence* seq,int n, int off, int le)
  {
  int t;
  if (seq->ALIGNMENTS[n][0]-off <le) le = seq->ALIGNMENTS[n][0]-off;

  for (t=1;t<=le;t++)
	{
	u[t]=seq->ALIGNMENTS[n][t+off];
	y[t]=seq->y[t];
	settato[t]=0;
	}
  length = le;
  return 0;
  }


  int ext_ali()
   {
   char* pi=des;
   int cont=0;
   while (pi[0]!='.' && pi[0]!=' '  && cont<512)
	{
	cont++;
	strncat(ali,pi,1);
	pi++;
	}
   return 1;
   }



int load_alignments()
	{
  	filebuf inbuf;
	char temp[MAX_T];
  	char fname[1024];


	strcpy(fname, alidir);


	strcat(fname, ali);
//	cout << fname << "\n" << flush;

	if (inbuf.open(fname, ios::in) != 0) {
    	  istream is(&inbuf);
	  is >> alignments_loaded;
//	cout << alignments_loaded << "\n" << flush;
	  if (alignments_loaded) {
		  ALIGNMENTS=new int*[alignments_loaded];
		  SCORES=new int[alignments_loaded];
		  INTERS=new int[alignments_loaded];
	  }
	  for (int quanti=0;quanti<alignments_loaded;quanti++)
		{
		is>>temp;
      		int ali_length = strlen(temp);
		ALIGNMENTS[quanti]=new int[ali_length+1];
		ALIGNMENTS[quanti][0]=ali_length;
		for (int i=1;i<=ali_length;i++)
			{
			if (temp[i-1]=='.' || temp[i-1]=='-')
			  ALIGNMENTS[quanti][i]=-1;
			else
			  ALIGNMENTS[quanti][i]=translateU[temp[i-1]-'A'];
			}
		}
  	} else {
	#ifdef _TRAIN_INFO
	 cout << "alignent of " << fname << " not exists\n"; 
	#endif
    	 alignments_loaded=0;
  	}
  	inbuf.close();

	return alignments_loaded;
	}

int unload_alignments()
	{
	if (alignments_loaded) {
	  for (int quanti=0;quanti<alignments_loaded;quanti++)
		delete[] ALIGNMENTS[quanti];
	  delete[] ALIGNMENTS;
	  delete[] SCORES;
	  delete[] INTERS;
	  alignments_loaded = 0;
	  return 1;
	  }
	else return 0;
	}

int Aps(int al,int offs)
	{
	int somma=0;
	int minl=length;
	int lali=ALIGNMENTS[al][0];
	if (lali-offs<minl)
		minl=lali-offs;
	INTERS[al]=0;
	for (int c=1;c<=minl;c++) {
            if (ALIGNMENTS[al][c+offs] != -1)
		INTERS[al]++;
	    if (u[c]==ALIGNMENTS[al][c+offs])
		somma++;
	    }
	return somma;
	}


int best_alignments()
	{
	if (alignments_loaded==0) return -1;
//	int score=0;
	Float score=0.0;
	int offset=0;
	int right=length;
	if (ALIGNMENTS[0][0]>right)
		right=ALIGNMENTS[0][0];
	int start=0;

	for (int o=0;o<right;o++)
	  {
	  start=0;
	  int parz=0;
	  int al_found=0;
	  while	(start<alignments_loaded)
		{
//		parz += Aps(start,o);
//	  	start++;
		int apsp=Aps(start,o);
		if (INTERS[start]) al_found++;
		parz+=apsp*INTERS[start];
	  	start++;
		}

	  if (al_found==0) al_found=1;

//	  if (parz >score)
	  if ((Float)parz/(Float)al_found >score)
		{
//		score = parz;
		score = (Float)parz/(Float)al_found;
		offset=o;
		}
	  }
	for (int a=0;a<alignments_loaded;a++)
		SCORES[a]=Aps(a,offset);
	return offset;
	}

int best_alignments(int AL)
	{
	if (alignments_loaded<AL) return -1;
	int score=0;
	int offset=0;
	int right=length;
	if (ALIGNMENTS[AL][0]>right)
		right=ALIGNMENTS[AL][0];
	int start=0;

	for (int o=0;o<right;o++)
	  {
	  int parz=0;
	  parz = Aps(AL,o);
	  parz *= INTERS[AL];

	  if (parz >score)
		{
		score = parz;
		offset=o;
		}
	  }
	SCORES[AL]=Aps(AL,offset);
	return offset;
	}


#ifdef NEVER




int generate_profile() {

    if (HePl) return 1;

	Float* Pres= new Float[length+1];

	HeP = new Float[20*(length+1)];

	memset(HeP,0,20*(length+1)*sizeof(Float));
	memset(Pres,0,(length+1)*sizeof(Float));
	int t;

	for (t=1;t<=length;t++) {
		int aa = u[t];
		if (aa != -1) {
			HeP[20*(t-1)+aa]++;
			Pres[t-1]++;
			}
	   	}

      load_alignments();
      int offset=0;
      int AL= alignments_loaded;

	for (int a=0;a<AL;a++) {
	   offset = best_alignments(a);
	   int sum=Aps(a,offset);
	   Float qual = (Float)sum/(Float)INTERS[a];

	   int minl=length;
	   int lali=ALIGNMENTS[a][0];
	   if (lali-offset<minl)
		minl=lali-offset;

	   if ((qual > 0.2) && INTERS[a] >= 10) {

	     for (t=1;t<=minl;t++) {
		int aa = ALIGNMENTS[a][t+offset];
		if (aa != -1) {
			HeP[20*(t-1)+aa]++;
			Pres[t-1]++;
			}
	   	}
	     }
	   }

	for (t=1;t<length;t++) {
		Float total=0.0;
		int aa;
		for (aa=0;aa<20;aa++) {
			HeP[20*(t-1)+aa] += belief*frequencies[aa];
			total += HeP[20*(t-1)+aa];
		}
		for (aa=0;aa<20;aa++) {
			if (total)
				HeP[20*(t-1)+aa] /= total;
			}
		}

	unload_alignments();
	HePl=1;
	delete[] Pres;

	return 0;

}


#endif


Float profile_entropy() {
Float sum=0;
for (int t=0;t<length;t++) {
	for (int aa=0;aa<20;aa++) {
		if (HeP[20*t+aa])
			sum -= HeP[20*t+aa]*(Float)log(HeP[20*t+aa]);
		if (HeP[20*t+aa]<0 || HeP[20*t+aa]>1) {
			cout << HeP[20*t+aa] << " " << aa << " " << t <<" "<< length<<"\n" << flush;
			}
		}
	}
return sum;
}


int generate_profile() {

    if (HePl) return 1;

	Float* Pres= new Float[length+1];

	HeP = new Float[20*(length+1)];


	int* flags;
	int* offsets;
	int* minls;

	memset(HeP,0,20*(length+1)*sizeof(Float));
	memset(Pres,0,(length+1)*sizeof(Float));
	int t, a;

	for (t=1;t<=length;t++) {
		int aa = u[t];
		if (aa != -1) {
			HeP[20*(t-1)+aa]++;
			Pres[t-1]++;
			}
	   	}

      load_alignments();
      int offset=0;
      int AL= alignments_loaded;

	  if (AL) {
	flags=new int[AL];
	offsets=new int[AL];
	minls=new int[AL];
	memset(flags,0,AL*sizeof(int));
	memset(offsets,0,AL*sizeof(int));
	memset(minls,0,AL*sizeof(int));
	  }

	for (a=0;a<AL;a++) {
	   offset = best_alignments(a);
	   int sum=Aps(a,offset);
//	   Float qual = (Float)sum/(Float)INTERS[a];

	   int minl=length;
	   int lali=ALIGNMENTS[a][0];
	   if (lali-offset<minl)
		minl=lali-offset;

//	   if ((qual > 0.2) && INTERS[a] >= 10) {
	   if (1) {
		   flags[a]=1;
		   offsets[a]=offset;
		   minls[a]=minl;
		   for (t=1;t<=minl;t++) {
			   int aa = ALIGNMENTS[a][t+offset];
			   if (aa != -1) {
				   HeP[20*(t-1)+aa]++;
				   Pres[t-1]++;
			   }
		   }
	   }
	}

	for (t=1;t<=length;t++) {
		Float total=0.0;
		int aa;
		for (aa=0;aa<20;aa++) {
			HeP[20*(t-1)+aa] += belief*frequencies[aa];
			total += HeP[20*(t-1)+aa];
		}
		for (aa=0;aa<20;aa++) {
			if (total)
				HeP[20*(t-1)+aa] /= total;
		}
	}


Float* weights;
if (AL) {
	weights=new Float[AL];
}

//cout << profile_entropy() << "\n" << flush;

for (int iter=0;iter<1;iter++) {

// We can now think to weight each sequence based on its information content

	if (AL) {
		memset(weights, 0, AL*sizeof(Float));
	}

	for (a=0;a<AL;a++) {
		if (!flags[a]) {
			continue;
		}


		Float sum=0;

		for (t=1;t<=minls[a];t++) {
			int aa = ALIGNMENTS[a][t+offsets[a]];
			if (aa != -1) {
				sum -= (Float)log(HeP[20*(t-1)+aa]);
			}
		}
		sum /= (Float)INTERS[a];
		weights[a]=sum;
	}



// Now we recompute the profile matrix

	memset(HeP,0,20*(length+1)*sizeof(Float));

	for (t=1;t<=length;t++) {
		int aa = u[t];
		if (aa != -1) {
			HeP[20*(t-1)+aa]++;
		}
   	}

	for (a=0;a<AL;a++) {
		if (flags[a]) {
			for (t=1;t<=minls[a];t++) {
				int aa = ALIGNMENTS[a][t+offsets[a]];
				if (aa != -1) {
					HeP[20*(t-1)+aa] += weights[a];
				}
	   		}
		}
	}

	for (t=1;t<=length;t++) {
		Float total=0.0;
		int aa;
		for (aa=0;aa<20;aa++) {
			HeP[20*(t-1)+aa] += belief*frequencies[aa];
			total += HeP[20*(t-1)+aa];
		}
		for (aa=0;aa<20;aa++) {
			if (total)
				HeP[20*(t-1)+aa] /= total;
		}
	}
//cout << profile_entropy() << "\n" << flush;
}

//cout << "\n"<<flush;

unload_alignments();
HePl=1;

delete[] Pres;
if (AL) {
	delete[] flags;
	delete[] offsets;
	delete[] minls;
	delete[] weights;
}

return 0;

}













void unload_profile() {
	HePl = 0;
	delete[] HeP;
	}


  void write(ostream& os)
    {
      int t;
      os << des << "\n";
      for (t=1; t<=length; t++) {
		  if (u[t]>=0) {
			os << Utranslate[u[t]];
		  }
		  else os<<"X";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y[t]==-1) os << "-";
	else if (y[t]==0) os << "H";
	else if (y[t]==1) os << "G";
	else if (y[t]==2) os << "I";
	else if (y[t]==3) os << "E";
	else if (y[t]==4) os << "B";
	else if (y[t]==5) os << "S";
	else if (y[t]==6) os << "T";
	else if (y[t]==7) os << ".";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y_pred[t]==-1) os << "-";
	else if (y_pred[t]==0) os << "H";
	else if (y_pred[t]==1) os << "G";
	else if (y_pred[t]==2) os << "I";
	else if (y_pred[t]==3) os << "E";
	else if (y_pred[t]==4) os << "B";
	else if (y_pred[t]==5) os << "S";
	else if (y_pred[t]==6) os << "T";
	else if (y_pred[t]==7) os << ".";
      }
      os << "\n\n";
    };

  void write_u(ostream& os)
    {
      int t;
//      os << des << "\n";
      for (t=1; t<=length; t++) {
	  os << Utranslate[u[t]];
      }
      os << "\n";
    };


  void write_predictions(ostream& os)
    {
      int t;
      os << des << "\n";
      for (t=1; t<=length; t++) {
		  if (u[t]>=0) {
			os << Utranslate[u[t]];
		  }
		  else os<<"X";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y[t]==-1) os << "-";
	else if (y[t]==0) os << "H";
	else if (y[t]==1) os << "G";
	else if (y[t]==2) os << "I";
	else if (y[t]==3) os << "E";
	else if (y[t]==4) os << "B";
	else if (y[t]==5) os << "S";
	else if (y[t]==6) os << "T";
	else if (y[t]==7) os << ".";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y_pred[t]==-1) os << "-";
	else if (y_pred[t]==0) os << "H";
	else if (y_pred[t]==1) os << "G";
	else if (y_pred[t]==2) os << "I";
	else if (y_pred[t]==3) os << "E";
	else if (y_pred[t]==4) os << "B";
	else if (y_pred[t]==5) os << "S";
	else if (y_pred[t]==6) os << "T";
	else if (y_pred[t]==7) os << ".";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y_pred[t]==-1) os << "0 0 0";
	else os << alpha[t] << " " << beta[t] << " " << gamma[t];
	os << " ";
      }
      os << "\n\n";
    };




  void write3(ostream& os)
    {
      int t;
      os << des << "\n";
      for (t=1; t<=length; t++) {
		  if (u[t]>=0) {
			os << Utranslate[u[t]];
		  }
		  else os<<"X";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y[t]==-1) os << "-";
	else if (y[t]==0 || y[t]==1) os << "H";
	else if (y[t]==3 || y[t]==4) os << "E";
	else os << "C";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y_pred[t]==-1) os << "-";
	else if (y_pred[t]==0 || y_pred[t]==1) os << "H";
	else if (y_pred[t]==3 || y_pred[t]==4) os << "E";
	else os << "C";
      }
      os << "\n\n";
    };



  void write3_predictions(ostream& os)
    {
      int t;
      os << des << "\n";
      for (t=1; t<=length; t++) {
		  if (u[t]>=0) {
			os << Utranslate[u[t]];
		  }
		  else os<<"X";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y[t]==-1) os << "-";
	else if (y[t]==0) os << "H";
	else if (y[t]==1) os << "G";
	else if (y[t]==2) os << "I";
	else if (y[t]==3) os << "E";
	else if (y[t]==4) os << "B";
	else if (y[t]==5) os << "S";
	else if (y[t]==6) os << "T";
	else if (y[t]==7) os << ".";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y_pred[t]==-1) os << "-";
	else if (y_pred[t]==0) os << "H";
	else if (y_pred[t]==1) os << "G";
	else if (y_pred[t]==2) os << "I";
	else if (y_pred[t]==3) os << "E";
	else if (y_pred[t]==4) os << "B";
	else if (y_pred[t]==5) os << "S";
	else if (y_pred[t]==6) os << "T";
	else if (y_pred[t]==7) os << ".";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y_pred[t]==-1) os << "0 0 0";
	else os << alpha[t] << " " << beta[t] << " " << gamma[t];
	os << " ";
      }
      os << "\n\n";
    };


};

















class SymbolicSolvaccSequence: public SymbolicSequence {
public:

  int off;
  int m_format; //1: 9-line, 2: 2-line 



// This is the version taking .adataset format as input

  SymbolicSolvaccSequence(istream& is, char* the_alidir) : SymbolicSequence()
    {
      m_format = 1; //9-line 
	  strcpy(alidir, the_alidir);
      belief=0.0;
      HePl = 0;

      char c[128];

	  int t;
      is >> des;
      is >> length;

      strcpy(ali,des);

      alignments_loaded=0;


      if (length >= MAX_T-1) {
//	FAULT("Sequence too long (" << length << ") MAX_T is " << MAX_T);
      }


      u = new int[length+1];

	racc  = new Float[length+1];
	acc  = new int[length+1];
	ss  = new int[length+1];


      y = new int[length+1];
      y_pred = new int[length+1];

	  alpha = new Float[length+1];

      u[0]=y[0]=acc[0]=ss[0]=0;
	  racc[0]=0;


	  // Reads AA
      for (t=1; t<=length; t++) {
		is >> c;
		u[t]=translateU[c[0]-'A'];
      }

      // Reads SS
      for (t=1; t<=length; t++) {
		is >> c;
		if (c[0]=='-')
			ss[t]=-1;
		else if (c[0]=='H' || c[0]=='G' || c[0]=='I')
			ss[t]=0;
		else if (c[0]=='E' || c[0]=='B')
			ss[t]=1;
		else
			ss[t]=2;

      }
      // Reads (and dumps) beta partners
      for (t=1; t<=2*length; t++) {
		  is >> c;
	  }

      // Reads Relative Accessibility
      for (t=1; t<=length; t++) {
		is >> racc[t];
		racc[t] /= 100.0;
      }
      // Reads (and dumps) xyz
      for (t=1; t<=3*length; t++) {
		  is >> c;
	  }


    };

  //  SymbolicSolvaccSequence::SymbolicSolvaccSequence(istream& is, char* the_alidir, int format) : SymbolicSequence()
 SymbolicSolvaccSequence(istream& is, char* the_alidir, int format) : SymbolicSequence()
    {
      m_format = format;  
      strcpy(alidir, the_alidir);
      belief=0.0;
      HePl = 0;
      char c[128];
      int t;
      is >> des; //des: is the alignment file name
      //construct the alignment file = alignment_dir + seq_name(alignment_file_name)
      strcpy(ali,des);

      //read sequence
      char tmp[20000];
      is >> tmp; 
      length = strlen(tmp);

      alignments_loaded=0;
      //allocate memory: most of them are not necessary for format 2. 
      u = new int[length+1];
      racc  = new Float[length+1];
      acc  = new int[length+1];
      ss  = new int[length+1];
      y = new int[length+1];
      y_pred = new int[length+1];
      alpha = new Float[length+1];
      for(int i = 0; i < length + 1; i++)
      {
      	u[i]=y[i]=acc[i]=ss[i]= 0;
      	racc[i]= y_pred[i] = 0;
        alpha[i] = 0;
      }

      //translate amino acid to integer
      for (t=1; t<=length; t++) {
	u[t]=translateU[tmp[t-1]-'A'];
      }

      //Reads AA
  //    for (t=1; t<=length; t++) {
//		is >> c;
//		u[t]=translateU[c[0]-'A'];
 //     }


// Reads Relative Accessibility
//      for (t=1; t<=length; t++) {
//		is >> racc[t];
//		racc[t] /= 100.0;
//     }

    };
 /*
 //Old version taking a modified .dataset (frasc) format as input

  SymbolicSolvaccSequence(istream& is, char* the_alidir) : SymbolicSequence()
    {
	  strcpy(alidir, the_alidir);
      belief=0.0;
      HePl = 0;

      char c[4];

	  int t;
      is >> des;

      strcpy(ali,des);
	  ext_ali();
//      strncat(ali,des,4);

      alignments_loaded=0;

      is >> length;

//	  cerr << ali << " " << length << " ";

      if (length >= MAX_T-1) {
	FAULT("Sequence too long (" << length << ") MAX_T is " << MAX_T);
      }


      u = new int[length+1];

	racc  = new Float[length+1];
	acc  = new int[length+1];
	ss  = new int[length+1];


      y = new int[length+1];
      y_pred = new int[length+1];

	  alpha = new Float[length+1];

      u[0]=y[0]=acc[0]=ss[0]=0;
	  racc[0]=0;


      for (t=1; t<=length; t++) {
		char c[4];
		is >> c;
		u[t]=translateU[c[0]-'A'];
      }
      for (t=1; t<=length; t++) {
		is >> acc[t];
      }
      for (t=1; t<=length; t++) {
		is >> racc[t];
		y[t] = (int)(racc[t]*0.1);
		racc[t] /= 100.0;
      }
      for (t=1; t<=length; t++) {
		char c[4];
		is >> c;
		ss[t]=translateU[c[0]-'A'];;
      }

    };
*/

  ~SymbolicSolvaccSequence() {
      delete[] racc;
      delete[] acc;
      delete[] acc;
      delete[] u;
      delete[] y;
      delete[] y_pred;
  }





  void write(ostream& os)
    {
      int t;
      os << des << " " << length << "\n";
      for (t=1; t<=length; t++) {
		  if (u[t]>=0) {
			os << Utranslate[u[t]]<<"\t";
		  }
		  else os<<"X\t";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
		os<<acc[t]<<"\t";
		}
      os << "\n";
      for (t=1; t<=length; t++) {
		os<<racc[t]<<"\t";
	  }
      os << "\n";
      for (t=1; t<=length; t++) {
		os<<ss[t]<<"\t";
		}
      os << "\n";
      for (t=1; t<=length; t++) {
		if (y[t]==-1) os << "-\t";
		else os<<y[t]<<"\t";
		}
      os << "\n";
      for (t=1; t<=length; t++) {
		if (y_pred[t]==-1) os << "-\t";
		else os<<y_pred[t]<<"\t";
		}
      os << "\n\n";
    };


  void write_flat(ostream& os)
    {
      int t;
      os << des << " " << length << "\n";
      for (t=1; t<=length; t++) {
		  if (u[t]>=0) {
			os << Utranslate[u[t]];
		  }
		  else os<<"X";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
		if (y[t]==-1) os << "-";
		else os<<y[t];
		}
      os << "\n";
      for (t=1; t<=length; t++) {
		if (y_pred[t]==-1) os << "-";
		else os<<y_pred[t];
		}
      os << "\n\n";
    };



  void write_flat2(ostream& os)
    {
      int t;
      os << des << " " << length << "\n";
      for (t=1; t<=length; t++) {
		  if (u[t]>=0) {
			os << Utranslate[u[t]];
		  }
		  else os<<"X";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
		if (y[t]==-1) os << "-";
		else os<<y[t];
		}
      os << "\n";
      for (t=1; t<=length; t++) {
		if (y_pred[t]==-1) os << "-";
		else os<<y_pred[t];
		}
      os << "\n";
      for (t=1; t<=length; t++) {
		os<<alpha[t]<<"\t";
		}      os << "\n\n";
    };


  void write_e(ostream& os)
    {
      int t;
      os << des << " " << length << "\n";
      for (t=1; t<=length; t++) {
		  if (u[t]>=0) {
			os << Utranslate[u[t]]<<"\t";
		  }
		  else os<<"X\t";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
		os<<acc[t]<<"\t";
		}
      os << "\n";
      for (t=1; t<=length; t++) {
		os<<racc[t]<<"\t";
		}
      os << "\n";
      for (t=1; t<=length; t++) {
		if (y[t]==-1) os << "-\t";
		else os<<ss[t]<<"\t";
		}
      os << "\n\n";
    };



  void write_u(ostream& os)
    {
      int t;
//      os << des << "\n";
      for (t=1; t<=length; t++) {
	  os << Utranslate[u[t]];
      }
      os << "\n";
    };

};











class InputSequence: public Sequence {
public:
  int off;
  char* s;

//  SymbolicSequence() : Sequence() {}

  InputSequence(istream& is) : Sequence()
    {
      HePl = 0;
      char tmp[MAX_T];

      strcpy(ali,"");
      alignments_loaded=0;
	  is >> tmp;
      length = strlen(tmp);
      if (length >= MAX_T-1) {
	    FAULT("Sequence too long (" << length << ") MAX_T is " << MAX_T);
      }
//      u = new int[length+1];
//      u[0]=0;

      s=new char[length+1];
	  strcpy(s,tmp);

//      for (t=1; t<=length; t++) {
//	    u[t]=translateU[tmp[t-1]-'A'];
//      }
    };

  InputSequence(int l=MAX_T) : Sequence()
  {
      HePl = 0;
      alignments_loaded=0;

      length = l;

//      u = new int[length+1];
//	  memset(u,0,(length+1)*sizeof(int));
//      u[0]=0;

	  s = new char[length+1];
  }

  ~InputSequence() {
      delete[] s;
  }



int read_u(istream& is)
    {
      HePl = 0;
      char tmp[MAX_T];
	  strcpy(tmp,"");
      is >> tmp;
      length = strlen(tmp);
      if (length >= MAX_T-1) {
		FAULT("Sequence too long (" << length << ") MAX_T is " << MAX_T);
      }
//      u[0]=0;
//      for (t=1; t<=length; t++) {
//		u[t]=translateU[tmp[t-1]-'A'];
//      }
	strcpy(s,tmp);

	  return length;
    };


  int copy_alignment(Sequence* seq,int n, int off, int le)
  {
	return 0;
  }


  int ext_ali()
   {
	return 0;
   }



int load_alignments()
	{
	return 0;
	}

int unload_alignments()
	{
	return 0;
	}

int Aps(int al,int offs)
	{
	return 0;
	}


int best_alignments()
	{
	return 0;
	}

int best_alignments(int AL)
	{
	return 0;
	}

int generate_profile() {
	return 0;
}


void unload_profile() {

	}


  void write(ostream& os)
    {

    };

  void write_u(ostream& os)
    {
      int t;
//      os << des << "\n";
      for (t=0; t<length; t++) {
	  os << s[t];
      }
      os << "\n";
    };


  void write_predictions(ostream& os)
    {

    };

};




class DataSet {
public:
  int length;
  int NU;
  int NY;
  int totalNFrames;
  int total_alpha;
  int total_beta;
  int total_gamma;

  SymbolicSolvaccSequence** sequence;


  DataSet(int the_length) {
	  NU=20;
	  NY=3;
      totalNFrames=0;
      total_alpha=0;
      total_beta=0;
      total_gamma=0;
	  length=the_length;
      sequence = new SymbolicSolvaccSequence*[length];
  }

//format 1: title+9-line, 2: title + 2-line
  DataSet(istream& is, char* alidir, int format = 1)
    {
      is >> length;
      is >> NU;
      is >> NY;
      totalNFrames=0;
      total_alpha=0;
      total_beta=0;
      total_gamma=0;

      sequence = new SymbolicSolvaccSequence*[length];


      for (int p=0; p<length; p++) {
		 if (format == 1)
		 {
		 	 sequence[p] = new SymbolicSolvaccSequence(is, alidir);
		 }
		 else if (format == 2)
		 {
			 
		 	 sequence[p] = new SymbolicSolvaccSequence(is, alidir, format);
		 }
		 else
		 {
			cout << "data set format is not supported.\n"; 
			exit(1); 
		 }

		  for (int t=1; t<=sequence[p]->length; t++) {
			  if (sequence[p]->y[t]>=0) {
				  totalNFrames++;

				  if (sequence[p]->y[t]==0) {
					  total_alpha++;
				  }
				  if (sequence[p]->y[t]==1) {
					  total_beta++;
				  }
			  }
		  }
	  }
  };



  void write(ostream& os)
    {
      os << length << " " << NU << " " << NY << "\n";
      for (int p=0; p<length; p++) {
	sequence[p]->write(os);
      }
    };



  void write(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
	ostream os(&outbuf);
	this->write(os);
      } else {
	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };


  void write_flat(ostream& os)
    {
      os << length << " " << NU << " " << NY << "\n";
      for (int p=0; p<length; p++) {
	sequence[p]->write_flat(os);
      }
    };

  void write_flat(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
	ostream os(&outbuf);
	this->write_flat(os);
      } else {
	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };

  void write_flat2(ostream& os)
    {
      os << length << " " << NU << " " << NY << "\n";
      for (int p=0; p<length; p++) {
	sequence[p]->write_flat2(os);
      }
    };

  void write_flat2(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
	ostream os(&outbuf);
	this->write_flat2(os);
      } else {
	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };




  void write_e(ostream& os)
    {
      os << length << " " << NU << " " << NY << "\n";
      for (int p=0; p<length; p++) {
	sequence[p]->write_e(os);
      }
    };

  void write_e(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
	ostream os(&outbuf);
	this->write_e(os);
      } else {
	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };



  void write_predictions(ostream& os)
    {
      os << length << " " << NU << " " << NY << "\n";
      for (int p=0; p<length; p++) {
	sequence[p]->write_predictions(os);
      }
    };



  void write_predictions(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
	ostream os(&outbuf);
	this->write_predictions(os);
      } else {
	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };



  void write3(ostream& os)
    {
      os << length << " " << NU << " " << NY << "\n";
      for (int p=0; p<length; p++) {
	sequence[p]->write3(os);
      }
    };



  void write3(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
	ostream os(&outbuf);
	this->write3(os);
      } else {
	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };


  void write3_predictions(ostream& os)
    {
      os << length << " " << NU << " " << NY << "\n";
      for (int p=0; p<length; p++) {
	sequence[p]->write3_predictions(os);
      }
    };



  void write3_predictions(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
	ostream os(&outbuf);
	this->write3_predictions(os);
      } else {
	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };













  void set_belief(Float be) {
	  for (int p=0;p<length;p++) {
		sequence[p]->set_belief(be);
	  }
  };
};





class AlignmentsSet {
public:
  int length;
  int totalNFrames;
  InputSequence** sequence;


  AlignmentsSet(istream& is)
    {
      is >> length;
      totalNFrames=0;

      sequence = new InputSequence*[length];

      for (int p=0; p<length; p++) {
		sequence[p] = new InputSequence(is);
	    totalNFrames+=sequence[p]->length;
		}
    };



  void write(ostream& os)
    {
      os << length << "\n";
      for (int p=0; p<length; p++) {
	    sequence[p]->write(os);
      }
    };



  void write(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
	ostream os(&outbuf);
	this->write(os);
      } else {
	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };

/*
  void write_predictions(ostream& os)
    {
      os << length << " " << NU << " " << NY << "\n";
      for (int p=0; p<length; p++) {
	sequence[p]->write_predictions(os);
      }
    };



  void write_predictions(char* fname)
    {
      filebuf outbuf;
      if (outbuf.open(fname, ios::out) != 0) {
	ostream os(&outbuf);
	this->write_predictions(os);
      } else {
	FAULT("Failed to write to file " << fname);
      }
      outbuf.close();
    };
*/
};




#endif // Sequence_h
