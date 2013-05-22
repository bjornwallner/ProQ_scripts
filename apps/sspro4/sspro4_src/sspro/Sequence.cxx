#include "Sequence.h"
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

Sequence::Sequence(istream& is, char* the_alidir)
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


Sequence::~Sequence() {
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


 void Sequence::write(ostream& os)
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

SymbolicSequence::SymbolicSequence(istream& is, char* the_alidir) : Sequence()
    {
	m_format = 0; 
	  strcpy(alidir,the_alidir);
	  belief=0.0;
      HePl = 0;
      char tmp[MAX_T];
      int t;
      is >> des;
      is >> tmp;

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
	else if (tmp[t-1]=='E')
	  y[t]=1;
	else
	  y[t]=2;
      }
    };


  SymbolicSequence::SymbolicSequence(istream& is, char* the_alidir, int format)
    {
      strcpy(alidir,the_alidir);
      m_format = format; 
      if (m_format == 4)
      {
	readSpecialFormat(is);
	return; 
      }
      belief = 0.0; 
      HePl=0;

      int t;
      is >> des;

      strcpy(ali, "");
      ext_ali(); 

      alignments_loaded=0;
      is >> length;

      if (length >= MAX_T-1) {
		FAULT("Sequence too long (" << length << ") MAX_T is " << MAX_T);
      }

      u = new int[length+1];
      y = new int[length+1];
      b1 = new int[length+1];
      b2 = new int[length+1];
      accessibility = new int[length+1];
      cx = new Float[length+1];
      cy = new Float[length+1];
      cz = new Float[length+1];
      y_pred = new int[length+1];


      alpha = new Float[length+1];
      beta = new Float[length+1];
      gamma = new Float[length+1];


      entropy = new Float[length+1];
      present = new Float[length+1];
      settato = new int[length+1];

      u[0]=y[0]=b1[0]=b2[0]=0;
      cx[0]=cy[0]=cz[0]=0.0;
	  
      char ch;
      int b;

      for (t=1; t<=length; t++) {
		is >> ch;
		u[t]=translateU[ch-'A'];
      }

      for (t=1; t<=length; t++) {
		settato[t] = 0; 
		is >> ch;
		if (ch=='-')
			y[t]=-1;
		else if (ch=='H' || ch=='G' || ch=='I')
			y[t]=0;
		else if (ch=='E' || ch=='B')
			y[t]=1;
		else
			y[t]=2;
      }

      for (t=1; t<=length; t++) {
		is >> b;
		if (b<0) b=0;
		b1[t]=b;
      }

      for (t=1; t<=length; t++) {
		is >> b;
		if (b<0) b=0;
		b2[t]=b;
      }

      for (t=1; t<=length; t++) {
		is >> accessibility[t];
      }

      for (t=1; t<=length; t++) {
		is >> cx[t] >> cy[t] >> cz[t];
      }

	if (format == 1)
	{
		return; 
	}
	if (format == 5)
	{
		readFormat5(is); 
		return; 
	}
	

	//read neighbor information from fine contact map prediction
	//read 8 angstrom
      for (t=1; t<=length; t++) 
      {
	//read index of aa (start from 1)
	int index;
	is >> index;
	if (index < 1 || index > length)
	{
		cout << "index of 3d information is out of boundary.\n";
		exit(1);
	} 
	//read aa
	char aa;
	is >> aa;
	int taa = translateU[aa - 'A'];
	if (u[t] != taa)
	{
		cout << "index " << t << " u[t]=" << u[t] << " aa=" << taa << " " << aa << endl;  
		cout << "amino acid in 3d doesn't match with the aa in sequence.\n";
		exit(1); 
	}
	//read number of neighbors
	int num;
	is >> num;
	std::vector<int> neighbors;
	if (num > 0)
	{
		//read neighbors
		for (int j = 0; j < num; ++j)
		{
			int nei;
			is >> nei;
			neighbors.push_back(nei); 
		}
		
	}
	//store the neighbors
	m_8aNeighbors.push_back(neighbors); 
      } //end of read 3d info of 8A
	
	//read 12 angstrom
      for (t=1; t<=length; t++) 
      {
	//read index of aa (start from 1)
	int index;
	is >> index;
	if (index < 1 || index > length)
	{
		cout << "index of 3d information is out of boundary.\n";
		exit(1);
	} 
	//read aa
	char aa;
	is >> aa;
	int taa = translateU[aa - 'A'];
	if (u[t] != taa)
	{
		cout << "index " << t << " u[t]=" << u[t] << " aa=" << taa << " " << aa << endl;  
		cout << "amino acid in 3d doesn't match with the aa in sequence.\n";
		exit(1); 
	}
	//read number of neighbors
	int num;
	is >> num;
	std::vector<int> neighbors;
	if (num > 0)
	{
		//read neighbors
		for (int j = 0; j < num; ++j)
		{
			int nei;
			is >> nei;
			neighbors.push_back(nei); 
		}
		
	}
	//store the neighbors
	m_12aNeighbors.push_back(neighbors); 
      } //end of read 3d info of 12A

	//read 16 angstrom
      for (t=1; t<=length; t++) 
      {
	//read index of aa (start from 1)
	int index;
	is >> index;
	if (index < 1 || index > length)
	{
		cout << "index of 3d information is out of boundary.\n";
		exit(1);
	} 
	//read aa
	char aa;
	is >> aa;
	int taa = translateU[aa - 'A'];
	if (u[t] != taa)
	{
		cout << "index " << t << " u[t]=" << u[t] << " aa=" << taa << " " << aa << endl;  
		cout << "amino acid in 3d doesn't match with the aa in sequence.\n";
		exit(1); 
	}
	//read number of neighbors
	int num;
	is >> num;
	std::vector<int> neighbors;
	if (num > 0)
	{
		//read neighbors
		for (int j = 0; j < num; ++j)
		{
			int nei;
			is >> nei;
			neighbors.push_back(nei); 
		}
		
	}
	//store the neighbors
	m_16aNeighbors.push_back(neighbors); 
      } //end of read 3d info of 16A
	if (m_format == 2)
	{
		return; 
	}
	if (m_format == 3)
	{
		//read in the predicted secondary structure from old sspro
		y_sspro = new int[length + 1]; 
      		for (t=1; t<=length; t++) {
			is >> ch;
			if (ch=='-')
				y_sspro[t]=-1;
			else if (ch=='H' || ch=='G' || ch=='I')
				y_sspro[t]=0;
			else if (ch=='E' || ch=='B')
				y_sspro[t]=1;
			else
				y_sspro[t]=2;
      		}
	}

  };


void SymbolicSequence::readFormat5(istream& is)
{
	
	if (m_format != 5)
	{
		cout << "format 5 is wrong in readFormat5()\n";
		exit(1); 
	}
	int t; 
	//read 12 angstrom
      for (t=1; t<=length; t++) 
      {
	//read index of aa (start from 1)
	int index;
	is >> index;
	if (index < 1 || index > length)
	{
		cout << "index of 3d information is out of boundary.\n";
		exit(1);
	} 
	//read aa
	char aa;
	is >> aa;
	int taa = translateU[aa - 'A'];
	if (u[t] != taa)
	{
		cout << "index " << t << " u[t]=" << u[t] << " aa=" << taa << " " << aa << endl;  
		cout << "amino acid in 3d doesn't match with the aa in sequence.\n";
		exit(1); 
	}
	//read number of neighbors
	int num;
	is >> num;
	std::vector<int> neighbors;
	if (num > 0)
	{
		//read neighbors
		for (int j = 0; j < num; ++j)
		{
			int nei;
			is >> nei;
			neighbors.push_back(nei); 
		}
		
	}
	//store the neighbors
	m_12aNeighbors.push_back(neighbors); 
      } //end of read 3d info of 12A

	char ch; 
		//read in the predicted secondary structure from old sspro
		y_sspro = new int[length + 1]; 
      		for (t=1; t<=length; t++) {
			is >> ch;
			if (ch=='-')
				y_sspro[t]=-1;
			else if (ch=='H' || ch=='G' || ch=='I')
				y_sspro[t]=0;
			else if (ch=='E' || ch=='B')
				y_sspro[t]=1;
			else
				y_sspro[t]=2;
      		}
}
 
  void SymbolicSequence::readSpecialFormat(istream& is)
    {
	if (m_format != 4)
	{
		cout << "the sequence is not in format 4.\n"; 
		exit(1); 
	}
      belief = 0.0; 
      HePl=0;
      int t;
      //read sequence name
      is >> des; //des name is sequence name, also the alignment file name. 
      strcpy(ali, "");
      ext_ali(); 
      alignments_loaded=0;

      //read sequence length
      is >> length;
      if (length >= MAX_T-1) {
		FAULT("Sequence too long (" << length << ") MAX_T is " << MAX_T);
      }

      //initialize variables, many of them are not used acutally. 
      u = new int[length+1];
      y = new int[length+1]; //not used
      b1 = new int[length+1]; //not used
      b2 = new int[length+1]; //not used
      accessibility = new int[length+1]; //not used
      cx = new Float[length+1]; //not used
      cy = new Float[length+1]; //not used
      cz = new Float[length+1]; //not used
      y_pred = new int[length+1];
      alpha = new Float[length+1];//not used
      beta = new Float[length+1]; //not used
      gamma = new Float[length+1]; //not used
      entropy = new Float[length+1]; //not used
      present = new Float[length+1]; //not used
      settato = new int[length+1]; //not used
      for (t = 0; t <= length; t++)
      {
	u[t]=y[t]=b1[t]=b2[t]=accessibility[t]=0;
        cx[t]=cy[t]=cz[t]=0;
        y_pred[t]=0;
        alpha[t]=beta[t]=gamma[t]= 0; 
	entropy[t]=present[t]=0;
        settato[t]=0; 
      }

      u[0]=y[0]=b1[0]=b2[0]=0;
      cx[0]=cy[0]=cz[0]=0.0;
	  
      char ch;
      int b;
	//read sequence
      for (t=1; t<=length; t++) {
		is >> ch;
		u[t]=translateU[ch-'A'];
      }

	//read neighbor information from fine contact map prediction
	//read 12 angstrom
      for (t=1; t<=length; t++) 
      {
	//read index of aa (start from 1)
	int index;
	//read index of amino acid
	is >> index;
	if (index < 1 || index > length)
	{
		cout << "index of 3d information is out of boundary.\n";
		exit(1);
	} 
	char aa;
	//read aa
	is >> aa;
	int taa = translateU[aa - 'A'];
	if (u[t] != taa)
	{
		cout << "index " << t << " u[t]=" << u[t] << " aa=" << taa << " " << aa << endl;  
		cout << "amino acid in 3d doesn't match with the aa in sequence.\n";
		exit(1); 
	}
	int num;
	//read number of neighbors
	is >> num;
	std::vector<int> neighbors;
	if (num > 0)
	{
		//read neighbors
		for (int j = 0; j < num; ++j)
		{
			int nei;
			is >> nei;
			neighbors.push_back(nei); 
		}
		
	}
	//store the neighbors
	m_12aNeighbors.push_back(neighbors); 
      } //end of read 3d info of 12A

	//read 16 angstrom
      for (t=1; t<=length; t++) 
      {
	//read index of aa (start from 1)
	int index;
	is >> index;
	if (index < 1 || index > length)
	{
		cout << "index of 3d information is out of boundary.\n";
		exit(1);
	} 
	//read aa
	char aa;
	is >> aa;
	int taa = translateU[aa - 'A'];
	if (u[t] != taa)
	{
		cout << "index " << t << " u[t]=" << u[t] << " aa=" << taa << " " << aa << endl;  
		cout << "amino acid in 3d doesn't match with the aa in sequence.\n";
		exit(1); 
	}
	//read number of neighbors
	int num;
	is >> num;
	std::vector<int> neighbors;
	if (num > 0)
	{
		//read neighbors
		for (int j = 0; j < num; ++j)
		{
			int nei;
			is >> nei;
			neighbors.push_back(nei); 
		}
		
	}
	//store the neighbors
	m_16aNeighbors.push_back(neighbors); 
      } //end of read 3d info of 16A

	//read in the predicted secondary structure from old sspro
	y_sspro = new int[length + 1]; 
      	for (t=1; t<=length; t++) {
		is >> ch;
		if (ch=='-')
			y_sspro[t]=-1;
		else if (ch=='H' || ch=='G' || ch=='I')
			y_sspro[t]=0;
		else if (ch=='E' || ch=='B')
			y_sspro[t]=1;
		else
			y_sspro[t]=2;
      	}

  }; 
int 	SymbolicSequence::getSS(int index)
  {

     if (index >= 1 && index <= length)
     {
	if (m_format != 3 && m_format != 4 && m_format != 5)
		//return true secondary structure
          	return y[index]; 
	else
		//return predicted secondary structure if format is 13 -line or format 4
		return y_sspro[index]; 
     }
     else
     {
	cout << "index error in SymbolicSequence::getSS\n";
	exit(1); 
     }
  }

  Float SymbolicSequence::getX(int index)
  {
     if (index >= 1 && index <= length)
     {
          return cx[index]; 
     }
     else
     {
	cout << "index error in SymbolicSequence::getX\n";
	exit(1); 
     }
  }

  Float  SymbolicSequence::getY(int index)
  {
     if (index >= 1 && index <= length)
     {
          return cy[index]; 
     }
     else
     {
	cout << "index error in SymbolicSequence::getY\n";
	exit(1); 
     }
  }

  Float SymbolicSequence::getZ(int index)
  {
     if (index >= 1 && index <= length)
     {
          return cz[index]; 
     }
     else
     {
	cout << "index error in SymbolicSequence::getZ\n";
	exit(1); 
     }
  } 

  SymbolicSequence::SymbolicSequence(int l=MAX_T) : Sequence()
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

 SymbolicSequence:: ~SymbolicSequence() {
      delete[] u;
      delete[] y;

      delete [] b1;
      delete [] b2;
      delete [] accessibility;
      delete [] cx;
      delete [] cy;
      delete [] cz;

      delete[] y_pred;
      delete[] alpha;
      delete[] beta;
      delete[] gamma;

      delete[] entropy;
      delete[] present;
      delete[] settato;
	if (m_format == 3 || m_format == 4 || m_format == 5) //for 13-line format data set or format 4
	{
		delete [] y_sspro; 
	}
  }


void SymbolicSequence::set_belief(Float be) {
	belief=be;
	};

int SymbolicSequence::read_u(istream& is)
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



  int SymbolicSequence::copy_alignment(Sequence* seq,int n, int off, int le)
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

  int SymbolicSequence::ext_ali()
   {
   char* pi=des;
   int cont=0;
   while (pi[0]!='.' && pi[0]!=' ' && pi[0]!= '\0'  && cont<512)
	{
	cont++;
	strncat(ali,pi,1);
	pi++;
	}
   return 1;
   }

int SymbolicSequence::load_alignments()
	{
  	filebuf inbuf;
	char temp[4096];
  	char fname[1024];

	strcpy(fname, alidir);
  	strcat(fname, ali);

	if (inbuf.open(fname, ios::in) != 0) {
    	  istream is(&inbuf);
	  is >> alignments_loaded;
//	cout << alignments_loaded << "\n" << flush;
	  ALIGNMENTS=new int*[alignments_loaded];
	  SCORES=new int[alignments_loaded];
	  INTERS=new int[alignments_loaded];
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
	 cout << "alignment file: " << fname << "does not exist." << endl; 
	#endif
    	 alignments_loaded=0;
  	}
  	inbuf.close();

	return alignments_loaded;
	}

int SymbolicSequence::unload_alignments()
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



int SymbolicSequence::Aps(int al,int offs)
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



int SymbolicSequence::best_alignments()
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

int SymbolicSequence::best_alignments(int AL)
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

int SymbolicSequence::generate_profile() {

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

Float SymbolicSequence::profile_entropy() {
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

int SymbolicSequence::generate_profile() {

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

	flags=new int[AL];
	offsets=new int[AL];
	minls=new int[AL];
	memset(flags,0,AL*sizeof(int));
	memset(offsets,0,AL*sizeof(int));
	memset(minls,0,AL*sizeof(int));

	for (a=0;a<AL;a++) {
	   offset = best_alignments(a);
	   int sum=Aps(a,offset);
	   Float qual = (Float)sum/(Float)INTERS[a];

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

Float* weights=new Float[AL];

//cout << profile_entropy() << "\n" << flush;

for (int iter=0;iter<1;iter++) {

// We can now think to weigh each sequence based on its information content

memset(weights, 0, AL*sizeof(Float));

for (a=0;a<AL;a++) {
	if (!flags[a])
		{continue;}

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
delete[] flags;
delete[] offsets;
delete[] minls;
delete[] weights;

return 0;

}



void SymbolicSequence::unload_profile() {
	HePl = 0;
	delete[] HeP;
	}



  void SymbolicSequence::write(ostream& os)
    {
      int t;
      os << des << "\n";
      for (t=1; t<=length; t++) {
	os << Utranslate[u[t]];
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y[t]==-1) os << "-";
	else if (y[t]==0) os << "H";
	else if (y[t]==1) os << "E";
	else if (y[t]==2) os << "C";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y_pred[t]==-1) os << "-";
	else if (y_pred[t]==0) os << "H";
	else if (y_pred[t]==1) os << "E";
	else if (y_pred[t]==2) os << "C";
      }
      os << "\n\n";
    };

  void SymbolicSequence::write_u(ostream& os)
    {
      int t;
//      os << des << "\n";
      for (t=1; t<=length; t++) {
	  os << Utranslate[u[t]];
      }
      os << "\n";
    };


  void SymbolicSequence::write_predictions(ostream& os)
    {
      int t;
      os << des << "\n";
      for (t=1; t<=length; t++) {
	os << Utranslate[u[t]];
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y[t]==-1) os << "-";
	else if (y[t]==0) os << "H";
	else if (y[t]==1) os << "E";
	else if (y[t]==2) os << "C";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y_pred[t]==-1) os << "-";
	else if (y_pred[t]==0) os << "H";
	else if (y_pred[t]==1) os << "E";
	else if (y_pred[t]==2) os << "C";
      }
      os << "\n";
      for (t=1; t<=length; t++) {
	if (y_pred[t]==-1) os << "0 0 0";
	else os << alpha[t] << " " << beta[t] << " " << gamma[t];
	os << " ";
      }
      os << "\n\n";
    };



InputSequence::InputSequence(istream& is) : Sequence()
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

      s=new char[length+1];
	  strcpy(s,tmp);

    };



  InputSequence::InputSequence(int l=MAX_T) : Sequence()
  {
      HePl = 0;
      alignments_loaded=0;

      length = l;

	  s = new char[length+1];
  }

  InputSequence::~InputSequence() {
      delete[] s;
  }



int InputSequence::read_u(istream& is)
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


  int InputSequence::copy_alignment(Sequence* seq,int n, int off, int le)
  {
	return 0;
  }


  int InputSequence::ext_ali()
   {
	return 0;
   }



int InputSequence::load_alignments()
	{
	return 0;
	}

int InputSequence::unload_alignments()
	{
	return 0;
	}

int InputSequence::Aps(int al,int offs)
	{
	return 0;
	}


int InputSequence::best_alignments()
	{
	return 0;
	}

int InputSequence::best_alignments(int AL)
	{
	return 0;
	}

int InputSequence::generate_profile() {
	return 0;
}


void InputSequence::unload_profile() {

	}


  void InputSequence::write(ostream& os)
    {

    };

  void InputSequence::write_u(ostream& os)
    {
      int t;
      for (t=0; t<length; t++) {
	  os << s[t];
      }
      os << "\n";
    };


  void InputSequence::write_predictions(ostream& os)
    {

    };


  DataSet::DataSet(istream& is, char* alidir, int format)
    {
      is >> length;
      is >> NU;
      is >> NY;
      totalNFrames=0;
      total_alpha=0;
      total_beta=0;
      total_gamma=0;

      sequence = new SymbolicSequence*[length];

	//cout << "length: " << length << " format: " << format << endl; 


      for (int p=0; p<length; p++) {
		if (format == 0)
		{
			//3-line old format dataset
			sequence[p] = new SymbolicSequence(is, alidir);
		}
		else  //9-line or 12 line format, 11 line or more
		{
			//9-line new format dataset
			sequence[p] = new SymbolicSequence(is, alidir, format); 
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

	    if (sequence[p]->y[t]==2) {
		total_gamma++;
		}

	    }
	}
      }
	


    };



  void DataSet::write(ostream& os)
    {
      os << length << " " << NU << " " << NY << "\n";
      for (int p=0; p<length; p++) {
	sequence[p]->write(os);
      }
    };



  void DataSet::write(char* fname)
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


  void DataSet::write_predictions(ostream& os)
    {
      os << length << " " << NU << " " << NY << "\n";
      for (int p=0; p<length; p++) {
	sequence[p]->write_predictions(os);
      }
    };



  void DataSet::write_predictions(char* fname)
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


  void DataSet::set_belief(Float be) {

	  for (int p=0;p<length;p++) {
		sequence[p]->set_belief(be);
	  }
  };




AlignmentsSet::AlignmentsSet(istream& is)
{
      is >> length;
      totalNFrames=0;

      sequence = new InputSequence*[length];

      for (int p=0; p<length; p++) {
		sequence[p] = new InputSequence(is);
	    totalNFrames+=sequence[p]->length;
		}
    };

  void AlignmentsSet::write(ostream& os)
    {
      os << length << "\n";
      for (int p=0; p<length; p++) {
	    sequence[p]->write(os);
      }
    };

  void AlignmentsSet::write(char* fname)
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

Float SymbolicSequence::calcDistance(Float x1, Float y1, Float z1, Float x2, Float y2, Float z2)
{
	double total = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2); 
	return sqrt(total); 
}

bool SymbolicSequence::isContacted(int index1, int index2, int threshold)
{
	if (index1 < 1 || index2 < 1 || index1 > length || index2 > length)
	{
		return false; 
	}
	if (m_format == 0)
	{
		cout << "format 0 sequence doesn't have contact informaton.\n";
		exit(1); 
	}

	if (m_format == 1) //true 3d
	{	

		float distance = calcDistance(getX(index1), getY(index1), getZ(index1), 
						getX(index2), getY(index2),getZ(index2) ); 
		if (distance < threshold)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else if (m_format == 2 || m_format == 3 || m_format == 4 || m_format == 5) //using predicted 3d information
	{
	//	cout << "check contact  " << m_8aNeighbors.size() << " " << m_12aNeighbors.size()
	//	<< " " << m_16aNeighbors.size() <<   endl; 
		if (threshold != 8 && threshold != 12 && threshold != 16)
		{
			return false; 
		}
		if (threshold == 8 && m_format == 4)
		{
			cout << "format 4 sequence doesn't have 8-A neighbors. exit.\n";
			exit(1); 
		}
		if (threshold != 12 && m_format == 5)
		{
			cout << "format 5 sequence doesn't have neighbors for this threshold. exit.\n";
			exit(1); 
		}
		std::vector<int>  neighbors;
		switch (threshold)
		{
			case 8: neighbors = m_8aNeighbors[index1 - 1]; break;
			case 12: neighbors = m_12aNeighbors[index1 - 1]; break;
			case 16: neighbors = m_16aNeighbors[index1 - 1]; break;
			default: cout << "threshold wrong in Sequence::isContacted\n"; exit(1); 	
		}
		//check if index2 is the neighbor of index 1; 
		bool bFound = false; 
		for (int i = 0; i < neighbors.size(); ++i)
		{
			if (index2 == neighbors[i])
			{
				bFound = true;
				break; 
				//cout << "size " <<  neighbors.size() << endl; 
				//cout << "neighbor..." << endl;
			}
		}
		return bFound; 
	}
	else
	{
		return false; 
	}
}
