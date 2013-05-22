#ifndef Sequence_h
#define Sequence_h 1

#include <iostream.h>
#include "General.h"
#include <vector>
#include <math.h>

#define MAX_T 8192

class Sequence {
public:
  char ali[512];
  char alidir[512];
  int length;
  int* u; //sequence
  int* y; //secondary structure?
  int* y_pred;

//add other information for sspro
	int *b1; //binding site 1
	int *b2; //binding site 2
	int *accessibility; //solvent accessibility
	Float * cx; //x coordinate
	Float * cy; //y coordinate
	Float * cz; //z coordinate
//End of addition, Nov 20, 2003, Jianlin Cheng

  Float* alpha;
  Float* beta;
  Float* gamma;

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

  Sequence(istream& is, char* the_alidir);

  ~Sequence(); 

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
  virtual void write(ostream& os);
  virtual void write_u(ostream& os)=0;

//get secondary structure of amino acid(for format 0, 1, 2: return true secondary structure, format 3: return predicted ss.
//index: the index of amino acid, start from 1; 
  virtual int 	getSS(int index){return 0;};
//index: the index of amino acid, start from 1; 
  virtual Float getX(int index){return 0;};
//index: the index of amino acid, start from 1; 
  virtual Float getY(int index){return 0;};
//index: the index of amino acid, start from 1; 
  virtual Float getZ(int index){return 0;}; 

//check if two amino acids are contacted or not 
//parameters: index of aa1, index of aa2, threshold. indext start from 1
  virtual bool isContacted(int index1, int index2, int threshold){return false;}; 


};

class SymbolicSequence: public Sequence {

private:
  	int* y_sspro; 
	//for 3d information
	int m_format; 
	//store the neighbors for each amino acid. 
	std::vector< std::vector<int> > m_8aNeighbors; 
	std::vector< std::vector<int> > m_12aNeighbors; 
	std::vector< std::vector<int> > m_16aNeighbors; 
private:
	//this function read foramt 4, only be used by prediction program. 
	//format: seq_name, seq_length, seq, predicted_contacts(2 lines), predicted_secondary_structure
	void readSpecialFormat(istream& is); 
	void readFormat5(istream& is); 
public:
  int off;
  char des[512];
  char des2[512];
  char des3[512];

// this constructor is used for the old format data set. 
  SymbolicSequence(istream& is, char* the_alidir); 

//this constructor is used for the 9-line new format
  SymbolicSequence(istream& is, char* the_alidir, int format); 

  SymbolicSequence(int l); 

  ~SymbolicSequence(); 

void set_belief(Float be); 

int read_u(istream& is);

  int copy_alignment(Sequence* seq,int n, int off, int le);

  int ext_ali();

int load_alignments();

int unload_alignments();

int Aps(int al,int offs);

int best_alignments();

int best_alignments(int AL);

#ifdef NEVER
int generate_profile(); 
#endif

Float profile_entropy(); 

int generate_profile(); 

void unload_profile(); 

void write(ostream& os);

void write_u(ostream& os);

void write_predictions(ostream& os);

//add for 3d information
  virtual int 	getSS(int index);
  virtual Float getX(int index);
  virtual Float getY(int index);
  virtual Float getZ(int index); 
  virtual bool isContacted(int index1, int index2, int threshold); 

private:
   Float calcDistance(Float x1, Float y1, Float z1, Float x2, Float y2, Float z2);
};

class InputSequence: public Sequence {
public:
  int off;
  char* s;

  InputSequence(istream& is); 

  InputSequence(int l); 

  ~InputSequence(); 

int read_u(istream& is);
  int copy_alignment(Sequence* seq,int n, int off, int le);

  int ext_ali();

int load_alignments();

int unload_alignments();

int Aps(int al,int offs);

int best_alignments();

int best_alignments(int AL);

int generate_profile(); 

void unload_profile(); 

  void write(ostream& os);

  void write_u(ostream& os);

  void write_predictions(ostream& os);

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
  SymbolicSequence** sequence;

  DataSet(istream& is, char* alidir, int format);

  void write(ostream& os);

  void write(char* fname);

  void write_predictions(ostream& os);

  void write_predictions(char* fname);

  void set_belief(Float be); 
};


class AlignmentsSet {
public:
  int length;
  int totalNFrames;
  InputSequence** sequence;
  AlignmentsSet(istream& is);
  void write(ostream& os);
  void write(char* fname);
};


#endif // Sequence_h
