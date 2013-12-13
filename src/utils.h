/*
    Copyright (C) 2009 Stefan Enroth

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/



#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
// needed to add these to compile on Uppmax ( linux, 64 bit)
#include <string.h>
#include <limits.h>
#include <iomanip>
#include <vector>
#include <iterator>
#include <map>
// itoa & random stuff
#include <ctime>
#include <sys/time.h>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cfloat>
#include <math.h>



/* this defines how the data is stored, easy to change in the Makefile if 
 * larger numbers are needed. 
 *  
 * 
 *
 */


#ifdef INT_REP // unsigned interger

typedef unsigned int storageType; 
const int MAXNR = UINT_MAX; 
const int CHUNK_SIZE = 32*USHRT_MAX;

#elif defined (DBL_REP) //doubles

typedef double storageType; 
const int MAXNR = DBL_MAX;
const int CHUNK_SIZE = 32*USHRT_MAX;

#elif defined (SHR_REP) // unsigned short

typedef unsigned short storageType; 
const int MAXNR = USHRT_MAX;
const int CHUNK_SIZE = 32*USHRT_MAX;

#else   // default, unsigned short

typedef unsigned short storageType; 
const int MAXNR = USHRT_MAX;
const int CHUNK_SIZE = 32*USHRT_MAX;

#endif

using namespace std;

/*
 * spits the current string on delimiters
 *
 */ 
void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters = " ");

/*
 * as Tokenize but returns the 'next' position to start from.
 *
 */ 
string::size_type readAndTokenize(const string& str,
				  vector<string>& tokens,
				  const string& delimiters = " ",
				  string::size_type firstPos = 0);

/***********************************************************
 *
 * Representation of the sequences and their properties. 
 *
 *
 **********************************************************/

class fsHolder
{
 public:
  int minc;
  string nameR;
  string nameF;
  string nameC;
  fstream * fsR;
  fstream * fsF;
  fstream * fsC;
  bool readOnly; 
  
  
  // constructor
  fsHolder(){
    //cerr<<"fsH constr."<<endl;
    minc = 0;
    nameR.assign("R.bin");
    nameF.assign("F.bin");
    nameC.assign("C.bin");
    readOnly = true;
  }
  // param. constr.
  fsHolder(string r,string f,string c,int low){
    //cerr<<"fsH param. constr."<<endl;
    minc = low;
    nameR.assign(r);
    nameF.assign(f);
    nameC.assign(c);
  }
  fsHolder(string r,string f,string c,int low,bool ro){
    //cerr<<"fsH param. constr."<<endl;
    minc = low;
    nameR.assign(r);
    nameF.assign(f);
    nameC.assign(c);
    readOnly = ro;
  }

  //copy constr.
  fsHolder(const fsHolder & fsH){
    //cerr<<"fsH copy. constr."<<endl;
    minc  = fsH.minc;
    nameR = fsH.nameR;
    nameF = fsH.nameF;
    nameC = fsH.nameC;
    readOnly = fsH.readOnly;
    if(fsH.readOnly)
      fsR = new fstream(nameR.c_str(),ios::in | ios::binary);
    else
      fsR = new fstream(nameR.c_str(),ios::out | ios::in | ios::trunc | ios::binary);
    if(fsH.readOnly)
      fsF = new fstream(nameF.c_str(),ios::in | ios::binary);
    else
      fsF = new fstream(nameF.c_str(),ios::out | ios::in | ios::trunc | ios::binary);
    if(strcmp(nameC.c_str(),"") != 0)
      if(fsH.readOnly)
	fsC = new fstream(nameC.c_str(),ios::in | ios::binary);
      else
	fsC = new fstream(nameC.c_str(),ios::out | ios::in | ios::trunc | ios::binary);
  }
  // destruct.
  ~fsHolder()
    {
      //cerr<<"fsH destr."<<endl;
      fsR->close();
      delete fsR;
      fsF->close();
      delete fsF;
      if(strcmp(nameC.c_str(),"") != 0)
	{
	  fsC->close();
	  delete fsC;
	}
    }
};


/*
 * used to keep track of some statistics of each sequence.
 *  
 */
struct seqStats{
  int minPos;
  int maxPos;
  int truncR;
  int truncF;
  int truncC;
  int countF;
  int countR;
};

/*
 *
 * used to store data from a single parsed input-line. 
 *
 */
struct inputLine{
  int strand;      // 1/-1
  int pos;         // 1-based leftmost position.
  int len;         // length of mapped fragment, determined by seq. 
  string seq;      // name of reference sequence.  
  string name;     // name of the BED line.
  int mapped;      // was this read mapped? (1/0)
  int header;      // is this line a header line (1/0), obsoletes the above. 
  storageType val; // if input specified values are used. 
};


/************************************************************
 *
 * Binary file accessors. 
 * 
 ************************************************************/


/************************************************************
 *
 * Calculate min & max position in a given file (or stream).
 *
 * writes parameters to the given int-pointers.
 * returns 0, -1 on error. 
 * 
 ************************************************************/

int getMinMaxC(const char* fname, int* min, int* max);
int getMinMaxF(ifstream * inf, int* min, int* max);

/************************************************************
 *
 * Access the given file and retrieve the signal in the given
 * coordinates. 
 *
 * Fails (returns -1) if the coordinates are illegal.   
 *
 *
 ************************************************************/

int getSignalC(const char* fname, int from, int to, storageType *signal);
int getSignalF(ifstream * inf, int from, int to, storageType * signal);


/************************************************************
 *
 * Access a given binary file, retrives all signals given a 
 * text-file with two columns start & stop. 
 * writes all failed access attempts to a second file. 
 *
 * puts the successful reads in another text/binary file, 
 * assumes all are of equal length (hard to sort out that file 
 * otherwise.)
 * 
 * returns the number of successful reads.
 *
 ************************************************************/

int getSignalsBin(const char * bin_file,const char* q_file, const char* out_file,const char* err_file);
int getSignalsTxt(const char* bin_file,const char* q_file, const char* out_file,const char* err_file,int avg,const char* err_symb);

/************************************************************
 *
 * Access a given binary file, retrives all signals given a 
 * text-file with two columns start & strand. 
 * writes all failed access attempts to a second file. 
 *
 * puts the avearage of successful reads split on strand in 
 * a text file. signals are read from the given postions +/- offset. 
 * the move parameter shifts the read window up or down.  
 * 
 * returns the total number of successful reads.
 *
 ************************************************************/

int getSignalsFP(const char* bin_file,const char* fc_bin_file,const char* q_file, const char* out_file,
		 const char* err_file,int offset, int progress, int co, int move);

int parseQueryLine(string str,inputLine*);

int queryControl(ifstream * inf,map<string,seqStats> * seqmap,bool verbose);

int queryControlBED(ifstream * inf,map<string,seqStats> * seqmap,bool verbose, string type);

int getMultipleSignalsFP(ifstream * infq,map<string,fsHolder>,string ofStub,
			 string eFile,int offset,int progress,int co,int move,
			 int both,int nq, int all);

int getMultipleSignalsFPfromBED(ifstream * infq,map<string,fsHolder>,string ofStub,
				string eFile,int offset,int progress,int co,int move,
				int both,int nq, int all,string type,int avg,
								string errSymb,int useOrg,int nrBins);

#endif
