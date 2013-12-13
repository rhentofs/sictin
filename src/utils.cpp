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

#include "utils.h"
#include "bed_utils.h"


void Tokenize(const string& str,
	      vector<string>& tokens,
	      const string& delimiters)
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
}


/*
 * returns the 'next' position to start from.
 *
 */ 
string::size_type readAndTokenize(const string& str,
				  vector<string>& tokens,
				  const string& delimiters,
				  string::size_type firstPos)
{
  
  // Skip delimiters at beginning.
  firstPos = str.find_first_not_of(delimiters + "\n",firstPos);
  // Set the end.
  string::size_type lastPos  = str.find_first_of("\n", firstPos);
  if(lastPos == string::npos)
    return(str.length() + 1);
  string::size_type pos     = str.find_first_of(delimiters + "\n", firstPos+1);
  if(pos == string::npos)
    return(str.length() + 1);
  //cerr<<"RAT: "<<firstPos<<" "<<lastPos<<" "<<pos<<endl;
  while (pos < lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(firstPos, pos - firstPos));
      //cerr<<str.substr(firstPos, pos - firstPos)<<" ";
      // Skip delimiters.  Note the "not_of"
      firstPos = str.find_first_not_of(delimiters + "\n", pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters + "\n", firstPos);
    } 
  // Push the last one.  
  tokens.push_back(str.substr(firstPos, pos - firstPos));
  //cerr<<str.substr(firstPos, pos - firstPos)<<endl;
  //cout<<"RAT: "<<tokens[1]<<" "<<tokens[2]<<" "<<tokens[3]<<" "<<tokens[4]<<endl;
  return(lastPos);
}

/************************************************************
 * 
 * Random stuff. this should be uniform in the given span...
 * 
 *
 ************************************************************/

int uniform_range(int low,int high)
{
  int range = high - low + 1;
  int ret = low + (int)((double)range * rand()/(RAND_MAX + 1.0));
  return (ret);
}


/************************************************************
 *
 * Binary file accessors. 
 * 
 ************************************************************/


/************************************************************
 *
 * Calculate min & max position in a given file.
 *
 * writes parameters to the given int-pointers.
 * returns 0, -1 on error. 
 * 
 ************************************************************/

int getMinMaxC(const char* fname, int* min, int* max)
{
  // open the file and get some initial parameters. 
  ifstream inf;
  inf.open(fname,ios::binary| ios::ate);
  if(!inf.is_open())
    return(-1);
  *max  = ((int)inf.tellg() - (int)sizeof(int))/(int)sizeof(storageType);
  inf.seekg(0,ios::beg);
  inf.read((char*)min,sizeof(int));
  *max += (*min-1);
  inf.close();
  return(0);
}

// moves around the file-pointer some, disregards the effects. 
//(could be done better with ios::cur and stuff) 
int getMinMaxF(ifstream * inf, int* min, int* max)
{
  inf->seekg(0,ios::end);
  *max  = ((int)inf->tellg() - (int)sizeof(int))/(int)sizeof(storageType);
  inf->seekg(0,ios::beg);
  inf->read((char*)min,sizeof(int));
  *max += (*min-1);
  return(0);
}

/************************************************************
 *
 * Access the given file and retrieve the signal in the given
 * coordinates. 
 *
 * Fails (returns -1) if the coordinates are illegal.   
 *
 *
 ************************************************************/

int getSignalC(const char* fname, int from, int to, storageType *signal)
{
  int min,max;
  if(getMinMaxC(fname,&min,&max) < 0)
  {
	  cerr<<"getMinMaxC < 0"<<endl;
    return(-1);
  }
  if(min>from)
  {
	  cerr<<"min "<<min<<" > "<<" from "<<from<<endl;
    return(-1);
  }
  if(max<to)
  {
	  cerr<<"max "<<max<<" < "<<" to "<<to<<endl;
    return(-1);
  }
  ifstream inf;
  inf.open(fname,ios::binary);
  inf.seekg(sizeof(storageType)*(from - min)+sizeof(int),ios::beg);
  inf.read((char*)signal,sizeof(storageType)*(to-from +1));
  inf.close();
  return(0);
}

int getSignalF(ifstream * inf, int from, int to, storageType * signal)
{
  int min,max;
  if(getMinMaxF(inf,&min,&max) < 0)
    return(-1);
  if(min>from)
    return(-1);
  if(max<to)
    return(-1);
  inf->seekg(sizeof(storageType)*(from - min)+sizeof(int),ios::beg);
  inf->read((char*)signal,sizeof(storageType)*(to-from +1));
  return(0);
}

/************************************************************
 *
 * Access a given binary file, retrives all signals from a 
 * text-file with two columns start & stop. 
 * writes all failed access attempts to a second file. 
 *
 * puts the successful reads in another binary file, 
 * assumes all are of equal length (hard to sort out that file 
 * otherwise.)
 * 
 * returns the number of successful reads.
 *
 ************************************************************/

int getSignalsBin(const char * bin_file,const char* q_file, const char* out_file,const char* err_file)
{
  //open the binaries
  ifstream inf;
  ofstream outf;
  inf.open(bin_file,ios::binary);
  outf.open(out_file,ios::binary);

  // open the text-files
  ifstream qinf;
  ofstream errf;
  qinf.open(q_file);
  errf.open(err_file);
  
  // do some error checking. 
  if(!inf.is_open() || !outf.is_open() || !qinf.is_open()|| !errf.is_open())
    {
      cerr<<"Some files could not be open"<<endl;
      inf.close();
      outf.close();
      qinf.close();
      errf.close();
    }
  vector<string> data;
  string line;
  storageType *tmpSignals;
  int from,to,success = 0;
  getline(qinf,line);
  while(!qinf.eof())
    {
      data.clear();
      Tokenize (line,data," \t");
      from = atoi(data[0].c_str());
      to   = atoi(data[1].c_str());
      tmpSignals = new storageType[to - from + 1];
      if(getSignalF(&inf, from,to,tmpSignals)<0)
	{
	  errf<<data[0]<<" "<<data[0]<<endl;
	}
      else
	{ 
	  outf.write((char*)tmpSignals,sizeof(storageType)*(to - from + 1));
	  success++;
	}
      delete[] tmpSignals;
      getline(qinf,line); 
    }
  if(success == 0) // no queries worked.
    success = -1;
  inf.close();
  outf.close();
  qinf.close();
  errf.close();
  return(success);
}

int getSignalsTxt(const char* bin_file,const char* q_file, const char* out_file,const char* err_file,int avg,const char* err_symb)
{
  //open the binary
  ifstream inf;
  inf.open(bin_file,ios::binary);
  
  // open the text-files
  ifstream qinf;
  ofstream errf;
  ofstream outf;
  qinf.open(q_file);
  errf.open(err_file);
  outf.open(out_file);
   
  // do some error checking. 
  if(!inf.is_open() || !outf.is_open() || !qinf.is_open()|| !errf.is_open())
    {
      cerr<<"Some files could not be open"<<endl;
      inf.close();
      outf.close();
      qinf.close();
      errf.close();
    }
  vector<string> data;
  string line;
  storageType *tmpSignals;
  int from,to,success = 0;
  getline(qinf,line);
  while(!qinf.eof())
    {
      data.clear();
      Tokenize (line,data," \t");
      from = atoi(data[0].c_str());
      to   = atoi(data[1].c_str());
      
      // not a valid range.
      if((to - from + 1) < 1)
	{
	  errf<<data[0]<<"\t"<<data[1]<<endl;
	  if (strcmp(err_symb,"") != 0)
	    {
	      outf<<err_symb<<endl;
	    }
	}
      else
	{
	  tmpSignals = new storageType[to - from + 1];
	  if(getSignalF(&inf, from,to,tmpSignals)<0)
	    {
	      errf<<data[0]<<"\t"<<data[1]<<endl;
	      if (strcmp(err_symb,"") != 0)
		{
		  outf<<err_symb<<endl;
		}
	    }
	  else
	    { 
	      if (strcmp(err_symb,"") == 0)
		outf<<">"<<data[0]<<" "<<data[1]<<endl;
	      if(avg ==0)
		{
		  for (int i = 0; i <= to-from ;i++)
		    outf<<tmpSignals[i]<<" ";
		  outf<<endl;
		}
	      else
		{
		  double tmpsum = 0;
		  for (int i = 0; i <= to-from ;i++)
		    tmpsum += tmpSignals[i];
		  if((to-from+1) > 0)
		    tmpsum = (double)tmpsum/(to-from+1);
		  outf<<tmpsum<<endl;
		}
	      success ++;
	    }
	  delete[] tmpSignals;
	}	
      getline(qinf,line);
    }
  if(success == 0) // no queries worked.
    success = -1;
 
  inf.close();
  outf.close();
  qinf.close();
  errf.close();
  return(success);
}

int getSignalsFP(const char* bin_file,const char* fc_bin_file,const char* q_file, 
		 const char* out_file,const char* err_file,
		 int offset, int progress, int co, int move)
{
  //open the binaries
  ifstream inf;
  inf.open(bin_file,ios::binary);
  
  bool fc = (strcmp(fc_bin_file,"") != 0);
  ifstream fcinf;
  if(fc)
    fcinf.open(fc_bin_file,ios::binary);
  
  // open the text-files
  ifstream qinf;
  ofstream errf;
  ofstream outf;
  qinf.open(q_file);
  errf.open(err_file);
  outf.open(out_file);
   
  // do some error checking. 
  if(!inf.is_open())
    {
      cerr<<"infile could not be open"<<endl;
      fcinf.close();
      inf.close();
      outf.close();
      qinf.close();
      errf.close();
      return(-1);
    } 
  if(fc && !fcinf.is_open())
    {
      cerr<<"fc infile could not be open"<<endl;
      fcinf.close();
      inf.close();
      outf.close();
      qinf.close();
      errf.close();
      return(-1);
    } 
  if(!outf.is_open())
    {
      cerr<<"outfile could not be open"<<endl;
      fcinf.close();
      inf.close();
      outf.close();
      qinf.close();
      errf.close();
      return(-1);
    } 
  if(!qinf.is_open())
    {
      cerr<<"queryfile could not be open"<<endl;
      fcinf.close();
      inf.close();
      outf.close();
      qinf.close();
      errf.close();
      return(-1);
    } 
  if(!errf.is_open())
    {
      cerr<<"errorfile could not be open"<<endl;
      fcinf.close();
      inf.close();
      outf.close();
      qinf.close();
      errf.close();
      return(-1);
    } 
  
  int nlines = -1;
  string line;
  if(progress != 0)
    {
      // count the number of lines.
      while(!qinf.eof())
	{
	  nlines++;
	  getline(qinf,line);
	}
      cerr<<"Input consists of "<<nlines<<" queries."<<endl;
      // reset the file-pointer to the beginning.
      qinf.clear();
      qinf.seekg(0,ios::beg);
    }
  vector<string> data;
  storageType * tmpSignals;
  tmpSignals = new storageType[2*offset + 1]; 
  storageType * tmpSignals_fc;
  tmpSignals_fc = new storageType[2*offset + 1]; 
  // use this for Positive oriented signals.
  //unsigned int SignalsP[2*offset + 1];
  double SignalsP[2*offset + 1];
  for (int i = 0; i < (2*offset + 1);i++)
    SignalsP[i] = 0;
  // use this for Negative oriented signals.
  //unsigned int SignalsN[2*offset + 1];
  double SignalsN[2*offset + 1];
  for (int i = 0; i < (2*offset + 1);i++)
    SignalsN[i] = 0;
  int successP = 0;
  int successN = 0;
  int line_cnt = 1;
  int pos;
  bool passCo = false;
  getline(qinf,line);
  while(!qinf.eof())
    {
      line_cnt++;
      if((progress != 0) && ((line_cnt % 1000) == 0)) 
	cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<100*line_cnt/nlines<<" % complete.\r";
      
      data.clear();
      Tokenize (line,data," \t");
      pos = atoi(data[0].c_str());
      if(getSignalF(&inf, pos-offset+move,pos + offset+move,tmpSignals)<0 || 
	 (fc ? getSignalF(&fcinf, pos-offset+move,pos + offset+move,tmpSignals_fc)<0 : false))
	{
	  errf<<data[0]<<"\t"<<data[1]<<endl;
	}
      else
	{
	  // include this signal only if some position is higher that 'co'
	  for (int i = 0; i < (2*offset + 1);i++)
	    {
	      passCo = tmpSignals[i] >= co;
	      if(passCo)
		break;
	    }
	  if(passCo)
	    {
	      if(data[1] == "+")
		{
		  for (int i = 0; i < (2*offset + 1);i++)
		    if(fc)
		      {
			SignalsP[i] += (1.0 + (double)tmpSignals[i])/(1.0 + (double)tmpSignals_fc[i]);
		      }
		    else
		      SignalsP[i] += tmpSignals[i];
		  successP++;
                }
	      else
		{ 
		  for (int i = 0; i < (2*offset + 1);i++)
		    if(fc)
		      {
			if(tmpSignals_fc[2*offset -i] > 0)
			  SignalsN[i] += (1.0 + (double)tmpSignals[2*offset - i])/(1.0 + (double)tmpSignals_fc[2*offset - i]);
		      }
		    else
		      SignalsN[i] += tmpSignals[2*offset - i];
		  successN++;
                }
	    }
	}
      getline(qinf,line); 
    }
  if(successP == 0 && successN == 0) // no queries worked.
    {
      successP = -1;
      successN = -1;
    }
  else
    {
      outf<<">+ "<<successP<<endl;
      if(successP != 0)
	{	      
	  outf<<(double)SignalsP[0]/successP;
	  for (int i = 1; i < (2*offset + 1);i++)
	    outf<<" "<<(double)SignalsP[i]/successP;
	  outf<<endl;
	}
      else
	{
	  outf<<0;
	  for (int i = 1; i < (2*offset + 1);i++)
	    outf<<" "<<0;
	  outf<<endl;
	}     
      outf<<">- "<<successN<<endl;
      if(successN != 0)
	{	      
	  outf<<(double)SignalsN[0]/successN;
	  for (int i = 1; i < (2*offset + 1);i++)
	    outf<<" "<<(double)SignalsN[i]/successN;
	  outf<<endl;
	}
      else
	{
	  outf<<0;
	  for (int i = 1; i < (2*offset + 1);i++)
	    outf<<" "<<0;
	  outf<<endl;
	}     
    }
  if(progress != 0)
    cerr<<" 100 % complete"<<endl;
  cerr<<successP + successN<<" queries succesfully executed."<<endl;
  delete[] tmpSignals; 
  delete[] tmpSignals_fc;
 inf.close();
  fcinf.close();
  outf.close();
  qinf.close();
  errf.close();
  return(successP + successN);
}

int parseQueryLine(string str,inputLine* res)
{
  vector<string> data;
  data.clear();
  Tokenize (str,data," \t");
  if(data.size() != 0)
    {
      // these are the only applicable fields.
      res->seq    = data[0];
      res->strand = (data[2]== "+" ? 1:-1);
      res->pos    = atoi(data[1].c_str());
    }
  return(data.size());
}
/*
 *
 * checks the query file, ie which sequences are to be queried and how many 
 * queries on each. 
 *
 */

  
int queryControl(ifstream * inf,map<string,seqStats> * seqmap,bool verbose)
{
  int nlines = 0;
  seqStats tmpStat;
  string line;
  inputLine tmp;
  pair<map<string,seqStats>::iterator,bool> ret;
  // initialize tmpStat
  tmpStat.truncR = 0;
  tmpStat.truncF = 0;
  tmpStat.truncC = 0;
  
  while(!(inf->eof()))
    { 
      if((nlines % 1000) == 0) 
	cerr<<"-\r";
      if((nlines % 2000) == 0) 
	cerr<<"/\r";
      if((nlines % 3000) == 0) 
	cerr<<"|\r";
      if((nlines % 4000) == 0) 
	cerr<<"\\\r";
      if((nlines % 5000) == 0) 
	cerr<<"-\r";
      
      getline(*inf,line); 
      if(parseQueryLine(line,&tmp) < 1) // skip possible empty lines 
	continue;
      
      nlines++;
      // try and insert into seqmap.
      tmpStat.minPos = tmp.pos;
      tmpStat.maxPos = tmp.pos; 
      tmpStat.countF = 0;
      tmpStat.countR = 0;
      if(tmp.strand == 1) // forward strand
	tmpStat.countF++;
      else
	tmpStat.countR++;
      ret = seqmap->insert(pair<string,seqStats> (tmp.seq,tmpStat));
      
      if(!ret.second) // already a key with this value.
	{
	  if(tmpStat.minPos < (*seqmap)[tmp.seq].minPos)
	    (*seqmap)[tmp.seq].minPos = tmpStat.minPos;
	  
	  if(tmpStat.maxPos > (*seqmap)[tmp.seq].maxPos)
	    (*seqmap)[tmp.seq].maxPos = tmpStat.maxPos;
	  
	  if(tmp.strand == 1) // forward strand
	    (*seqmap)[tmp.seq].countF++;
	  else
	    (*seqmap)[tmp.seq].countR++;
	}
    }
  // reset the file-stream  
  inf->clear(); 
  inf->seekg(0,ios::beg);
  return(nlines);
}

int queryControlBED(ifstream * inf,map<string,seqStats> * seqmap,bool verbose,string type)
{
  int nlines = 0;
  seqStats tmpStat;
  string line;
  inputLine tmp;
  pair<map<string,seqStats>::iterator,bool> ret;
  // initialize tmpStat
  tmpStat.truncR = 0;
  tmpStat.truncF = 0;
  tmpStat.truncC = 0;
  
  while(!(inf->eof()))
    { 
      if((nlines % 1000) == 0) 
	cerr<<"-\r";
      if((nlines % 2000) == 0) 
	cerr<<"/\r";
      if((nlines % 3000) == 0) 
	cerr<<"|\r";
      if((nlines % 4000) == 0) 
	cerr<<"\\\r";
      if((nlines % 5000) == 0) 
	cerr<<"-\r";
      
      getline(*inf,line); 
      if(parseBEDline(line,&tmp,0,1,2,5) < 1) // skip possible empty lines
	continue;
      // prepare the query accoding to which position that should be used?
      if(strcmp(type.c_str(),"c") == 0)
	tmp.pos += (int)(((double)(tmp.len-1)/2.0) + 0.5);
      if(strcmp(type.c_str(),"s") == 0)
	if(tmp.strand == -1) // negative strand. flip. 
	  tmp.pos += (tmp.len-1);
      if(strcmp(type.c_str(),"e") == 0)
	if(tmp.strand == 1) // positive strand
	  tmp.pos += (tmp.len-1);
      
      nlines++;
      // try and insert into seqmap.
      tmpStat.minPos = tmp.pos;
      tmpStat.maxPos = tmp.pos; 
      tmpStat.countF = 0;
      tmpStat.countR = 0;
      if(tmp.strand == 1) // forward strand
	tmpStat.countF++;
      else
	tmpStat.countR++;
      ret = seqmap->insert(pair<string,seqStats> (tmp.seq,tmpStat));
      
      if(!ret.second) // already a key with this value.
	{
	  if(tmpStat.minPos < (*seqmap)[tmp.seq].minPos)
	    (*seqmap)[tmp.seq].minPos = tmpStat.minPos;
	  
	  if(tmpStat.maxPos > (*seqmap)[tmp.seq].maxPos)
	    (*seqmap)[tmp.seq].maxPos = tmpStat.maxPos;
	  
	  if(tmp.strand == 1) // forward strand
	    (*seqmap)[tmp.seq].countF++;
	  else
	    (*seqmap)[tmp.seq].countR++;
	}
    }
  // reset the file-stream  
  inf->clear(); 
  inf->seekg(0,ios::beg);
  return(nlines);
}


int getMultipleSignalsFP(ifstream * infq,map<string,fsHolder> seqFiles,string ofStub,string eFile,
			 int offset,int progress,int co,int move,int both,int nq,int all)

{
  
  string preR = "R_";
  string preF = "F_";
  string preC = "C_";
  string preRall = "R_all_";
  string preFall = "F_all_";
  string preCall = "C_all_";

  bool allFilesOk = true;
  ofstream ofs[4];
  ofstream aofs[3];
  
  // the order of the output files is assumed to be FRC below. indexes (0-2) are used.
  ofs[0].open((preF + ofStub).c_str(),ios::trunc);
  if(!ofs[0].good())
    {
      cerr<<"Could not open "<<(preF + ofStub).c_str()<<endl;
      allFilesOk = false;
    }
  ofs[1].open((preR + ofStub).c_str(),ios::trunc);
  if(!ofs[1].good())
    {
      cerr<<"Could not open "<<(preR + ofStub).c_str()<<endl;
      allFilesOk = false;
    }
  
  if(both > 0)
    {  
      ofs[2].open((preC + ofStub).c_str(),ios::trunc);
      if(!ofs[2].good())
	{
	  cerr<<"Could not open "<<(preC + ofStub).c_str()<<endl;
	  allFilesOk = false;
	}
    }
  
  if(all == 1) // each data point should be written. open the extra output files. 
    {
      aofs[0].open((preFall + ofStub).c_str(),ios::trunc);
      if(!aofs[0].good())
	{
	  cerr<<"Could not open "<<(preFall + ofStub).c_str()<<endl;
	  allFilesOk = false;
	}
      aofs[1].open((preRall + ofStub).c_str(),ios::trunc);
      if(!aofs[1].good())
	{
	  cerr<<"Could not open "<<(preRall + ofStub).c_str()<<endl;
	  allFilesOk = false;
	}
      
      if(both > 0)
	{  
	  aofs[2].open((preCall + ofStub).c_str(),ios::trunc);
	  if(!aofs[2].good())
	    {
	      cerr<<"Could not open "<<(preCall + ofStub).c_str()<<endl;
	      allFilesOk = false;
	    }
	}
    }
  
  // the error-file, where failed queries will be recorded.
  ofs[3].open(eFile.c_str(),ios::trunc);
  if(!ofs[3].good())
    {
      cerr<<"Could not open "<<eFile<<endl;
      allFilesOk = false;
    }
  
  


  if(!allFilesOk)
    {
      cerr<<"All necessary output files could no be opened. Aborting."<<endl;
      for(int i = 0;i<4;i++)
	ofs[i].close();
      for(int i = 0;i<3;i++)
	aofs[i].close();
      return(-1);
    }

  /*
   *  process queries and store +/- signals for each type R,F and (possibly) C
   * 
   */ 
  int line_cnt = 0;
  bool passCo = false;
  string line;
  inputLine tmp;
  storageType * tmpSignals;
  tmpSignals = new storageType[2*offset + 1]; 
  // use this to store the signals, FP,RP,CP,FN,RN,CN ordered.
  unsigned int signals[6][2*offset + 1];
  // counters. FP,RP,CP,FN,RN,CN ordered.
  int succCnt[6] = {0,0,0,0,0,0};
  int sigReturn;
  for (int j = 0; j < 6;j++)
    for (int i = 0; i < (2*offset + 1);i++)
      signals[j][i] = 0;
  
  while(!infq->eof())
    {
      line_cnt++;
      if((progress != 0) && ((line_cnt % 1000) == 0)) 
	cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<100*line_cnt/nq<<" % complete.\r";
      // read and parse one query.
      getline(*infq,line); 
      if(parseQueryLine(line,&tmp) < 1) // skip possible empty lines 
	continue;
      for (int sigType = 0;sigType<=2;sigType++)
	{
	  if((sigType == 2) && (both == 0)) // combined but no combined wanted
	    break;
	  // Forward (0) 
	  if(sigType == 0)
	    sigReturn = getSignalF((ifstream*)seqFiles[tmp.seq].fsF,
				   tmp.pos-offset+move,tmp.pos+offset+move,tmpSignals);
	  // Reverse (1) 
	  if(sigType == 1)
	    sigReturn = getSignalF((ifstream*)seqFiles[tmp.seq].fsR,
				   tmp.pos-offset+move,tmp.pos+offset+move,tmpSignals);
	  // Combined (2) 
	  if(sigType == 2)
	    sigReturn = getSignalF((ifstream*)seqFiles[tmp.seq].fsC,
				   tmp.pos-offset+move,tmp.pos+offset+move,tmpSignals);
	  
	  if(sigReturn < 0) // failed query. write to error-file
	    {
	      ofs[3]<<tmp.seq<<"\t"<<tmp.pos<<"\t"<<(tmp.strand == 1 ? "+":"-")<<endl;
	    }
	  else
	    {
	      passCo = false;
	      if(co == 0) // are we to look at cut-offs?
		passCo = true;
	      else
		{
		  // run through until we find atleast one occurance above co
		  for (int i = 0; i < (2*offset + 1);i++)
		    {
		      passCo = tmpSignals[i] >= co;
		      if(passCo)
			break;
		    }
		}
	      if(passCo)
		{
		  if(tmp.strand == 1)
		    {
		      for (int i = 0; i < (2*offset + 1);i++)
			signals[0+sigType][i] += tmpSignals[i];
		      succCnt[0+sigType]++;
		      if(all == 1)
			{
			  for (int i = 0; i < (2*offset + 1);i++)
			    if(i == 0)
			      aofs[sigType]<<tmpSignals[i];
			    else
			      aofs[sigType]<<"\t"<<tmpSignals[i];
			  aofs[sigType]<<endl;
			}
		    }
		  else
		    { 
		      for (int i = 0; i < (2*offset + 1);i++)
			signals[3+sigType][i] += tmpSignals[2*offset - i];
		      succCnt[3+sigType]++;
		      if(all == 1)
			{
			  for (int i = 0; i < (2*offset + 1);i++)
			    if(i == 0)
			      aofs[sigType]<<tmpSignals[2*offset - i];
			    else
			      aofs[sigType]<<"\t"<<tmpSignals[2*offset - i];
			  aofs[sigType]<<endl;
			}
		    }
		}
	      
	    }
	}
    } // while()
  // clear the infq in case the caller expects it to be ok.
  infq->clear();
  
  /*
   * write out results.
   */
  for (int sigType = 0;sigType<=2;sigType++)
    {
      if((sigType == 2) && (both == 0)) // combined but no combined wanted
	break;
      // sense strand
      ofs[sigType]<<">+ "<<succCnt[0+sigType]<<endl;
      ofs[sigType]<<(double)signals[sigType][0]/(succCnt[0+sigType] == 0 ? 1 : succCnt[0+sigType]);
      for (int i = 1; i < (2*offset + 1);i++)
	ofs[sigType]<<" "<<(double)signals[sigType][i]/(succCnt[0+sigType] == 0 ? 1 : succCnt[0+sigType]);
      ofs[sigType]<<endl;
      // anti sense strand
      ofs[sigType]<<">- "<<succCnt[3+sigType]<<endl;
      ofs[sigType]<<(double)signals[3+sigType][0]/(succCnt[3+sigType] == 0 ? 1 : succCnt[3+sigType]);
      for (int i = 1; i < (2*offset + 1);i++)
	ofs[sigType]<<" "<<(double)signals[3+sigType][i]/(succCnt[3+sigType] == 0 ? 1 : succCnt[3+sigType]);
      ofs[sigType]<<endl;
    }
  
  for(int i = 0;i<4;i++)
    ofs[i].close();
  for(int i = 0;i<3;i++)
    aofs[i].close();
  return(1);
}



int getMultipleSignalsFPfromBED(ifstream * infq,map<string,fsHolder> seqFiles,string ofStub,string eFile,
				int offset,int progress,int co,int move,int both,int nq,int all,string type,
								int avg,string errSymb,int useOrg,int nrBins)

{
  
  string preR = "R_";
  string preF = "F_";
  string preC = "C_";
  string preRall = "R_all_";
  string preFall = "F_all_";
  string preCall = "C_all_";

  bool allFilesOk = true;
  ofstream ofs[4];
  ofstream aofs[3];
  
  // the order of the output files is assumed to be FRC below. indexes (0-2) are used.
  ofs[0].open((preF + ofStub).c_str(),ios::trunc);
  if(!ofs[0].good())
    {
      cerr<<"Could not open "<<(preF + ofStub).c_str()<<endl;
      allFilesOk = false;
    }
  ofs[1].open((preR + ofStub).c_str(),ios::trunc);
  if(!ofs[1].good())
    {
      cerr<<"Could not open "<<(preR + ofStub).c_str()<<endl;
      allFilesOk = false;
    }
  
  if(both > 0)
    {  
      ofs[2].open((preC + ofStub).c_str(),ios::trunc);
      if(!ofs[2].good())
		{
		  cerr<<"Could not open "<<(preC + ofStub).c_str()<<endl;
		  allFilesOk = false;
		}
    }
  
  if(all == 1) // each data point should be written. open the extra output files. 
    {
      aofs[0].open((preFall + ofStub).c_str(),ios::trunc);
      if(!aofs[0].good())
		{
		  cerr<<"Could not open "<<(preFall + ofStub).c_str()<<endl;
		  allFilesOk = false;
		}
      aofs[1].open((preRall + ofStub).c_str(),ios::trunc);
      if(!aofs[1].good())
		{
		  cerr<<"Could not open "<<(preRall + ofStub).c_str()<<endl;
		  allFilesOk = false;
		}
      
      if(both > 0)
		{  
		  aofs[2].open((preCall + ofStub).c_str(),ios::trunc);
		  if(!aofs[2].good())
			{
			  cerr<<"Could not open "<<(preCall + ofStub).c_str()<<endl;
			  allFilesOk = false;
			}
		}
    }
  
  // the error-file, where failed queries will be recorded.
  ofs[3].open(eFile.c_str(),ios::trunc);
  if(!ofs[3].good())
    {
      cerr<<"Could not open "<<eFile<<endl;
      allFilesOk = false;
    }
  
  


  if(!allFilesOk)
    {
      cerr<<"All necessary output files could no be opened. Aborting."<<endl;
      for(int i = 0;i<4;i++)
		ofs[i].close();
      for(int i = 0;i<3;i++)
		aofs[i].close();
	  return(-1);
    }

  /*
   *  process queries and store +/- signals for each type R,F and (possibly) C
   * 
   */ 
  int line_cnt = 0;
  bool passCo = false;
  string line;
  inputLine tmp;
  storageType * tmpSignals;
  tmpSignals = new storageType[2*offset + 1]; 
  storageType * tmpSignals2; // to be used if 'useOrg' is set.
  storageType * tmpSignals3; // to be used if 'nrBins' is set.
  int orgStart,orgEnd; 
  // use this to store the signals, FP,RP,CP,FN,RN,CN ordered.
  unsigned int signals[6][2*offset + 1];
  // counters. FP,RP,CP,FN,RN,CN ordered.
  int succCnt[6] = {0,0,0,0,0,0};
  int sigReturn,sigReturn2,sigReturn3;
  for (int j = 0; j < 6;j++)
    for (int i = 0; i < (2*offset + 1);i++)
      signals[j][i] = 0;
  // for the binned results.
  int binSize;
  int binCntAdd;
  int binSum;
  int binCnt;
  int binCurPos;
  double * binValues[6];
  for(int i = 0;i<6;i++)
	{
	  binValues[i] = new double[nrBins + 2*offset]; // used to store the averages and initial/ending values
	  for(int j = 0;j<(nrBins + 2*offset);j++)
		binValues[i][j] = 0.0;
	}
  int nrBin;
  while(!infq->eof())
    {
      line_cnt++;
      if((progress != 0) && ((line_cnt % 1000) == 0)) 
		cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<100*line_cnt/nq<<" % complete.\r";
      // read and parse one query.
      getline(*infq,line); 
      if(parseBEDline(line,&tmp,0,1,2,5) < 1) // skip possible empty lines
		continue;
      // prepare the query accoding to which position that should be used?
      if(useOrg == 1 || nrBins > 0) // original region needed
		{
		  orgStart = tmp.pos;
		  orgEnd   = tmp.pos + tmp.len -1;
		}
	  
      if(strcmp(type.c_str(),"c") == 0)
		tmp.pos += (int)(((double)(tmp.len-1)/2.0) + 0.5);
      if(strcmp(type.c_str(),"s") == 0)
		if(tmp.strand == -1) // negative strand. flip. 
		  tmp.pos += (tmp.len-1);
      if(strcmp(type.c_str(),"e") == 0)
		if(tmp.strand == 1) // positive strand
		  tmp.pos += (tmp.len-1);
	  
      for (int sigType = 0;sigType<=2;sigType++)
		{
		  if((sigType == 2) && (both == 0)) // combined but no combined wanted
			break;
		  if(nrBins == 0) // the "normal" case 
			{
			  // Forward (0) 
			  if(sigType == 0)
				sigReturn = getSignalF((ifstream*)seqFiles[tmp.seq].fsF,
									   tmp.pos-offset+move,tmp.pos+offset+move,tmpSignals);
			  // Reverse (1) 
			  if(sigType == 1)
				sigReturn = getSignalF((ifstream*)seqFiles[tmp.seq].fsR,
									   tmp.pos-offset+move,tmp.pos+offset+move,tmpSignals);
			  // Combined (2) 
			  if(sigType == 2)
				sigReturn = getSignalF((ifstream*)seqFiles[tmp.seq].fsC,
									   tmp.pos-offset+move,tmp.pos+offset+move,tmpSignals);
			}
		  if(all == 1 && useOrg == 1) // write out orignal regions in every position
			{
			  tmpSignals2 = new storageType[tmp.len]; 
			  
			  // Forward (0) 
			  if(sigType == 0)
				sigReturn2 = getSignalF((ifstream*)seqFiles[tmp.seq].fsF,
										orgStart+move,orgEnd+move,tmpSignals2);
			  // Reverse (1) 
			  if(sigType == 1)
				sigReturn2 = getSignalF((ifstream*)seqFiles[tmp.seq].fsR,
										orgStart+move,orgEnd+move,tmpSignals2);
			  // Combined (2) 
			  if(sigType == 2)
				sigReturn2 = getSignalF((ifstream*)seqFiles[tmp.seq].fsC,
										orgStart+move,orgEnd+move,tmpSignals2);
			}
		  
		  if(nrBins > 0) // here we need original regions + offset
			{
			  tmpSignals3 = new storageType[tmp.len + 2*offset]; 
			
			  // Forward (0) 
			  if(sigType == 0)
				sigReturn3 = getSignalF((ifstream*)seqFiles[tmp.seq].fsF,
										orgStart+move-offset,orgEnd+move+offset,tmpSignals3);
			  // Reverse (1) 
			  if(sigType == 1)
				sigReturn3 = getSignalF((ifstream*)seqFiles[tmp.seq].fsR,
										orgStart+move-offset,orgEnd+move+offset,tmpSignals3);
			  // Combined (2) 
			  if(sigType == 2)
				sigReturn3 = getSignalF((ifstream*)seqFiles[tmp.seq].fsC,
										orgStart+move-offset,orgEnd+move+offset,tmpSignals3);
			}
		  
		  
		  // this whole part needs to be split depending on if 'useOrg' or nrBins' are set, the offset-queries or the original could fail independently. 
		  if(nrBins == 0)
			{
			  if(sigReturn < 0) // failed query. write to error-file
				{
				  ofs[3]<<tmp.seq<<"\t"<<tmp.pos<<"\t"<<(tmp.strand == 1 ? "+":"-")<<endl;
				  if(all == 1 && useOrg == 0) 
					{
					  if(avg == 1)
						{
						  aofs[sigType]<<errSymb<<endl; // write error-symbol to "all"-files if '-a' and '-avg' were set. 
						}
					  else
						{
						  for (int i = 0; i < (2*offset + 1);i++) // write (matrix preserving) error-symbols to "all"-files if '-avg' is not set. 
							if(i == 0)
							  aofs[sigType]<<errSymb;
							else
							  aofs[sigType]<<"\t"<<errSymb;
						  aofs[sigType]<<endl;
						}
					}
				}
			  else // signals retrieved ok.
				{
				  passCo = false;
				  if(co == 0) // are we to look at cut-offs?
					passCo = true;
				  else
					{
					  // run through until we find atleast one occurance above co
					  for (int i = 0; i < (2*offset + 1);i++)
						{
						  passCo = tmpSignals[i] >= co;
						  if(passCo)
							break;
						}
					}
				  if(passCo)
					{
					  for (int i = 0; i < (2*offset + 1);i++)
						signals[(tmp.strand == 1 ? 0 : 3)+sigType][i] += tmpSignals[(tmp.strand == 1 ? i : 2*offset-i)];
					  succCnt[(tmp.strand == 1 ? 0 : 3)+sigType]++;
					  if(all == 1 && avg  == 0 && useOrg == 0)
						{
						  for (int i = 0; i < (2*offset + 1);i++)
							if(i == 0)
							  aofs[sigType]<<tmpSignals[(tmp.strand == 1 ? i : 2*offset-i)];
							else
							  aofs[sigType]<<"\t"<<tmpSignals[(tmp.strand == 1 ? i : 2*offset-i)];
						  aofs[sigType]<<endl;
						}
					}
				  if(all == 1 && avg == 1 && useOrg == 0) // no point in splitting on direction (+/-) when averaging. 
					{
					  double tmpSum = 0.0;
					  for (int i = 0; i < (2*offset + 1);i++)
						tmpSum+=(double)tmpSignals[i];
					  tmpSum = tmpSum/(2.0*(double)offset + 1.0);
					  aofs[sigType]<<tmpSum<<endl;
					}
				}
			  /*
			   * this whole part is the above repeated but for the special case when useOrg is set. 
			   */
			  if(all == 1 && useOrg == 1) // if useOrg is set, nrBins cannot be set since nrBins auto turns useOrg off. right. 
				{
				  if(sigReturn2 < 0) // failed query. write to error-file
					aofs[sigType]<<errSymb<<endl;
				  else
					{
					  passCo = false;
					  if(co == 0) // are we to look at cut-offs?
						passCo = true;
					  else
						{
						  // run through until we find atleast one occurance above co
						  for (int i = 0; i < tmp.len;i++)
							{
							  passCo = tmpSignals2[i] >= co;
							  if(passCo)
								break;
							}
						}
					  if(passCo)
						{
						  if(avg == 0)
							{
							  for (int i = 0; i < tmp.len;i++)
								if(i == 0)
								  aofs[sigType]<<tmpSignals2[(tmp.strand == 1 ? i : tmp.len -i -1)];
								else
								  aofs[sigType]<<"\t"<<tmpSignals2[(tmp.strand == 1 ? i : tmp.len -i -1)];
							  aofs[sigType]<<endl;
							}
						  
						  if(avg == 1) // no point in splitting on direction (+/-) when averaging. 
							{
							  double tmpSum = 0.0;
							  for (int i = 0; i < tmp.len;i++)
								tmpSum+=(double)tmpSignals2[i];
							  tmpSum = tmpSum/(double)tmp.len;
							  aofs[sigType]<<tmpSum<<endl;
							}
						  
						}
					}
				  delete[] tmpSignals2;
				} // all & useOrg
			}
		  /*
		   * binned results
		   */
		  if(nrBins > 0) // if useOrg is set, nrBins cannot be set since nrBins auto turns useOrg off. right. 
			{
			  if(sigReturn3 < 0) // failed query. write to error-file
				{
				  ofs[3]<<tmp.seq<<"\t"<<tmp.pos<<"\t"<<(tmp.strand == 1 ? "+":"-")<<endl;
				  if(all > 0)
					{
					  for(int i = 0;i<(nrBins + 2*offset);i++) // write out 'errSymb' at each position to facilitate matrix processing later on. 
						if(i == 0)
						  aofs[sigType]<<errSymb;
						else
						  aofs[sigType]<<"\t"<<errSymb;
					  aofs[sigType]<<endl;
					}
				}
			  else
				{
				  passCo = false;
				  if(co == 0) // are we to look at cut-offs?
					passCo = true;
				  else
					{
					  // run through until we find atleast one occurance above co
					  for (int i = 0; i < (tmp.len+2*offset);i++)
						{
						  passCo = tmpSignals3[i] >= co;
						  if(passCo)
							break;
						}
					}
				  if(passCo) 
					{
					  if (tmp.len < nrBins) // too short region that should be binned
						{
						  ofs[3]<<tmp.seq<<"\t"<<tmp.pos<<"\t"<<(tmp.strand == 1 ? "+":"-")<<endl;
						  if(all == 1)
							aofs[sigType]<<errSymb<<endl;
						}
					  else
						{
						  succCnt[(tmp.strand == 1 ? 0 : 3)+sigType]++;
						  for(int i = 0;i<offset;i++)
							{
							  binValues[(tmp.strand == 1 ? 0 : 3)+sigType][i] += tmpSignals3[(tmp.strand == 1 ? i : tmp.len+2*offset-i-1)];
							  if(all==1)
								aofs[sigType]<<tmpSignals3[(tmp.strand == 1 ? i : tmp.len + 2*offset -i -1)]<<"\t";
							}

						  // decide on binsize, loop over the regions and calculate average signals,write out if 'all' is set.
						  binSize = tmp.len/nrBins;
						  binCntAdd = tmp.len - nrBins*binSize; // these many bins should have 1bp extra in length.
						  
						  binCurPos = offset;
						  nrBin = 0;
						  while(binCurPos < (tmp.len+offset))
							{
							  //cerr<<"["<<binCurPos<<","<<(binCurPos + binSize + (binCntAdd > 0 ? 1 : 0))<<") "<<tmp.len<<endl;
							  binSum = 0;
							  binCnt = 0;
							  for(int i = binCurPos; i < (binCurPos + binSize + (binCntAdd > 0 ? 1 : 0));i++)
								{
								  binSum += tmpSignals3[(tmp.strand == 1 ? i : tmp.len+2*offset -i -1)];
								  binCnt++;
								}

							  binCntAdd--;
							  if(all == 1)
								if(binCurPos == offset)
								  aofs[sigType]<<(double)(binSum)/(double)(binCnt);
								else
								  aofs[sigType]<<"\t"<<(double)(binSum)/(double)(binCnt);
							  binCurPos += binCnt;
							  binValues[(tmp.strand == 1 ? 0 : 3)+sigType][offset + nrBin] += (double)(binSum)/(double)(binCnt); 
							  nrBin++;
							}
						  
						  // get the last bases after the bins.
						  for(int i = 0;i<offset;i++)
							{
							  binValues[(tmp.strand == 1 ? 0 : 3)+sigType][nrBin+offset+i] += (double)tmpSignals3[(tmp.strand == 1 ?offset+tmp.len+i : offset-i-1)];
							  if(all == 1)
								aofs[sigType]<<"\t"<<tmpSignals3[(tmp.strand == 1 ? offset+tmp.len+i : offset-i-1)];
							}
						  if(all == 1)
							aofs[sigType]<<endl;
						}
					  /*
						for(int i = 0;i<6;i++)
						{
						for(int j = 0;j<(nrBins + 2*offset);j++)
						cerr<<binValues[i][j]<<" ";
						cerr<<endl;
						}
					  */
					}
				}
			  delete[] tmpSignals3;
			}
		} // sigType
    } // while()
  // clear the infq in case the caller expects it to be ok.
  infq->clear();
  delete[] tmpSignals;
  
  /*
   * write out results.
   */
  for (int sigType = 0;sigType<=2;sigType++)
    {
      if((sigType == 2) && (both == 0)) // combined but no combined wanted
		break;
      // sense strand
      ofs[sigType]<<">+ "<<succCnt[0+sigType]<<endl;
      if(nrBins == 0)
		{
		  ofs[sigType]<<(double)signals[sigType][0]/(succCnt[0+sigType] == 0 ? 1 : succCnt[0+sigType]);
		  for (int i = 1; i < (2*offset + 1);i++)
			ofs[sigType]<<" "<<(double)signals[sigType][i]/(succCnt[0+sigType] == 0 ? 1 : succCnt[0+sigType]);
		}
	  else
		{
		  ofs[sigType]<<(double)binValues[sigType][0]/(succCnt[0+sigType] == 0 ? 1 : succCnt[0+sigType]);
		  for (int i = 1; i < (2*offset + nrBins);i++)
			ofs[sigType]<<" "<<(double)binValues[0+sigType][i]/(succCnt[0+sigType] == 0 ? 1 : succCnt[0+sigType]);
		}
	  ofs[sigType]<<endl;
      // anti sense strand
	  ofs[sigType]<<">- "<<succCnt[3+sigType]<<endl;
      if(nrBins == 0)
		{
		  ofs[sigType]<<(double)signals[3+sigType][0]/(succCnt[3+sigType] == 0 ? 1 : succCnt[3+sigType]);
		  for (int i = 1; i < (2*offset + 1);i++)
			ofs[sigType]<<" "<<(double)signals[3+sigType][i]/(succCnt[3+sigType] == 0 ? 1 : succCnt[3+sigType]);
		}
	  else
		{
		  ofs[sigType]<<(double)binValues[3+sigType][0]/(succCnt[3+sigType] == 0 ? 1 : succCnt[3+sigType]);
		  for (int i = 1; i < (2*offset + nrBins);i++)
			ofs[sigType]<<" "<<(double)binValues[3+sigType][i]/(succCnt[3+sigType] == 0 ? 1 : succCnt[3+sigType]);
		}
	  ofs[sigType]<<endl;
    }
  
  for(int i = 0;i<6;i++)
	delete binValues[i];
  for(int i = 0;i<4;i++)
    ofs[i].close();
  for(int i = 0;i<3;i++)
    aofs[i].close();
  return(1);
}
