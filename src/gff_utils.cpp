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

#include "gff_utils.h"

int parseGFFline(string line,inputLine * res,string delim,int seqCol,int startCol,
		 int endCol,int strandCol,int valCol)
{
  vector<string> data;
  float tmpVal,tmp; 
  data.clear();
  Tokenize (line,data,delim);
  if(data.size() == 0)
    return(data.size());
  res->header = 0;     // no headers
  res->mapped = 1;     // only mapped fragments
  res->seq    = data[seqCol];
  res->strand = (data[strandCol]== "+" ? 1:-1);
  res->pos    = atoi(data[startCol].c_str());
  res->len    = atoi(data[endCol].c_str()) - atoi(data[startCol].c_str()) + 1;
  
  if(valCol != -1)
    {
      tmp = atof(data[valCol].c_str());
      if(fabs(tmp)>MAXNR)
	tmp=(tmp>0?1:-1)*(double)MAXNR;
      tmpVal=(storageType)tmp;
    }
  else
    {
      tmpVal = 1.0;
    }
  res->val = (storageType)tmpVal;
  return(data.size());
}

/*
 * returns the number of lines with mapped reads 
 */
int initControlGFF(ifstream * inf,map<string,seqStats> * seqmap,string delim,int seqCol,int startCol,
		   int endCol,int strandCol,int valCol,int both,bool verbose)
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
      if((nlines % 10000) == 0) 
	cerr<<"-\r";
      if((nlines % 20000) == 0) 
	cerr<<"/\r";
      if((nlines % 30000) == 0) 
	cerr<<"|\r";
      if((nlines % 40000) == 0) 
	cerr<<"\\\r";
      if((nlines % 50000) == 0) 
	cerr<<"-\r";
      
      getline(*inf,line); 
      if(parseGFFline(line,&tmp,delim,seqCol,startCol,endCol,strandCol,valCol) < 1) // skip possible empty lines 
	continue;

      nlines++;
      // try and insert into seqmap.
      tmpStat.minPos = max(tmp.pos-both,0); // allow zero-based coordinates. 
      tmpStat.maxPos = tmp.pos + tmp.len-1 + both; 
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

