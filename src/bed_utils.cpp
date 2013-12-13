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

#include "bed_utils.h"

int parseBEDline(string line,inputLine * res,int seqCol,int startCol,
				 int endCol,int strandCol, int nameCol)
{
  vector<string> data;
  float tmpVal = 1.0,tmp; 
  data.clear();
  Tokenize (line,data," \t"); // need to split on both " " and "\t" to be able to parse headers.
  if(data.size() == 0)
    return(data.size());
  res->header = 0;    
  if(strcmp((data[0]).c_str(),"track") == 0) // header line, should start with "track name=..."
    {
      res->header = 1;
      return(data.size());
    }
  res->mapped = 1;     // only mapped fragments
  res->seq    = data[seqCol];
  res->strand = (data[strandCol]== "+" ? 1:-1);
  res->pos    = atoi(data[startCol].c_str()) +1; // BED format is 0-based. change to 1-based.
  res->len    = atoi(data[endCol].c_str()) - atoi(data[startCol].c_str()); // crds are [x,y) so no need for +1 in length.
  res->name   = data[nameCol];
  res->val = (storageType)tmpVal;
  return(data.size());
}

/*
 * returns the number of lines with mapped reads 
 */
int initControlBED(ifstream * inf,map<string,seqStats> * seqmap,int seqCol,int startCol,
		   int endCol,int strandCol,int both,bool verbose)
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
      if(parseBEDline(line,&tmp,seqCol,startCol,endCol,strandCol) < 1) // skip possible empty lines 
	continue;
      
      if(tmp.header)
	continue;

      nlines++;
      // try and insert into seqmap.
      tmpStat.minPos = max(tmp.pos-both,0); // allow 0-based coordinates. 
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

