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

#include "sam_utils.h"

int parseSAMline(string line,inputLine * res)
{
  vector<string> data;
  data.clear();
  /*
   * 2 FLAG 
   * 3 RNAME 
   * 4 POS, 1-based, leftmost 
   * 10 SEQ
   */
  Tokenize (line,data,"\t");
  if(data.size() == 0)
    return(data.size());
  res->header = 0;
  if((data[0]).substr(0,1) == "@") // header line.
    {
      res->header = 1;
      return(data.size());
    }
  // the already prepared fields 
  res->seq = data[2];
  res->pos = atoi((data[3]).c_str());
  res->len = (data[9]).length();
  // parse the FLAG
  unsigned int tmpFlag = atoi((data[1]).c_str());
  res->strand = 1;
  if((tmpFlag & SAM_REV_STRAND) == SAM_REV_STRAND)
    res->strand = -1;
  res->mapped = 1;
  if((tmpFlag & SAM_UNMAPPED_SEQ) == SAM_UNMAPPED_SEQ)
    {
      res->mapped = 0;
    }
  res->val = (storageType)1; // always 1 for SAM-files.
  return(data.size());
}

/*
 * returns the number of lines with mapped reads 
 */
int initControlSAM(ifstream * inf,map<string,seqStats> * seqmap,int both,bool verbose)
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
      if(parseSAMline(line,&tmp) < 1) // skip possible empty lines 
	continue;
      if(!tmp.header)
	{
	  if(tmp.mapped == 1) // query is mapped
	    {
	      //cout<<"data: "<<tmp.strand<<" "<<tmp.pos<<" "<<tmp.len<<" "<<tmp.seq<<" "<<tmp.mapped<<endl;
	      nlines++;
	      // try and insert into seqmap.
	      tmpStat.minPos = max(tmp.pos-both,0); // SAM supposed to be 1-based. but allow zero anyway. 
	      tmpStat.maxPos = tmp.pos + tmp.len-1 + both; // this will be larger that needed, but ok.
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
	}
    }
  
  // reset the file-stream  
  inf->clear(); 
  inf->seekg(0,ios::beg);
  return(nlines);
}

