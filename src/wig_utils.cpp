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

#include "wig_utils.h"

int parseWIGline(string line,inputLine * res,string* currName,int* currStep,int* currPos,int* currSpan,
		 int* fixed, int scale)
{
  vector<string> data,detailData;
  data.clear();
  detailData.clear();
  Tokenize (line,data," \t\r");
  
  if(data.size() == 0)
    return(data.size());
  res->header = 0;
  if(strcmp(data[0].c_str(),"track") ==0) // track line, ignore
    {
      res->header = 1;
      return(data.size());
    }
  if(strcmp(data[0].c_str(),"fixedStep") ==0) // line specifying type of data
    {
      *fixed = 1;
      *currSpan  = 1; // default for wig-tracks.
      for(int i=1;i<data.size();i++) // parse the parts of this line
	{
	  detailData.clear();
	  Tokenize(data[i],detailData,"="); 
	  if(strcmp(detailData[0].c_str(),"chrom") ==0)
	    currName->assign(detailData[1]);
	  else if(strcmp(detailData[0].c_str(),"start") ==0)
	    *currPos = atoi(detailData[1].c_str());
	  else if(strcmp(detailData[0].c_str(),"step") ==0)
	    *currStep = atoi(detailData[1].c_str());
	  else if(strcmp(detailData[0].c_str(),"span") ==0)
	    *currSpan = atoi(detailData[1].c_str());
	}
      res->header = 1;
      return(data.size());
    }

  if(strcmp(data[0].c_str(),"variableStep") ==0) 
    {
      *fixed = 0;
      *currSpan  = 1; // default for wig-tracks.
      for(int i=1;i<data.size();i++) // parse the parts of this line
	{
	  detailData.clear();
	  Tokenize(data[i],detailData,"="); 
	  if(strcmp(detailData[0].c_str(),"chrom") ==0)
	    currName->assign(detailData[1]);
	  else if(strcmp(detailData[0].c_str(),"step") ==0)
	    *currStep = atoi(detailData[1].c_str());
	  else if(strcmp(detailData[0].c_str(),"span") ==0)
	    *currSpan = atoi(detailData[1].c_str());
	}
      res->header = 1;
      return(data.size());
    }
  
  // data line.
  res->seq = *currName;
  if(*fixed)
    {
      res->pos = *currPos;
      // update position with the current step. 
      *currPos += *currStep; 
      res->val = (storageType)((double)scale * atof(data[0].c_str()));
    }
  else
    {
      // in variable step, just use the position given in the file
      res->pos = atoi(data[0].c_str());
      res->val = (storageType)((double)scale * atof(data[1].c_str()));
    }
  res->len = *currSpan;
  res->strand = 1; // no strand info in wig-files.
  res->mapped = 1;
  return(data.size());
}

/*
 * returns the number of lines with mapped reads 
 */
int initControlWIG(ifstream * inf,map<string,seqStats> * seqmap,int both,bool verbose,int scale)
{
  int nlines = 0;
  seqStats tmpStat;
  string line;
  inputLine tmp;
  string currName;
  int currPos,currStep,currSpan,fixed;
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
      if(parseWIGline(line,&tmp,&currName,&currStep,&currPos,&currSpan,&fixed,scale) < 1) // skip possible empty lines 
	continue;
      if(!tmp.header)
	{
	  if(tmp.mapped == 1) // query is mapped
	    {
	      //cout<<"data: "<<tmp.strand<<" "<<tmp.pos<<" "<<tmp.len<<" "<<tmp.seq<<" "<<tmp.mapped<<endl;
	      nlines++;
	      // try and insert into seqmap.
	      tmpStat.minPos = max(tmp.pos-both,0); // WIG supposed to be 1-based. but allow zero anyway. 
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
	}else{ //header
	//cout<<"Header found:"<<endl;
	//cout<<currName<<"\t"<<currPos<<"\t"<<currStep<<"\t"<<currSpan<<"\t"<<fixed<<endl;
      }
    }
  
  // reset the file-stream  
  inf->clear(); 
  inf->seekg(0,ios::beg);
  return(nlines);
}

