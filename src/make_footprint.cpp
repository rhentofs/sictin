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
#include <dirent.h>

//using namespace std;

/************************************************************
 *
 * main 
 *
 ************************************************************/

int main(int argc, char* argv[]) 
{
  // parameters that you can set. 
  string path     = "";
  string stub     = "";
  string outfile  = "";
  string qfile    = "";
  string efile    = "";
  int    offset   = 5000;
  int    progress = 0;
  int    co       = 0;
  int    move     = 0;
  int    both     = 0;
  int    all      = 0;

  string errorLine =  "usage " + 
    string(argv[0]) + 
    " [parameters] \n" +
    "\t-ip <infile path> \n" +
    "\t-is <infile stub, eg '<stub>' part of 'R_<stub>chr1.bin'\n" +
    "\t-o <outfile stub, 'F_', 'R_' (and 'C_') prefixes will be added.>\n" +
    "\t-c <toggle combined ('C_'). Default off. > \n" + 
    "\t-q <tab separated query file, contains seq,position & strand (-/+) e.g. 'chr1 12300 +'>\n" + 
    "\t-e <errorfile, defaults to 'err.txt', failed queries end up here>\n" + 
    "\t-d <offset, how long distance +/- the positions to consider default:5000>\n" + 
    "\t-p <toggle progressmeter. default off>\n" + 
    "\t-s <score cutoff, at least on bp must have atleast this score to be included.  default 0>\n" +
    "\t-m <move offset, move all the queries this many bp, default 0>\n" +
    "\t-a <toggle write out all individual data points, default off.  Useful for e.g. heatmaps, adds 'all_' to the outfile stub. >\n" +
    " Writes 4 lines in each outfile \n" + 
    "\t1) >+ <number of successfully executed queries>\n" + 
    "\t2) space separated positive strand signal values \n" + 
    "\t3) >- <number of successfully executed queries>\n" + 
    "\t4) space separated negative strand signal values, reversed to direction of transcription. \n";
    
  bool fail = false;
  string failmessage = "";
  
  for (int i=1;i<argc;i++)
    {
      if(strcmp(argv[i],"-ip") == 0)
	path.assign(argv[++i]);
      else if (strcmp(argv[i],"-is") == 0)
	stub.assign(argv[++i]);
      else if (strcmp(argv[i],"-e") == 0)
	efile.assign(argv[++i]);
      else if(strcmp(argv[i],"-o") == 0)
	outfile.assign(argv[++i]);
      else if(strcmp(argv[i],"-c") == 0)
	both = 1;
      else if(strcmp(argv[i],"-q") == 0)
	qfile.assign(argv[++i]); 
      else if(strcmp(argv[i],"-d") == 0)
	offset = atoi(argv[++i]);
      else if(strcmp(argv[i],"-p") == 0)
	progress = 1;
      else if(strcmp(argv[i],"-s") == 0)
	co = atoi(argv[++i]);
      else if(strcmp(argv[i],"-m") == 0)
	move = atoi(argv[++i]);
      else if(strcmp(argv[i],"-a") == 0)
	all = 1;
      else
	{
	  failmessage.assign("Unknown argument: ");
	  failmessage.append(argv[i]);
	  failmessage.append("\n");
	  fail = true;
	}
    }
  
  if(strcmp(outfile.c_str(),"") == 0)
    {
       failmessage.append("outfile (-o) must be specified.\n");
       fail = true;
    }
  
  if(strcmp(qfile.c_str(),"") == 0)
    {
       failmessage.append("query file (-q) must be specified.\n");
       fail = true;
    }
  // make sure the path exists.
  if(strcmp(path.c_str(),"") != 0)
    {
      DIR *d;
      d = opendir(path.c_str());
      if(d)
	{
	  closedir(d);
	}
      else
	{
	  failmessage.append("The given path does not exist.\n");
	  fail = true;
	}
    }
  // default if not specified. 
  if(strcmp(efile.c_str(),"") == 0)
    efile.assign("err.txt");
  
  ifstream infq;
  infq.open(qfile.c_str());
  
  if(!infq)
    {
      failmessage.append("Could not open query-file(does the file exist?)\n");
      fail = true;
    }
  
  // are we ok so far? 
  if (fail)
    {
      cerr << endl << failmessage.c_str() << endl << errorLine << endl;
      return(-1);
    }
  
  map <string,seqStats> seqMap; 
  map<string,seqStats>::iterator it;
  
  // Read the queries and find out which sequences are to be queried.
  cout<<"Checking Queries."<<endl;
  int nq = queryControl(&infq,&seqMap,true);
  cout <<"Query Statistics:"<<endl;
  cout <<setw(10)<<"Name\t"<<setw(10)<<"minCrd\t"<<setw(10)<<"maxCrd\t"<<setw(10)<<"F_counts\t"<<setw(10)<<"R_counts\t"<<endl;
  for ( it=seqMap.begin() ; it != seqMap.end(); it++ )
    {
      cout <<setw(10)<< (*it).first << "\t" <<setw(10)<< (*it).second.minPos << "\t" << setw(10)<<(*it).second.maxPos<<"\t";
      cout <<setw(10)<< (*it).second.countF << "\t" <<setw(10)<< (*it).second.countR << endl;
    }
  string preR = "R_";
  string preF = "F_";
  string preC = "C_";
  
  // Create the file-pointers.
  map <string,fsHolder> seqFiles; 
  bool allFilesOk = true; 
  for ( it=seqMap.begin() ; it != seqMap.end(); it++ )
    {
      // there's a nasty side effect here: the files-pointers are opened in the copy-constructor. 
      seqFiles.insert(pair<string,fsHolder>((*it).first,* new fsHolder(path + preR + stub + (*it).first + ".bin",
								       path + preF + stub + (*it).first + ".bin",
								       (both > 0 ? path + preC + stub + (*it).first + ".bin" : ""),
								       0,true)));
      // check if the streams are ok.
      if(!seqFiles[(*it).first].fsR->good())
	{
	  cerr<<seqFiles[(*it).first].nameR<<" could not be opened."<<endl; 
	  allFilesOk = false;
	}
      if(!seqFiles[(*it).first].fsF->good())
	{
	  cerr<<seqFiles[(*it).first].nameF<<" could not be opened."<<endl; 
	  allFilesOk = false;
	}
      if((both > 0) && !seqFiles[(*it).first].fsC->good())
	{
	  cerr<<seqFiles[(*it).first].nameC<<" could not be opened."<<endl; 
	  allFilesOk = false;
	}
    }
  if(!allFilesOk)
    {
      cerr<<"All necessary input files could no be opened. Aborting."<<endl;
      goto exit_point;
    }
  
  // call the footprinter. 
  if(getMultipleSignalsFP(&infq,seqFiles,outfile,efile,
			  offset,progress,co,move,both,nq,all)<0)
    {
      cerr<<"Footprinting failed."<<endl;
      goto exit_point;
    }

 exit_point:
  infq.close();

  return(0);
}
  
  

  
