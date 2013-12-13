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


/************************************************************
 *
 * main 
 *
 ************************************************************/

int main(int argc, char* argv[]) 
{
  // parameters that the user can set. 
  string infile  = "";
  string outfile = "";
  string qfile   = "";
  string efile   = "";
  int from       = -1;
  int to         = -1;
  int avg        = 0;
  string errSymb = "";

  string errorLine =  "usage " + 
    string(argv[0]) + 
    " [parameters] \n" +
    "\t-i <infile> \n" + 
    "\t-o <outfile> \n" +
    "\t-q <query file, contains start & stop>\n" + 
    "\t-from <if single query, from this bp>\n" + 
    "\t-to <if signle query, to (including) this bp>\n"+
    "\t-e <errorfile, defaults to 'err.txt', failed queries end up here>\n"+
    "\t-errSymb <errorsymbol, disables headers in the 'outfile' and 'errorfile' completely. Writes the symbol instead in 'outfile'>\n"+
    "\t-avg <toggle average(mean) signal in each query only (def. no)>\n";

  for (int i=1;i<argc;i++)
    {
      if(strcmp(argv[i],"-i") == 0)
	infile.assign(argv[++i]);
      else if (strcmp(argv[i],"-e") == 0)
	efile.assign(argv[++i]);
      else if (strcmp(argv[i],"-errSymb") == 0)
	errSymb.assign(argv[++i]);
      else if(strcmp(argv[i],"-o") == 0)
	outfile.assign(argv[++i]);
      else if(strcmp(argv[i],"-q") == 0)
	qfile.assign(argv[++i]);
      else if(strcmp(argv[i],"-from") == 0)
	from = atoi(argv[++i]);
      else if(strcmp(argv[i],"-to") == 0)
	to = atoi(argv[++i]);
      else if(strcmp(argv[i],"-avg") == 0)
	avg = 1;
      else
	{
	  cerr<<"Unknown argument: "<<argv[i]<<endl<<errorLine<<endl;
	  return(-1);
	}
    }
  
  if(strcmp(infile.c_str(),"") == 0)
    {
      cerr<<"Infile (-i) must be specified"<<endl<<errorLine<<endl;
      return(-1);
    }
  
  if(strcmp(outfile.c_str(),"") == 0)
    {
       cerr<<"Outfile (-o) must be specified"<<errorLine<<endl;
       return(-1);
    }
  
  ofstream outf;
  outf.open(outfile.c_str());
  if(!outf.is_open())
    {
      cerr<<"Couldn't open outfile'"<<outfile<<"'"<<endl;
      return(-1);
    }

  if(strcmp(qfile.c_str(),"")==0 && (from == -1 || to == -1))
    {
      cerr<<"queryfile (-q) must be specified OR both 'from' and 'to'"<<endl<<errorLine<<endl;
      outf.close();
      return(-1);
    }


  
  // if not specified used the default name.
  if(strcmp(efile.c_str(),"") == 0)
    efile.assign("err.txt");
  
  // single query
  if (from >= 0 && to >= 0)
    {
      if(from > to)
	{
	  cerr<<"Error: 'from' must be less than 'to'"<<endl;
	  return(-1);
	}
      storageType *tmp;
      tmp = new storageType [to-from+1];
      if(getSignalC(infile.c_str(),from,to,tmp) <0) // something failed
	{
	  cerr<<"Error: getSignalC failed."<<endl;
	  delete[] tmp;
	  return(-1);
	}
      else // data was successfully read.
	{
	  if(avg == 0) // should the signal written as is?
	    {
	      if (strcmp(errSymb.c_str(),"") == 0)
		outf<<"> "<<from<<"-"<<to<<endl;
	      outf<<tmp[0];
	      for (int i = 1;i<(to-from+1);i++)
		{
		  outf<<" "<<tmp[i];
		}
	      outf<<endl;
	    }
	  else // should the signal be averaged?
	    {
	      double tmpsum = 0;
	      for (int i = 0;i<(to-from+1);i++)
		tmpsum += tmp[i];
	      tmpsum = (double)tmpsum/(to-from+1);
	      if (strcmp(errSymb.c_str(),"") == 0)
		outf<<"> "<<from<<"-"<<to<<endl;
	      outf<<tmpsum<<endl;
	    }
	}
      outf.close();
      delete[] tmp;
    }
  else // multiple queries are handled by getSignalsTxt
    {
      if(getSignalsTxt(infile.c_str(),qfile.c_str(), outfile.c_str(),efile.c_str(),avg,errSymb.c_str())<0)
	{
	  cerr<<"Error: getSignalsTxt failed."<<endl;
	  return(-1);
	}
    }
  return(0);
}
  
  

  
