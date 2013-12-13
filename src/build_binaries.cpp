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
#include "sam_utils.h"
#include "gff_utils.h"
#include "bed_utils.h"
#include "wig_utils.h"

#include <dirent.h>

/*
 * updateBlock, updates the given coordinates in the given file with the given value.
 * returns how many (if any) values that was truncated.
 */



int updateBlock(fstream * fs,int from, int to,int minc,storageType score)
{
  int truncCnt = 0;
  // declare and allocate memory.
  storageType *tmpSignal;
  tmpSignal = new storageType[to-from +1];
  // read the existing data
  fs->seekg(sizeof(storageType)*(from-minc)+sizeof(int),ios::beg);
  fs->read((char*)tmpSignal,sizeof(storageType)*(to-from +1));
  if(!fs->good())
    {
      cerr<<endl<<from<<"\t"<<minc<<endl;
      cerr<<"updatateBlock:read failed"<<endl;
      delete[] tmpSignal;
      throw ios_base::failure("read");
    }
  // update, truncate at MAXNR
  for (int i=0;i<=(to-from);i++)
    if(tmpSignal[i]<(MAXNR-score))
      {
		tmpSignal[i]+=score;
      }
	else
      {
		tmpSignal[i] = MAXNR;
		truncCnt++;
      }
  // write back.

  fs->seekg(sizeof(storageType)*(from-minc)+sizeof(int),ios::beg);
  fs->write((char*)tmpSignal,sizeof(storageType)*(to-from +1));
  // make sure we're finished.
  fs->flush();
  if(!fs->good())
    {
      if(fs->eof()) // not a problem, just reset the filestream and exit normally.
		fs->clear();
      else
		{
		  cerr<<endl<<from<<"\t"<<minc<<endl;
		  cerr<<"updatateBlock: write failed"<<endl;
		  delete[] tmpSignal;
		  throw ios_base::failure("write");
		}
    }
  // release used memory.
  delete[] tmpSignal;
  return(truncCnt);
}


/*
 * create a empty (full of zeroes) binary file with the given number of positions.
 */

void createEmptyBinary (fstream * fs,int minc,int maxc,bool progress)
{

  int nrPos = maxc-minc;
  int wrPos = 0;
  storageType tmpShort = 0;

  storageType * tmpSignal;
  /*
   * Write the start coordinate offset and fill with zeroes.
   * Large chunks when possible and then one at a time.
   */
  fs->write((char*)&minc,sizeof(int));
  tmpSignal = new storageType[CHUNK_SIZE];
  for (int i=0;i<CHUNK_SIZE;i++)
    tmpSignal[i] = 0;

  while ((nrPos - wrPos >= CHUNK_SIZE))
    {
      fs->write((char*)tmpSignal,CHUNK_SIZE*sizeof(storageType));
      if(!fs->good())
		{
		  cerr<<"createEmptyBinary: write failed"<<endl;
		  delete[] tmpSignal;
		  throw ios_base::failure("write");
		}
      wrPos += CHUNK_SIZE;
      if(progress)
		cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<(int)(100*((double)wrPos/(double)nrPos))<<" % complete.\r";
    }
  delete[] tmpSignal;
  for (;wrPos<=nrPos;wrPos++)
    {
      if(progress && (wrPos % 1000) == 0)
		cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<(int)(100*((double)wrPos/(double)nrPos))<<" % complete.\r";
      fs->write((char*)&tmpShort,sizeof(storageType));
      if(!fs->good())
		{
		  cerr<<"createEmptyBinary: write failed"<<endl;
		  delete[] tmpSignal;
		  throw ios_base::failure("write");
		}
    }
  // just to make sure.
  fs->flush();
  if(progress)
    cerr<<" 100 % complete."<<endl;
}


/************************************************************
 *
 * main
 *
 ************************************************************/

int main(int argc, char* argv[])
{
  // parameters that you can set.
  int both        = 0;   // default No, if > 0 prolong the fragments this much.
  string delim    = "\t ";
  string infile   = "";
  string outfile  = "";
  string path     = "";
  int seqPos      = -1;
  int scPos       = -1;
  int fsPos       = -1;
  int fePos       = -1;
  int so          = 0;  // default No.
  int SAM         = 0;  //input file format
  int GFF         = 0;  //input file format
  int BED         = 0;  //input file format
  int WIG         = 0;  //input file format
  int scale       = 1;
  int valPos      = -1;

  string errorLine = "usage " +
    string(argv[0]) +
    " [parameters]\n" +
    "\t-i <infile>\n" +
    "\t-SAM <infile format, outfiles will be created using the \n" +
    "\t\t sequence-names in the SAM-file adding R_, F_ (and C_)\n" +
    "\t\t prefixes, eg 'R_chr1.bin'>\n" +
    "\t-GFF <infile format,uses additional parameter below.>\n" +
    "\t\t-d <delimiter, default tab/whitespace>\n" +
    "\t\t-si <sequence indentifier column in the infile, default 1>\n" +
    "\t\t-fs <fragment start column in the infile, default 2>\n" +
    "\t\t-fe <fragment end column in the infile, default 3>\n" +
    "\t\t-sc <fragment strand column in the infile, default 4>\n" +
    "\t\t-vc <value column in the infile, if > 0 use the (integer) value in this column instead of '1'>\n" +
    "\t-BED <infile format, outfiles will be 1-based> \n" +
    "\t\t-si <sequence indentifier column in the infile, default 1>\n" +
    "\t\t-fs <fragment start column in the infile,0-based default 2>\n" +
    "\t\t-fe <fragment end column in the infile, open, default 3>\n" +
    "\t\t-sc <fragment strand column in the infile, default 6>\n" +
    "\t-WIG <infile format, outfiles will be 1-based> \n" +
    "\t\t-scale <scale values in the file with this number, usefull when integer values \n"+
    "\t\t        are used as storage and wig-file contails real numbers.>\n" +
    "\t-bl <create total ovelap signal (C_), if > 0 with this length. def 0>\n" +
    "\t-so <toggle StartOnly, ie only first base of reads will be recorded. def Off>\n " +
    "\t-o <outfile stub, R_, F_ (and C_) prefixes will be added. eg 'R_<stub>chr1.bin'>\n" +
    "\t-p <outfile path, where files are to be placed.>\n"
    ;


  bool fail = false;
  string failmessage = "";

  for (int i=1;i<argc;i++)
    {
      if(strcmp(argv[i],"-i") == 0)
		infile.assign(argv[++i]);
      else if(strcmp(argv[i],"-o") == 0)
		outfile.assign(argv[++i]);
      else if(strcmp(argv[i],"-p") == 0)
		path.assign(argv[++i]);
      else if(strcmp(argv[i],"-d") == 0)
		delim.assign(argv[++i]);
      else if(strcmp(argv[i],"-si") == 0)
		seqPos = atoi(argv[++i]) -1;// column X 0-based is X-1.
      else if(strcmp(argv[i],"-sc") == 0)
        scPos = atoi(argv[++i]) -1; // column X 0-based is X-1.
      else if(strcmp(argv[i],"-fs") == 0)
        fsPos = atoi(argv[++i]) -1; // column X 0-based is X-1.
      else if(strcmp(argv[i],"-fe") == 0)
        fePos = atoi(argv[++i]) -1; // column X 0-based is X-1.
      else if(strcmp(argv[i],"-vc") == 0)
        valPos = atoi(argv[++i]) -1; // column X 0-based is X-1.
      else if(strcmp(argv[i],"-bl") == 0)
		both = atoi(argv[++i]);
      else if(strcmp(argv[i],"-so") == 0)
		so = 1;
      else if(strcmp(argv[i],"-SAM") == 0)
		SAM = 1;
      else if(strcmp(argv[i],"-BED") == 0)
		BED = 1;
      else if(strcmp(argv[i],"-GFF") == 0)
		GFF = 1;
      else if(strcmp(argv[i],"-WIG") == 0)
		WIG = 1;
      else if(strcmp(argv[i],"-scale") == 0)
		scale = atoi(argv[++i]);
      else
		{
		  failmessage.assign("Unknown argument: ");
		  failmessage.append(argv[i]);
		  failmessage.append("\n");
		  fail = true;
		}
    }

  if(infile == "")
    {
      failmessage.append("infile (-i) must be specified.\n");
      fail = true;
    }
  // check and set the defaults for columns.
  if(seqPos == -1) // sequence name
    if(GFF)
      seqPos = 0;
    else if (BED)
      seqPos = 0;
  if(fsPos == -1) // fragment start
    if(GFF)
      fsPos = 1;
    else if (BED)
      fsPos = 1;
  if(fePos == -1) // fragment end
    if(GFF)
      fePos = 2;
    else if (BED)
      fePos = 2;
  if(scPos == -1) // strand info
    if(GFF)
      scPos = 3;
    else if (BED)
      scPos = 5;

  int nrSet = 0;
  if(GFF)
    nrSet++;
  if(BED)
    nrSet++;
  if(SAM)
    nrSet++;
  if(WIG)
    nrSet++;

  if(nrSet > 1)
    {
      failmessage.append("only one of -GFF,-SAM,-WIG and -BED can be set at the same time.\n");
      fail = true;
    }

  if(nrSet == 0)
    {
      failmessage.append("infile format (-GFF,-SAM,-WIG or -BED) must be specified.\n");
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
		  failmessage.append("Please create outfile path. It does not exist.\n");
		  fail = true;
		}
    }


  ifstream inf;
  inf.open(infile.c_str());

  if(!inf)
    {
      failmessage.append("Could not open infile (does the file exist?)\n.");
      fail = true;
    }
  // are we ok so far?
  if (fail)
    {
      cerr << endl << failmessage.c_str() << endl << errorLine << endl;
      return(-1);
    }

  /*
   * Parse and build.
   */
  map <string,seqStats> seqMap;
  map<string,seqStats>::iterator it;

  // Read the reference sequences and the range of each.
  cout<<"Checking input."<<endl;
  int nlines;
  if(SAM)
    nlines = initControlSAM(&inf,&seqMap,both,true);
  if(GFF)
    nlines = initControlGFF(&inf,&seqMap,delim,seqPos,fsPos,fePos,scPos,valPos,both,true);
  if(BED)
    nlines = initControlBED(&inf,&seqMap,seqPos,fsPos,fePos,scPos,both,true);
  if(WIG)
    nlines = initControlWIG(&inf,&seqMap,both,true,scale);
  int nrFails = 0;
  int maxFails = 10;
  cout<<"Input consists of "<<nlines<<" mapped fragments."<<endl;

  /*
   * decide on file names & write empty binaries.
   */
  string preR = "R_";
  string preF = "F_";
  string preC = "C_";
  string line;
  inputLine tmpLine;
  int line_cnt = 1;
  int from, to;
  map <string,fsHolder> seqFiles;
  map <string,fsHolder>::iterator fsIt;
  map <string,seqStats> seqStats;
  fsHolder *tmpHolder;
  // parameters needed by the WIG-format
  string currName;
  int currPos,currStep,currSpan,fixed;

  // print some stats.
  cout <<"Input Statistics::"<<endl;
  cout <<setw(10)<<"Name\t"<<setw(10)<<"minCrd\t"<<setw(10)<<"maxCrd\t"<<setw(10)<<"F_counts\t"<<setw(10)<<"R_counts\t"<<endl;
  for ( it=seqMap.begin() ; it != seqMap.end(); it++ )
    {
      cout <<setw(10)<< (*it).first << "\t" <<setw(10)<< (*it).second.minPos << "\t" << setw(10)<<(*it).second.maxPos<<"\t";
      cout <<setw(10)<< (*it).second.countF << "\t" <<setw(10)<< (*it).second.countR << endl;
    }

  int noBinsTotal = (both > 0 ? 3 : 2)*seqMap.size();
  double noBinsWritten = 0;
  cout << "Writing empty binaries:" <<endl;

  for ( it=seqMap.begin() ; it != seqMap.end(); it++ )
    {
      cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<100*noBinsWritten/noBinsTotal<<" % complete.\r";
      seqFiles.insert(pair<string,fsHolder>((*it).first,* new fsHolder(path + preR + outfile + (*it).first + ".bin",
																	   path + preF + outfile + (*it).first + ".bin",
																	   (both > 0 ? path + preC + outfile + (*it).first + ".bin" : ""),
																	   (*it).second.minPos,false)));
      try{
		createEmptyBinary(seqFiles[(*it).first].fsF,(*it).second.minPos,(*it).second.maxPos,false);
      }catch(ios_base::failure &f)
		{
		  cerr<<"Failure: '"<<f.what()<<"' aborting."<<endl;
		  goto exit_point;
		}
      cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<100*(++noBinsWritten)/noBinsTotal<<" % complete.\r";
      try{
		createEmptyBinary(seqFiles[(*it).first].fsR,(*it).second.minPos,(*it).second.maxPos,false);
      }catch(ios_base::failure &f)
		{
		  cerr<<"Failure: '"<<f.what()<<"' aborting."<<endl;
		  goto exit_point;
		}
      cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<100*(++noBinsWritten)/noBinsTotal<<" % complete.\r";
      if(both > 0)
		{
		  try{
			createEmptyBinary(seqFiles[(*it).first].fsC,(*it).second.minPos,(*it).second.maxPos,false);
		  }catch(ios_base::failure &f)
			{
			  cerr<<"Failure: '"<<f.what()<<"' aborting."<<endl;
			  goto exit_point;
			}
		  cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<100*(++noBinsWritten)/noBinsTotal<<" % complete.\r";
		}
    }

  cerr<<setw(10)<<100<<" % complete."<<endl;
  cout<<"Processing input:"<<endl;

  /*
   * process input-file and write out overlaps.
   */

  cerr<<"0 % complete.\r";
  while(!(inf.eof()))
    {
      if((line_cnt % 1000) == 0 )
		cerr<<setw(4)<<setfill(' ')<<setprecision(2)<<(int)(100*(double)line_cnt/(double)nlines)<<" % complete.\r";
      getline(inf,line);
      if(SAM)
		{
		  if(parseSAMline(line,&tmpLine) < 1) // skip possible empty lines
			{
			  continue;
			}
		  if(tmpLine.mapped == 0) // non-mapped read? skip.
			{
			  continue;
			}
		  if(tmpLine.header) // header line? skip.
			{
			  continue;
			}
		}
      if(GFF)
		{
		  if(parseGFFline(line,&tmpLine,delim,seqPos,fsPos,fePos,scPos,valPos) < 1) // skip possible empty lines and non-mapped.
			continue;
		}
      if(BED)
		{
		  if(parseBEDline(line,&tmpLine,seqPos,fsPos,fePos,scPos) < 1) // skip possible empty lines and non-mapped.
			continue;
		}
      if(WIG)
		{
		  if(parseWIGline(line,&tmpLine,&currName,&currStep,&currPos,&currSpan,&fixed,scale) < 1) // skip possible empty lines
			continue;
		}
      line_cnt++;
      // use tmpHolder so that we just have to do one lookup-call in seqFiles
      tmpHolder = &(seqFiles.find(tmpLine.seq))->second;
      if(tmpLine.strand == -1) // anti-sense strand
		{
		  from =  tmpLine.pos;
		  to   =  tmpLine.pos + tmpLine.len -1;
		  if(so == 1) // just using start crds
			from = to;

		  try
			{
			  seqMap[tmpLine.seq].truncR += updateBlock(tmpHolder->fsR,
														from,
														to,
														tmpHolder->minc,
														tmpLine.val);
			}
		  catch (ios_base::failure &f)
			{
			  cerr<<"Failure: '"<<f.what()<<"'\t";
			  if(++nrFails < maxFails)
				{
				  //seqFiles[tmpLine.seq].fsR->clear();
				  tmpHolder->fsR->clear();
				  cerr<<"Ignoring line "<<line_cnt<<" and continuing."<<endl;
				}
			  else
				{
				  goto exit_point;
				}
			}
		  if(both>0)
			{
			  from = max(to - both+1,0);
			  if(so == 1) // just using start crds
				from = to;
			  try
				{
				  seqMap[tmpLine.seq].truncC += updateBlock(tmpHolder->fsC,
															from,
															to,
															tmpHolder->minc,
															tmpLine.val);
				}
			  catch (ios_base::failure &f)
				{
				  cerr<<"Failure: '"<<f.what()<<"'\t";
				  if(++nrFails < maxFails)
					{
					  tmpHolder->fsC->clear();
					  cerr<<"Ignoring line "<<line_cnt<<" and continuing."<<endl;
					}
				  else
					{
					  goto exit_point;
					}
				}
			}
		} // anti sense strand

      if(tmpLine.strand == 1) // sense strand
		{
		  from =  tmpLine.pos;
		  to   =  tmpLine.pos + tmpLine.len -1;
		  if(so == 1) // just using start crds
			to = from;

		  try
			{
			  seqMap[tmpLine.seq].truncF += updateBlock(tmpHolder->fsF,
														from,
														to,
														tmpHolder->minc,
														tmpLine.val);
			}
		  catch (ios_base::failure &f)
			{
			  cerr<<"Failure: '"<<f.what()<<"'\t";
			  if(++nrFails < maxFails)
				{
				  //seqFiles[tmpLine.seq].fsF->clear();
				  tmpHolder->fsF->clear();
				  cerr<<"Ignoring line "<<line_cnt<<" and continuing."<<endl;
				}
			  else
				{
				  goto exit_point;
				}
			}
		  if(both>0)
			{
			  to = tmpLine.pos + both -1;
			  if(so == 1) // just using start crds
				to = from;
			  try
				{
				  seqMap[tmpLine.seq].truncC += updateBlock(tmpHolder->fsC,
															from,
															to,
															tmpHolder->minc,
															tmpLine.val);
				}
			  catch (ios_base::failure &f)
				{
				  cerr<<"Failure: '"<<f.what()<<"'\t";
				  if(++nrFails < maxFails)
					{
					  tmpHolder->fsC->clear();
					  cerr<<"Ignoring line "<<line_cnt<<" and continuing."<<endl;
					}
				  else
					{
					  goto exit_point;
					}
				}
			}
		} // sense strand


    }
  cerr<<setw(10)<<100<<" % complete."<<endl;

 exit_point:
  // echo some stats.
  cout <<"Signals were built using:"<<endl;
  cout <<setw(10)<<"Name\t"<<setw(10)<<"minCrd\t"<<setw(10)<<"maxCrd\t"<<setw(10)<<"F_counts\t"<<setw(10)<<"R_counts\t";
  cout <<setw(10)<<"Trunc_F\t"<<setw(10)<<"Trunc_R\t"<<setw(10)<<(both > 0 ? "Trunc_C" : "")<<endl;
  for ( it=seqMap.begin() ; it != seqMap.end(); it++ )
    {
      cout <<setw(10)<< (*it).first << "\t" <<setw(10)<< (*it).second.minPos << "\t" << setw(10)<<(*it).second.maxPos<<"\t";
      cout <<setw(10)<< (*it).second.countF << "\t" <<setw(10)<< (*it).second.countR << "\t";
      cout <<setw(10)<< (*it).second.truncF << "\t" << setw(10)<<(*it).second.truncR;
      if(both > 0)
		cout <<"\t"<< setw(10)<<(*it).second.truncC;
      cout <<endl;
    }

  string statFname = path + "buildStats.txt";
  bool writeStats = true;
  ofstream ofc;
  ofc.open(statFname.c_str(),ios::trunc);
  if (ofc.fail())
    {
      failmessage.clear();
      failmessage.append("ERROR: Output file \"");
      failmessage.append(statFname.c_str());
      failmessage.append("\" could not be created, skipping.\n");
      writeStats = false;
    }

  if(writeStats)
    {
      ofc <<"Name\t"<<"minCrd\t"<<"maxCrd\t"<<"F_counts\t"<<"R_counts\t";
      ofc <<"Trunc_F\t"<<"Trunc_R"<<(both > 0 ? "\tTrunc_C" : "")<<endl;
      for ( it=seqMap.begin() ; it != seqMap.end(); it++ )
		{
		  ofc << (*it).first << "\t" << (*it).second.minPos << "\t" << (*it).second.maxPos<<"\t";
		  ofc << (*it).second.countF << "\t" << (*it).second.countR << "\t";
		  ofc << (*it).second.truncF << "\t" << (*it).second.truncR;
		  if(both > 0)
			ofc <<"\t"<< (*it).second.truncC;
		  ofc <<endl;
		}
    }else{
    cerr<<failmessage.c_str()<<endl;
  }


  // close files.
  inf.close();
  ofc.close();

  return(0);
}




