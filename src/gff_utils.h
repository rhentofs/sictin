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



#ifndef __GFF_UTILS_H__
#define __GFF_UTILS_H__

#include "utils.h"


/*
 *  
 *  Extract the necessary information from a single line of the GFF-file
 *
 */
int parseGFFline(string str,inputLine*,string delim,int seqCol,int startCol,
		 int endCol,int strandCol,int valCol);


/*
 *
 * returns number of lines in infile + high/low coordinates of all seqs. 
 */
int initControlGFF(ifstream * inf,map<string,seqStats> * seqmap,string delim,int seqCol,int startCol,
		   int endCol,int strandCol,int valCol,int both,bool verbose);

#endif
