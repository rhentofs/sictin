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



#ifndef __WIG_UTILS_H__
#define __WIG_UTILS_H__

#include "utils.h"


/*
 *  
 *  Extract the necessary information from a single line of the WIG-file
 *
 */
int parseWIGline(string str,inputLine*,string* currName, 
		 int* currStep,int* currPos,int* currSpan,
		 int* fixed,int scale);


/*
 *
 * returns number of lines in infile + high/low coordinates of all seqs. 
 */
int initControlWIG(ifstream * inf,map<string,seqStats> * seqmap,int both,bool verbose,int scale);

#endif
