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



#ifndef __SAM_UTILS_H__
#define __SAM_UTILS_H__

#include "utils.h"

/*
 * Bits defined in the SAM-FLAG. 
 *
 */

#define SAM_PAIRED_READ      0x0001
#define SAM_MAPPED_PAIR      0x0002
#define SAM_UNMAPPED_SEQ     0x0004
#define SAM_UNMAPPED_MATE    0x0008
#define SAM_REV_STRAND       0x0010
#define SAM_STRAND_MATE      0x0020
#define SAM_PAIR_FIRST       0x0040
#define SAM_PAIR_SECOND      0x0080
#define SAM_NON_PRIMARY      0x0100
#define SAM_READ_QUAL_FAIL   0x0200
#define SAM_READ_PCR_OPT_DUP 0x0400


/*
 *  
 *  Extract the necessary information from a single line of the SAM-file
 *
 */
int parseSAMline(string str,inputLine*);


/*
 *
 * returns number of lines in infile + high/low coordinates of all seqs. 
 */
int initControlSAM(ifstream * inf,map<string,seqStats> * seqmap,int both,bool verbose);

#endif
