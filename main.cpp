//----------------------------------------------------------------
// Name        : main.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include "split_N_cigar.h"
//----------------------------------------------------------------
int main(int argc,char * argv[])
{
    //char command[]="/tools/split_N_cigar/split_N_cigar -i /home/remco/Desktop/test_split_N_cigar/2194.bam -r /home/remco/Desktop/test_split_N_cigar/GRCh38.primary_assembly.genome.fa -s 60 -R 2 -p 16 -c -t"; char * pCommand=command;
    //int argc=13; char * argv[16]; for(int i=0;i<argc;i++) argv[i]=strsep(&pCommand," ");

    return SplitNCigar::Run(argc,argv);
}
//----------------------------------------------------------------

