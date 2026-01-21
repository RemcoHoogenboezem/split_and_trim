//----------------------------------------------------------------
#ifndef SPLIT_N_CIGAR_H
#define SPLIT_N_CIGAR_H
//----------------------------------------------------------------
#include <map>
#include "fasta_file.h"
#include "bam_file.h"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class SplitNCigar
{
private:

    bool useClipping;
    bool useTrimming;

    uint8_t maxAlignmentScore;

    size_t nSupportReads;
    size_t minMismatches;
    size_t maxOverhang;

    const map<string,FastaEntry> & fastaEntries;
    const char * referenceSeq;

    BamFile inputBamFile,outputBamFile;

    map<hts_pos_t,size_t> intronBegin;
    map<hts_pos_t,size_t> intronEnd;

    SplitNCigar(bool useClipping,bool useTrimming,uint8_t maxAlignmentScore,size_t nSupportReads,size_t minMismatches,size_t maxOverhang,const map<string,FastaEntry> & fastaEntries,const string & inputBamFilename);

    void TrimBegin(bam1_t ** pRead,hts_pos_t endPos);
    void TrimEnd(bam1_t ** pRead,hts_pos_t endPos);
    void Run(size_t i);

public:

    static int Run(int argc,char *argv[]);
};
//----------------------------------------------------------------
#endif // SPLIT_N_CIGAR_H
