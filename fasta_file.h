//----------------------------------------------------------------
#ifndef FASTA_FILE_H
#define FASTA_FILE_H
//----------------------------------------------------------------
#include <map>
#include <string>
#include "hts.h"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class FastaEntry
{
private:
public:

    hts_pos_t len;
    char * seq;

    FastaEntry(void){}
    FastaEntry(hts_pos_t len,char * seq) : len(len),seq(seq){}
};
//----------------------------------------------------------------
class FastaFile
{
private:

    size_t totalSize; void * rawFile;

public:

    map<string,FastaEntry> entries;

    FastaFile(void);
    FastaFile(const string & filename);
    ~FastaFile(void);

    void Open(const string & filename);
};
//----------------------------------------------------------------
#endif
