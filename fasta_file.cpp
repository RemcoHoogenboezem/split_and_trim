//----------------------------------------------------------------
// Name        : fasta_file.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description : Simple fasta file implementation
//----------------------------------------------------------------
#include <stdexcept>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include "faidx.h"
#include "fasta_file.h"
//----------------------------------------------------------------
FastaFile::FastaFile(void) : rawFile(MAP_FAILED){}
//----------------------------------------------------------------
FastaFile::FastaFile(const string & filename) : rawFile(MAP_FAILED) {Open(filename);}
//----------------------------------------------------------------
FastaFile::~FastaFile(void)
{
    if(this->rawFile!=MAP_FAILED) munmap(this->rawFile,this->totalSize);
}
//----------------------------------------------------------------
void FastaFile::Open(const string & filename)
{
    string prefix=filename.substr(0,filename.find_last_of('.'));

    //----------------------------------------------------------------
    //Try to create a lock file (Try it for 2 minutes)
    //----------------------------------------------------------------

    int lockFileID;
    size_t attempts=120;
    string lockFilename=prefix+".raw.lock";

    for(;;)
    {
        lockFileID=open(lockFilename.c_str(),O_EXCL|O_CREAT|O_WRONLY,0664);
        if(lockFileID!=-1) break;
        if(--attempts==0) throw runtime_error(string("Error: Could not create lock file (Maybe it already exists): ")+lockFilename);
        sleep(1);
    }

    //----------------------------------------------------------------
    //Open fasta index
    //----------------------------------------------------------------

    faidx_t * fai=fai_load(filename.c_str());

    if(fai==nullptr)
    {
        close(lockFileID); remove(lockFilename.c_str());
        throw runtime_error(string("Error: Could not open fasta index: ")+filename);
    }

    int nSeq=faidx_nseq(fai); size_t totalSize=0;
    for(int i=0;i<nSeq;i++) totalSize+=size_t(faidx_seq_len64(fai,faidx_iseq(fai,i)));

    //----------------------------------------------------------------
    //Create new raw file to only hold the sequence
    //----------------------------------------------------------------

    map<string,FastaEntry> entries;

    string rawFilename=prefix+".raw";
    int rawFileID=open(rawFilename.c_str(),O_RDONLY);

    if(rawFileID==-1)
    {
        rawFileID=open(rawFilename.c_str(),O_RDWR|O_CREAT|O_TRUNC,0664);

        if(rawFileID==-1)
        {
            close(rawFileID);
            fai_destroy(fai);
            close(lockFileID); remove(lockFilename.c_str());
            throw runtime_error("Error: Could not create raw file");
        }

        if(ftruncate(rawFileID,__off_t(totalSize))==-1)
        {
            close(rawFileID);
            fai_destroy(fai);
            close(lockFileID); remove(lockFilename.c_str());
            throw runtime_error("Error: Could not truncate raw file");
        }

        void * rawFile=mmap(nullptr,totalSize,PROT_READ|PROT_WRITE,MAP_SHARED,rawFileID,0);

        close(rawFileID);

        if(rawFile==MAP_FAILED)
        {
            fai_destroy(fai);
            close(lockFileID); remove(lockFilename.c_str());
            throw runtime_error("Error: Could not memory map raw file");
        }

        char * entrySeq=(char*)rawFile;

        for(int i=0;i<nSeq;i++)
        {
            const char * entryName=faidx_iseq(fai,i);
            hts_pos_t entryLen; char * tempSeq=fai_fetch64(fai,entryName,&entryLen);

            if(tempSeq==nullptr)
            {
                fai_destroy(fai);
                close(lockFileID); remove(lockFilename.c_str());
                throw runtime_error("Error: Could not fetch sequence");
            }

            for(char *src=tempSeq,*srcEnd=tempSeq+entryLen,*dest=entrySeq;src<srcEnd;src++,dest++){*dest=toupper(*src);} free(tempSeq);

            entries.emplace(entryName,FastaEntry(entryLen,entrySeq));
            entrySeq+=entryLen;
        }

        fai_destroy(fai);

        if(this->rawFile!=MAP_FAILED) munmap(this->rawFile,this->totalSize);
        this->totalSize=totalSize; this->rawFile=rawFile;
        this->entries=entries;

        close(lockFileID); remove(lockFilename.c_str());

        return;
    }

    //----------------------------------------------------------------
    //Open existing raw file
    //----------------------------------------------------------------

    void * rawFile=mmap(nullptr,size_t(totalSize),PROT_READ,MAP_SHARED,rawFileID,0);

    close(rawFileID);

    if(rawFile==MAP_FAILED)
    {
        fai_destroy(fai);
        close(lockFileID); remove(lockFilename.c_str());
        throw runtime_error("Error: Could not memory map raw file");
    }

    char * entrySeq=(char*)rawFile;

    for(int i=0;i<nSeq;i++)
    {
        const char * entryName=faidx_iseq(fai,i); hts_pos_t entryLen=faidx_seq_len64(fai,entryName);
        entries.emplace(entryName,FastaEntry(entryLen,entrySeq));
        entrySeq+=entryLen;
    }

    fai_destroy(fai);

    if(this->rawFile!=MAP_FAILED) munmap(this->rawFile,this->totalSize);
    this->totalSize=totalSize; this->rawFile=rawFile;
    this->entries=entries;

    close(lockFileID); remove(lockFilename.c_str());
}
//----------------------------------------------------------------
