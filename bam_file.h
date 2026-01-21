//----------------------------------------------------------------
#ifndef BAM_FILE_H
#define	BAM_FILE_H
//----------------------------------------------------------------
#include <string>
#include <stdexcept>
#include "sam.h"
//----------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------
class BamFile
{
private:

    bool canRead;
    bool canWrite;

    htsFile * handle;
    hts_idx_t * index;

    hts_itr_t * readIterator;

    void Close(void);

public:

    bam_hdr_t * header;
    bam1_t * read;

    BamFile(void);
    BamFile(const string & filename);
    BamFile(const string & filename,const bam_hdr_t * header,const string & program,const string & version,const string & command);
    ~BamFile(void);

    void Open(const string & filename);
    void Create(const string & filename,const bam_hdr_t * header,const string & program,const string & version,const string & command);
    void SetRegion(const string & region);
    inline int ReadRegion(void) {int ret=hts_itr_next(handle->fp.bgzf,readIterator,read,handle); if(ret<-1) throw runtime_error("Error: Could not read bam file"); return ret;}
    inline bool Write(bam1_t * read) {return sam_write1(handle,header,read)>=0;}
};
//----------------------------------------------------------------
#endif // BAM_FILE_H
