//----------------------------------------------------------------
// Name        : bam_file.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include "bam_file.h"
//----------------------------------------------------------------
void BamFile::Close(void)
{
    if(canRead)
    {
        canRead=false;

        if(readIterator!=nullptr){sam_itr_destroy(readIterator);readIterator=nullptr;}
        bam_destroy1(read);
        hts_idx_destroy(index);index=nullptr;
        sam_hdr_destroy(header);header=nullptr;
        sam_close(handle);handle=nullptr;
        return;
    }

    if(canWrite)
    {
        canWrite=false;

        sam_hdr_destroy(header);header=nullptr;
        sam_close(handle);handle=nullptr;
        return;
    }
}
//----------------------------------------------------------------
BamFile::BamFile(void) : canRead(false),canWrite(false),readIterator(nullptr) {}
BamFile::BamFile(const string & filename) : canRead(false),canWrite(false),readIterator(nullptr) {Open(filename);}
BamFile::BamFile(const string & filename,const bam_hdr_t * header,const string & program,const string & version,const string & command) : canRead(false),canWrite(false),readIterator(nullptr) {Create(filename,header,program,version,command);}
BamFile::~BamFile(void) {Close();}
//----------------------------------------------------------------
void BamFile::Open(const string & filename)
{
    //----------------------------------------------------------------
    //Open bam file
    //----------------------------------------------------------------

    htsFile * handle=sam_open(filename.c_str(),"r");
    if(handle==nullptr) throw runtime_error(string("Error: Could not open bam file: ")+filename);

    //----------------------------------------------------------------
    //Read header
    //----------------------------------------------------------------

    bam_hdr_t * header=sam_hdr_read(handle);

    if(header==nullptr)
    {
        sam_close(handle);
        throw runtime_error(string("Error: Could not read bam header: ")+filename);
    }

    //----------------------------------------------------------------
    //Load index
    //----------------------------------------------------------------

    hts_idx_t * index=sam_index_load(handle,filename.c_str());

    if(index==nullptr)
    {
        sam_hdr_destroy(header);
        sam_close(handle);
        throw runtime_error(string("Error: Could not read bam index: ")+filename);
    }

    //----------------------------------------------------------------
    //Update object
    //----------------------------------------------------------------

    Close();

    this->handle=handle;
    this->header=header;
    this->index=index;
    read=bam_init1();

    canRead=true;
}
//----------------------------------------------------------------
void BamFile::Create(const string & filename,const bam_hdr_t * header,const string & program,const string & version,const string & command)
{
    //----------------------------------------------------------------
    //Create bam file
    //----------------------------------------------------------------

    htsFile * handle=hts_open(filename.c_str(),"wb");
    if(handle==nullptr) throw runtime_error(string("Error: Could not create bam file: ")+filename);

    //----------------------------------------------------------------
    //Create header
    //----------------------------------------------------------------

    bam_hdr_t * dupHeader=sam_hdr_dup(header);

    if(dupHeader==nullptr)
    {
        sam_close(handle);
        throw runtime_error(string("Error: Could not duplicate bam header: ")+filename);
    }

    //----------------------------------------------------------------
    //Add program line to bam header
    //----------------------------------------------------------------

    if(sam_hdr_add_line(dupHeader,"PG","ID",program.c_str(),"PN",program.c_str(),"VN",version.c_str(),"CL",command.c_str(),nullptr)==-1)
    {
        sam_hdr_destroy(dupHeader);
        sam_close(handle);
        throw runtime_error(string("Error: Could not add line to bam header: ")+filename);
    }

    //----------------------------------------------------------------
    //Write header to bam file
    //----------------------------------------------------------------

    if(sam_hdr_write(handle,dupHeader)<0)
    {
        sam_hdr_destroy(dupHeader);
        sam_close(handle);
        throw runtime_error(string("Error: Could not write bam header: ")+filename);
    }

    //----------------------------------------------------------------
    //Update object
    //----------------------------------------------------------------

    Close();

    this->handle=handle;
    this->header=dupHeader;

    canWrite=true;
}
//----------------------------------------------------------------
void BamFile::SetRegion(const string & region)
{
    //----------------------------------------------------------------
    //Check if the bam file is initialized
    //----------------------------------------------------------------

    if(canRead==false) throw runtime_error("Error: No open bam file for reading");

    //----------------------------------------------------------------
    //Create iterator
    //----------------------------------------------------------------

    hts_itr_t * readIterator=sam_itr_querys(index,header,region.c_str());
    if(readIterator==nullptr) throw runtime_error(string("Error: Could not init read iterator: ")+handle->fn);

    if(this->readIterator!=nullptr) sam_itr_destroy(this->readIterator);
    this->readIterator=readIterator;
}
//----------------------------------------------------------------
