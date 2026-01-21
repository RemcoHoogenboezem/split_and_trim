//----------------------------------------------------------------
// Name        : split_N_cigar.cpp
// Author      : Remco Hoogenboezem
// Version     :
// Copyright   :
// Description :
//----------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <getopt.h>
#include "fasta_file.h"
#include "split_N_cigar.h"
//----------------------------------------------------------------
SplitNCigar::SplitNCigar(bool useClipping,bool useTrimming,uint8_t maxAlignmentScore,size_t nSupportReads,size_t minMismatches,size_t maxOverhang,const map<string,FastaEntry> & fastaEntries,const string & inputBamFilename) : useClipping(useClipping),useTrimming(useTrimming),maxAlignmentScore(maxAlignmentScore),nSupportReads(nSupportReads),minMismatches(minMismatches),maxOverhang(maxOverhang),fastaEntries(fastaEntries)
{
    inputBamFile.Open(inputBamFilename);
}
//----------------------------------------------------------------
void SplitNCigar::TrimBegin(bam1_t ** pRead,hts_pos_t endPos)
{
    bam1_t * read=*pRead;
    size_t cigarLen=read->core.n_cigar;
    uint32_t * cigar=bam_get_cigar(read);
    hts_pos_t pos=read->core.pos; if(bam_cigar_op(cigar[0])==BAM_CSOFT_CLIP) pos-=bam_cigar_oplen(cigar[0]);

    hts_pos_t junctionPos; for(auto junction=intronEnd.lower_bound(pos);;junction++)
    {
        if(junction==intronEnd.end() || junction->first>=endPos) return;
        if(junction->second>=nSupportReads){junctionPos=junction->first; break; }
    }

    size_t queryIndex=0;
    uint8_t * seq=(uint8_t*)bam_get_seq(read);
    char * qual=(char*)bam_get_qual(read);

    size_t newCigarLen=0;
    uint32_t newCigar[256];

    size_t newQueryLen=0;
    char newSeq[2048]; char newQual[2048];

    size_t nBasesCovered=0;
    size_t nBasesToClip=0;
    size_t nMismatches=0;

    size_t i=0; for(;i<cigarLen;i++)
    {
        uint32_t op=bam_cigar_op(cigar[i]);
        uint32_t opLen=bam_cigar_oplen(cigar[i]);

        switch(op)
        {
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:

            for(;opLen>0 && pos<=junctionPos;opLen--,pos++,queryIndex++)
            {
                nBasesCovered++;nBasesToClip++; nMismatches+=seq_nt16_table[size_t(referenceSeq[pos])]!=bam_seqi(seq,queryIndex);
            }

            if(pos<=junctionPos) continue;
            if(nBasesCovered>maxOverhang || nMismatches<minMismatches) return;

            if(opLen!=0)
            {
                newCigar[newCigarLen++]=bam_cigar_gen(nBasesToClip,BAM_CHARD_CLIP);
                newCigar[newCigarLen++]=bam_cigar_gen(opLen,op);

                for(;opLen>0;opLen--,queryIndex++,newQueryLen++)
                {
                    newSeq[newQueryLen]=seq_nt16_str[bam_seqi(seq,queryIndex)]; newQual[newQueryLen]=qual[queryIndex];
                }
                goto COPY_REMAINDER;
            }

            op=bam_cigar_op(cigar[i+1]); if(op==BAM_CINS || op==BAM_CDEL)
            {
                nBasesToClip+=bam_cigar_oplen(cigar[++i]); newCigar[newCigarLen++]=bam_cigar_gen(nBasesToClip,BAM_CHARD_CLIP);
            }

            goto COPY_REMAINDER;

        case BAM_CINS:

            nBasesToClip+=opLen;
            queryIndex+=opLen;
            continue;

        case BAM_CDEL:

            if(pos<=junctionPos && junctionPos<pos+opLen)
            {
                if(nBasesCovered>maxOverhang || nMismatches<minMismatches) return;

                pos+=opLen;
                newCigar[newCigarLen++]=bam_cigar_gen(nBasesToClip,BAM_CHARD_CLIP);

                goto COPY_REMAINDER;
            }

            pos+=opLen;
            continue;

        case BAM_CSOFT_CLIP:

            for(;opLen>0;opLen--,pos++,queryIndex++)
            {
                nBasesCovered++;nBasesToClip++; nMismatches+=seq_nt16_table[size_t(referenceSeq[pos])]!=bam_seqi(seq,queryIndex);
            }

            if(pos>junctionPos)
            {
                if(nBasesCovered>maxOverhang || nMismatches<minMismatches) return;

                newCigar[newCigarLen++]=bam_cigar_gen(nBasesToClip,BAM_CHARD_CLIP);

                goto COPY_REMAINDER;
            }

            continue;
        }
    }
COPY_REMAINDER:

    for(i++;i<cigarLen;i++)
    {
        newCigar[newCigarLen++]=cigar[i];

        switch(bam_cigar_op(cigar[i]))
        {
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
        case BAM_CINS:
        case BAM_CSOFT_CLIP:

            uint32_t opLen=uint32_t(bam_cigar_oplen(cigar[i]));

            for(;opLen>0;opLen--,queryIndex++,newQueryLen++)
            {
                newSeq[newQueryLen]=seq_nt16_str[bam_seqi(seq,queryIndex)]; newQual[newQueryLen]=qual[queryIndex];
            }

            continue;
        }
    }

    bam1_t * newRead=bam_init1();
    size_t l_aux=bam_get_l_aux(read);

    if(bam_set1(newRead,strlen(bam_get_qname(read)),bam_get_qname(read),read->core.flag,read->core.tid,pos,read->core.qual,newCigarLen,newCigar,read->core.mtid,read->core.mpos,read->core.isize,newQueryLen,newSeq,newQual,l_aux)<0)
    {
        bam_destroy1(newRead);
        throw runtime_error("Error: Could not initialize new read");
    }

    memcpy(bam_get_aux(newRead),bam_get_aux(read),l_aux); newRead->l_data+=l_aux;
    bam_destroy1(read); *pRead=newRead;
}
//----------------------------------------------------------------
void SplitNCigar::TrimEnd(bam1_t ** pRead,hts_pos_t endPos)
{
    bam1_t * read=*pRead;
    hts_pos_t pos=read->core.pos;
    uint32_t cigarLen=read->core.n_cigar;
    uint32_t * cigar=bam_get_cigar(read);
    if(bam_cigar_op(cigar[cigarLen-1])==BAM_CSOFT_CLIP) endPos+=bam_cigar_oplen(cigar[cigarLen-1]);

    hts_pos_t junctionPos; for(auto junction=intronBegin.upper_bound(pos);;junction++)
    {
        if(junction==intronBegin.end() || junction->first>endPos) return;
        if(junction->second>=nSupportReads) {junctionPos=junction->first; break;}
    }

    size_t queryIndex=0;
    uint8_t * seq=(uint8_t*)bam_get_seq(read);
    char * qual=(char*)bam_get_qual(read);

    size_t newCigarLen=0;
    uint32_t newCigar[256];

    size_t newQueryLen=0;
    char newSeq[2048]; char newQual[2048];

    size_t nBasesCovered=0;
    size_t nBasesToClip=0;
    size_t nMismatches=0;

    uint32_t i=0; for(;i<cigarLen;i++)
    {
        uint32_t op=bam_cigar_op(cigar[i]);
        uint32_t opLen=bam_cigar_oplen(cigar[i]); uint32_t remainder=opLen;

        switch(op)
        {
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:

            for(;remainder>0 && pos<junctionPos;remainder--,pos++,queryIndex++,newQueryLen++)
            {
                newSeq[newQueryLen]=seq_nt16_str[bam_seqi(seq,queryIndex)]; newQual[newQueryLen]=qual[queryIndex];
            }

            if(pos<junctionPos){newCigar[newCigarLen++]=cigar[i]; continue;}
            if(remainder<opLen) newCigar[newCigarLen++]=bam_cigar_gen(opLen-remainder,op);

            for(;remainder>0;remainder--,pos++,queryIndex++)
            {
                nBasesCovered++; nBasesToClip++; nMismatches+=seq_nt16_table[size_t(referenceSeq[pos])]!=bam_seqi(seq,queryIndex);
            }

            goto CALC_REMAINDER;

        case BAM_CINS:

            newCigar[newCigarLen++]=cigar[i];

            for(;opLen>0;opLen--,queryIndex++,newQueryLen++)
            {
                newSeq[newQueryLen]=seq_nt16_str[bam_seqi(seq,queryIndex)];
                newQual[newQueryLen]=qual[queryIndex];
            }

            continue;

        case BAM_CDEL:

            if(pos<=junctionPos && junctionPos<=pos+opLen) {pos+=opLen; goto CALC_REMAINDER;}   //Also skip when the read ends with del hence junctionPos<=pos+opLen
            newCigar[newCigarLen++]=cigar[i];
            pos+=opLen;
            continue;

        case BAM_CSOFT_CLIP:

            if(pos+opLen>=junctionPos) return;

            newCigar[newCigarLen++]=cigar[i];

            for(;opLen>0;opLen--,queryIndex++,newQueryLen++)
            {
                newSeq[newQueryLen]=seq_nt16_str[bam_seqi(seq,queryIndex)]; newQual[newQueryLen]=qual[queryIndex];
            }

            continue;

        case BAM_CHARD_CLIP:

            newCigar[newCigarLen++]=cigar[i];
            continue;
        }
    }
CALC_REMAINDER:

    for(i++;i<cigarLen;i++)
    {
        uint32_t opLen=bam_cigar_oplen(cigar[i]);

        switch(bam_cigar_op(cigar[i]))
        {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
            case BAM_CSOFT_CLIP:

                for(;opLen>0;opLen--,pos++,queryIndex++)
                {
                    nBasesCovered++;nBasesToClip++; nMismatches+=seq_nt16_table[size_t(referenceSeq[pos])]!=bam_seqi(seq,queryIndex);
                }

                continue;

            case BAM_CINS:

                nBasesToClip+=opLen;
                continue;
        }
    }

    if(nBasesCovered>maxOverhang || nMismatches<minMismatches) return;
    if(nBasesToClip!=0) newCigar[newCigarLen++]=bam_cigar_gen(nBasesToClip,BAM_CHARD_CLIP);

    bam1_t * newRead=bam_init1();
    size_t l_aux=bam_get_l_aux(read);

    if(bam_set1(newRead,strlen(bam_get_qname(read)),bam_get_qname(read),read->core.flag,read->core.tid,read->core.pos,read->core.qual,newCigarLen,newCigar,read->core.mtid,read->core.mpos,read->core.isize,newQueryLen,newSeq,newQual,l_aux)<0)
    {
        bam_destroy1(newRead);
        throw runtime_error("Error: Could not initialize new read");
    }

    memcpy(bam_get_aux(newRead),bam_get_aux(read),l_aux); newRead->l_data+=l_aux;
    bam_destroy1(read); *pRead=newRead;
}
//----------------------------------------------------------------
void SplitNCigar::Run(size_t i)
{
    //----------------------------------------------------------------
    //First pass determine intron boundaries
    //----------------------------------------------------------------

    string targetName=inputBamFile.header->target_name[i];
    inputBamFile.SetRegion(targetName);

    bam1_t * read=inputBamFile.read;

    while(inputBamFile.ReadRegion()>=0)
    {
        hts_pos_t pos=read->core.pos;
        size_t cigarLen=size_t(read->core.n_cigar);
        uint32_t * cigar=bam_get_cigar(read);

        for(size_t i=0;i<cigarLen;i++)
        {
            switch(bam_cigar_op(cigar[i]))
            {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
            case BAM_CDEL:

                pos+=hts_pos_t(bam_cigar_oplen(cigar[i]));
                continue;

            case BAM_CREF_SKIP:

                intronBegin[pos]++;
                pos+=hts_pos_t(bam_cigar_oplen(cigar[i]));
                intronEnd[pos-1]++;
                continue;
            }
        }
    }

    //----------------------------------------------------------------
    //Second pass do the actual split and trim
    //----------------------------------------------------------------

    referenceSeq=fastaEntries.at(targetName).seq;
    inputBamFile.SetRegion(targetName);

    vector<bam1_t*> segments; segments.reserve(1024);   //Segment buffer

    while(inputBamFile.ReadRegion()>=0)
    {
        if(read->core.qual>maxAlignmentScore) read->core.qual=maxAlignmentScore;

        hts_pos_t startPos,pos; startPos=pos=read->core.pos;

        size_t cigarLen=size_t(read->core.n_cigar); if(cigarLen==0) continue;
        uint32_t * cigar=bam_get_cigar(read);

        size_t queryNameLen=strlen(bam_get_qname(read));
        const char * queryName=bam_get_qname(read);

        size_t queryIndex=0;
        uint8_t * seq=(uint8_t*)bam_get_seq(read);
        char * qual=(char*)bam_get_qual(read);

        size_t newCigarLen=0;
        uint32_t newCigar[256];

        size_t newQueryLen=0;
        char newSeq[2048];char newQual[2048];

        size_t l_aux=bam_get_l_aux(read);
        uint8_t * aux=bam_get_aux(read);

        int32_t sumOfSegments=0;

        bam1_t * segment;

        for(size_t i=0;i<cigarLen;i++)
        {
            uint32_t opLen=bam_cigar_oplen(cigar[i]);

            switch(bam_cigar_op(cigar[i]))
            {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:

                newCigar[newCigarLen++]=cigar[i];

                pos+=opLen;

                for(;opLen>0;opLen--,queryIndex++,newQueryLen++)
                {
                    newSeq[newQueryLen]=seq_nt16_str[bam_seqi(seq,queryIndex)]; newQual[newQueryLen]=qual[queryIndex];
                }

                continue;

            case BAM_CINS:
            case BAM_CSOFT_CLIP:

                newCigar[newCigarLen++]=cigar[i];

                for(;opLen>0;opLen--,queryIndex++,newQueryLen++)
                {
                    newSeq[newQueryLen]=seq_nt16_str[bam_seqi(seq,queryIndex)]; newQual[newQueryLen]=qual[queryIndex];
                }

                continue;

            case BAM_CDEL:

                newCigar[newCigarLen++]=cigar[i];
                pos+=opLen;
                continue;

            case BAM_CREF_SKIP:

                if(useClipping==true) newCigar[newCigarLen++]=bam_cigar_gen(opLen,BAM_CHARD_CLIP);

                segment=bam_init1();

                if(bam_set1(segment,queryNameLen,queryName,read->core.flag,read->core.tid,startPos,read->core.qual,newCigarLen,newCigar,read->core.mtid,read->core.mpos,read->core.isize,newQueryLen,newSeq,newQual,l_aux)<0)
                {
                    bam_destroy1(segment);
                    throw runtime_error("Error: Could not initialize new segment");
                }

                memcpy(bam_get_aux(segment),aux,l_aux); segment->l_data+=l_aux;

                if(useTrimming==true && startPos==read->core.pos) TrimBegin(&segment,pos-1);

                segments.push_back(segment); sumOfSegments+=segment->core.l_qseq;

                pos+=opLen; startPos=pos;

                newCigarLen=0;
                newQueryLen=0;

                if(useClipping==true) newCigar[newCigarLen++]=bam_cigar_gen(opLen,BAM_CHARD_CLIP);

                continue;
            }
        }

        segment=bam_init1();

        if(bam_set1(segment,queryNameLen,queryName,read->core.flag,read->core.tid,startPos,read->core.qual,newCigarLen,newCigar,read->core.mtid,read->core.mpos,read->core.isize,newQueryLen,newSeq,newQual,l_aux)<0)
        {
            bam_destroy1(segment);
            throw runtime_error("Error: Could not initialize new segment");
        }

        memcpy(bam_get_aux(segment),aux,l_aux); segment->l_data+=l_aux;

        if(useTrimming)
        {
            if(startPos==read->core.pos) TrimBegin(&segment,pos-1);
            TrimEnd(&segment,pos-1);
        }

        segments.push_back(segment); sumOfSegments+=segment->core.l_qseq;

        if(useClipping)
        {
            int32_t cumSum=0; for(size_t i=0,iEnd=segments.size()-1;i<iEnd;i++) //Correct clipping
            {
                segment=segments[i];
                cumSum+=segment->core.l_qseq;

                auto & right=bam_get_cigar(segment)[segment->core.n_cigar-1];
                right=bam_cigar_gen(bam_cigar_oplen(right)+sumOfSegments-cumSum,BAM_CHARD_CLIP);

                auto & left=bam_get_cigar(segments[i+1])[0];
                left=bam_cigar_gen(bam_cigar_oplen(left)+cumSum,BAM_CHARD_CLIP);
            }
        }

        for(size_t i=0,nSegments=segments.size();i<nSegments;i++)
        {
            segment=segments[i];
            bool ret=outputBamFile.Write(segment); bam_destroy1(segment);
            if(ret==false) throw runtime_error("Error: Could not write bam file");
        }

        segments.clear();
    }

    intronBegin.clear(); intronEnd.clear();
}
//----------------------------------------------------------------
#define VERSION "1.1"
//----------------------------------------------------------------
#define INPUT_BAM_FILE      'i'
#define OUTPUT_BAM_FILE     'o'
#define REFERENCE_FILE      'r'
#define MAX_ALIGNMENT_SCORE 's'
#define SUPPORT_READS       'R'
#define MIN_MISMATCHES      'm'
#define MAX_OVERHANG        'O'
#define THREADS             'p'
#define ENABLE_CLIPPING     'c'
#define ENABLE_TRIMMING     't'
#define HELP                'h'
#define SHORT_OPTIONS       "i:o:r:s:R:m:O:p:cth"
//----------------------------------------------------------------
struct option longOptions[] =
{
    {"input-bam-file",required_argument,nullptr,INPUT_BAM_FILE},
    {"output-bam-file",required_argument,nullptr,OUTPUT_BAM_FILE},
    {"reference-file",required_argument,nullptr,REFERENCE_FILE},
    {"max-alignment-score",required_argument,nullptr,MAX_ALIGNMENT_SCORE},
    {"support-reads",required_argument,nullptr,SUPPORT_READS},
    {"min-mismatches",required_argument,nullptr,MIN_MISMATCHES},
    {"max-overhang",required_argument,nullptr,MAX_OVERHANG},
    {"threads",required_argument,nullptr,THREADS},
    {"enable-clipping",no_argument,nullptr,ENABLE_CLIPPING},
    {"enable-trimming",no_argument,nullptr,ENABLE_TRIMMING},
    {"help",no_argument,nullptr,HELP}
};
//----------------------------------------------------------------
int SplitNCigar::Run(int argc,char *argv[])
{
    try
    {
        //----------------------------------------------------------------
        //Get input arguments
        //----------------------------------------------------------------

        bool showHelp=(argc==1);
        int option,optionIndex;

        bool useClipping=false;
        bool useTrimming=false;
        int maxAlignmentScore=255;
        size_t nSupportReads=2;
        size_t minMismatches=2;
        size_t maxOverhang=40;
        int nThreads=1;

        string command,inputBamFilename,outputPrefix,referenceFilename;

        command=string(argv[0]); for(int i=1;i<argc;i++) command+=string(" ")+string(argv[i]);

        while((option=getopt_long(argc,argv,SHORT_OPTIONS,longOptions,&optionIndex))>=0)
        {
            switch(option)
            {
            case INPUT_BAM_FILE: inputBamFilename=string(optarg); break;
            case OUTPUT_BAM_FILE: outputPrefix=string(optarg); outputPrefix=outputPrefix.substr(0,outputPrefix.find_last_of('.')); break;
            case REFERENCE_FILE: referenceFilename=string(optarg); break;
            case MAX_ALIGNMENT_SCORE: maxAlignmentScore=atoi(optarg); break;
            case SUPPORT_READS: nSupportReads=atoi(optarg); break;
            case MIN_MISMATCHES: minMismatches=atoi(optarg); break;
            case MAX_OVERHANG: maxOverhang=atoi(optarg); break;
            case THREADS: nThreads=atoi(optarg); break;
            case ENABLE_CLIPPING: useClipping=true; break;
            case ENABLE_TRIMMING: useTrimming=true; break;
            case HELP: showHelp=true; break;
            }
        }

        //----------------------------------------------------------------
        //Show help
        //----------------------------------------------------------------

        if(showHelp)
        {
            cerr << "split_N_cigar [options]"                                                                               << endl;
            cerr                                                                                                            << endl;
            cerr << "-i --input-bam-file <text>     Single input bam file (required)"                                       << endl;
            cerr << "-o --output-bam-file <text>    Output bam file (optional default _split.bam)"                          << endl;
            cerr << "-r --reference-file <text>     Reference file to remove overhang (required)"                           << endl;
            cerr << "-s --max-alignment-score <int> Maximum alignment score (optional default 255)"                         << endl;
            cerr << "-R --support-reads <int>       Number of supporting reads to use with trimming (optional default 2)"   << endl;
            cerr << "-m --min-mismatches <int>      Minimum number of mismatches to use with trimming (optional default 2)" << endl;
            cerr << "-O --max-overhang <int>        Maximum overhang in bases to use with trimming (optional default 40)"   << endl;
            cerr << "-p --threads <int>             Number of threads to use (optional default 1)"                          << endl;
            cerr << "-c --enable-clipping <void>    When introduced use hard clipping (optional default no clipping)"       << endl;
            cerr << "-t --enable-trimming <void>    When introduced use trimming (optional default not trimming)"           << endl;
            cerr << "-h --help <void>               This help"                                                              << endl;
            cerr                                                                                                            << endl;

            return 0;
        }

        //----------------------------------------------------------------
        //Check input arguments
        //----------------------------------------------------------------

        if(inputBamFilename.empty()) throw runtime_error("Error: Please specify an input bam file!");
        if(referenceFilename.empty()) throw runtime_error("Error: Please specify a reference fasta file!");

        if(outputPrefix.empty())outputPrefix=inputBamFilename.substr(0,inputBamFilename.find_last_of('.')) + string("_split");

        if(maxAlignmentScore<20) maxAlignmentScore=20;
        if(maxAlignmentScore>255) maxAlignmentScore=255;

        if(nThreads<1) nThreads=1;

        //----------------------------------------------------------------
        //Open fasta file
        //----------------------------------------------------------------

        FastaFile fastaFile(referenceFilename);

        size_t nTargets;
        string targetFiles;
        {
            //----------------------------------------------------------------
            //Test open input bam file
            //----------------------------------------------------------------

            BamFile bamFile(inputBamFilename);

            //----------------------------------------------------------------
            //Get number of targets and check if they exist in the fasta file
            //----------------------------------------------------------------

            string targetName;

            nTargets=bamFile.header->n_targets;
            for(size_t i=0;i<nTargets;i++)
            {
                targetName=bamFile.header->target_name[i];
                if(fastaFile.entries.count(targetName)==0) throw runtime_error(string("Error: Target name not in fasta file: ")+targetName);

                 targetFiles+=' '+outputPrefix+'_'+bamFile.header->target_name[i]+".bam";
            }

            targetFiles+=' '+outputPrefix+"_unmapped.bam";
        }

        //----------------------------------------------------------------
        //Start parallel code
        //----------------------------------------------------------------

        omp_set_num_threads(nThreads);

        #pragma omp parallel
        {
            SplitNCigar splitNCigar(useClipping,useTrimming,maxAlignmentScore,nSupportReads,minMismatches,maxOverhang,fastaFile.entries,inputBamFilename);   //Create splitNCigar object in each worker thread

            if(nThreads==1) splitNCigar.outputBamFile.Create(outputPrefix+".bam",splitNCigar.inputBamFile.header,"split_N_cigar",VERSION,command);

            #pragma omp for schedule(dynamic)
            for(size_t i=0;i<nTargets;i++)
            {
                if(nThreads>1) splitNCigar.outputBamFile.Create(outputPrefix+"_"+splitNCigar.inputBamFile.header->target_name[i]+".bam",splitNCigar.inputBamFile.header,"split_N_cigar",VERSION,command);

                splitNCigar.Run(i);
            }

            #pragma omp single
            {
                splitNCigar.inputBamFile.SetRegion("*");

                if(nThreads>1) splitNCigar.outputBamFile.Create(outputPrefix+"_unmapped.bam",splitNCigar.inputBamFile.header,"split_N_cigar",VERSION,command);

                bam1_t * read=splitNCigar.inputBamFile.read;
                while(splitNCigar.inputBamFile.ReadRegion()>=0) splitNCigar.outputBamFile.Write(read);
            }
        }

        //----------------------------------------------------------------
        //Merge all files and sort
        //----------------------------------------------------------------

        {
            char command[1048576];

            if(nThreads>1)
            {
                snprintf(command,1048576,"samtools cat --no-PG -o %s.bam%s",outputPrefix.c_str(),targetFiles.c_str()); system(command);
                snprintf(command,1048576,"rm%s",targetFiles.c_str()); system(command);
            }

            snprintf(command,1048576,"samtools sort -o %s.bam##idx##%s.bam.bai --threads %i --write-index %s.bam",outputPrefix.c_str(),outputPrefix.c_str(),nThreads,outputPrefix.c_str()); system(command);
        }
    }
    catch(const runtime_error & error)
    {
        cerr << error.what() << endl;
        return 1;
    }

    return 0;
}
//----------------------------------------------------------------

