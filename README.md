# split_and_trim
Simple program to split and trim spliced reads\
\
\
Dependencies:\
\
HTSLib should be installed in the same directory in the source code and compiled to genereate\
\
The program also requires additional system development libraries.\
On Debian/Ubuntu libbz2-dev, liblzma-dev, libcurl4-openssl-dev, libssl-dev, libdeflate-dev\
\
To build use:\
cmake -DCMAKE_BUILD_TYPE=Release -S ./\
make\
\
Usage:\
\
split_and_trim [options]\
|Short options|Long option|Type|Description|
|---|---|---|---|
|-i|--input-bam-file      |<text> |Single input bam file (required)|
|-o|--output-bam-file     |<text> |Output bam file (optional default _split.bam)|
|-r|--reference-file      |<text> |Reference file to remove overhang (required)|
|-s|--max-alignment-score |<int>  |Maximum alignment score (optional default 255)|
|-R|--support-reads       |<int>  |Number of supporting reads to use with trimming (optional default 2)|
|-m|--min-mismatches      |<int>  |Minimum number of mismatches to use with trimming (optional default 2)|
|-O|--max-overhang        |<int>  |Maximum overhang in bases to use with trimming (optional default 40)|
|-p|--threads             |<int>  |Number of threads to use (optional default 1)|
|-c|--enable-clipping     |<void> |When introduced use hard clipping (optional default no clipping)|
|-t|--enable-trimming     |<void> |When introduced use trimming (optional default not trimming)|
|-h|--help                |<void> |This help|
