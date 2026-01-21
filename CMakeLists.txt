cmake_minimum_required(VERSION 3.5)

include_directories(htslib/include)
link_directories(htslib/lib)

add_compile_options(-Wall -mavx2)

project(split_N_cigar LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(split_N_cigar main.cpp split_N_cigar.cpp split_N_cigar.h bam_file.cpp bam_file.h fasta_file.cpp fasta_file.h)

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
target_link_libraries (split_N_cigar PRIVATE libhts.a z m bz2 lzma curl crypto deflate Threads::Threads OpenMP::OpenMP_CXX)
