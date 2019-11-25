#!/usr/bin/env nextflow

params.mp_lib_path = "/home/humebc/projects/st_genome/mate_pair_libs/raw_reads"
params.pe_lib_path = "/home/humebc/projects/st_genome/paired_end_reads/raw_reads"

/* The channel that will contain the files mate pair and paired end raw seq 
files that will be consumed by the make_fastqc_pre_trimming process*/
mp_raw_read_for_fastqc_ch = Channel.fromPath(params.mp_lib_path)
pe_raw_read_for_fastqc_ch = Channel.fromPath(params.pe_lib_path)

process make_fastqc_pre_trim {
    input:
    file seq_file from mp_raw_read_for_fastqc_ch

    """
    echo $seq_file
    """
}
