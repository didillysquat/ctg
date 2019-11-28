#!/usr/bin/env nextflow

params.mp_read_paths = "/home/humebc/projects/st_genome/mate_pair_libs/raw_reads/*fastq.gz"
params.pe_read_paths = "/home/humebc/projects/st_genome/paired_end_reads/raw_reads/*fastq.gz"
params.fastqc_output_dir_pre_trim = "/home/humebc/projects/st_genome/paired_end_reads/fastqc"

/* The channel that will contain the files mate pair and paired end raw seq 
files that will be consumed by the make_fastqc_pre_trimming process*/
mp_raw_read_for_fastqc_ch = Channel.fromPath(params.mp_read_paths)
mp_test = Channel.fromPath(params.mp_read_paths)
pe_raw_read_for_fastqc_ch = Channel.fromPath(params.pe_read_paths)
/*
process make_fastqc_pre_trim {
    echo true

    input:
    file seq_file from mp_raw_read_for_fastqc_ch

    """
    fastqc -o $params.fastqc_output_dir_pre_trim $seq_file 
    """
}
*/
process make_fastqc_pre_trim_test {
    echo true

    input:
    file seq_file from mp_test 
    file seq_file from pe_raw_read_for_fastqc_ch

    """
    echo $seq_file 
    """
}
