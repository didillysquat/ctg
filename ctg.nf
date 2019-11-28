#!/usr/bin/env nextflow
params.basedir = "/home/humebc/projects/st_genome"

params.mp_dir = "${params.basedir}/mate_pair_libs"
params.pe_dir = "${params.basedir}/paired_end_reads"
// The sub sampled directory if subsampling will be used
params.mp_dir_sub = "${params.basedir}/mate_pair_libs_sub"
params.pe_dir_sub = "${params.basedir}/paired_end_reads_sub"


if (params.sub_sample){
	// The full size raw files that will be used as input to the seqtk
	// Channel.fromFilePairs("${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz").set{ch_seqtk_mp}
	Channel.fromFilePairs("${params.mp_dir}/raw_reads/M_18_1922_HUME-ST-35-45_AD002_Lane2_{R1,R2}.fastq.gz").set{ch_seqtk_mp}
	Channel.fromFilePairs("${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz").set{ch_seqtk_pe}
	
	// The directory for the sub_sampled raw reads to be output
	params.mp_dir_sub_raw_reads = "${params.mp_dir_sub}/raw_reads" 	
	params.pe_dir_sub_raw_reads = "${params.pe_dir_sub}/raw_reads" 	
	
	// Create new file objects to check their existence and make if not exist
	def mp_raw_read_dir = new File("${params.mp_dir_sub_raw_reads}")
	def pe_raw_read_dir = new File("${params.pe_dir_sub_raw_reads}")
	if (!mp_raw_read_dir.exists()){
		mp_raw_read_dir.mkdirs()
	}
	if (!pe_raw_read_dir.exists()){
		pe_raw_read_dir.mkdirs()
	}
	
	// Now do the subsampling using seqtk
	process sub_sample_mp{
		echo true
		tag "${seq_base}"
		publishDir 	path: "${params.mp_dir_sub}/raw_reads"

		input:
		tuple val(seq_base), file(reads) from ch_seqtk_mp

		output:
		file "*_sub.fastq.gz" into results

		script:
		
		read_0_out = reads[0].getName().replaceAll(".fastq.gz", "_sub.fastq")
		read_1_out = reads[1].getName().replaceAll(".fastq.gz", "_sub.fastq")
		
		"""
		seqtk sample -s100 ${reads[0]} 10000 > ${read_0_out}
		gzip ${read_0_out}
		seqtk sample -s100 ${reads[1]} 10000 > ${read_1_out}
		gzip ${read_1_out}
		"""
		
	}
	// Once seqtk is complete then will need to do the gz compression
	// Then will need to put this gz compressed file into the ch_fastqc_pr channel for the fastqc
}else{
	ch_fastqc_pe = Channel.fromPath("${params.mp_dir}/raw_reads/*.fastq.gz")
}

results.subscribe { println "value: $it" }
