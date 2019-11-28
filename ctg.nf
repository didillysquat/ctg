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
	Channel.fromPath("${params.mp_dir}/raw_reads/*.fastq.gz").set{ch_seqtk_mp}
	Channel.fromPath("${params.pe_dir}/raw_reads/*.fastq.gz").set{ch_seqtk_pe}
	
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
	// The mate pair libraries
	process sub_sample_mp{
		tag "${read}"
		publishDir 	path: "${params.mp_dir_sub}/raw_reads"

		input:
		file read from ch_seqtk_mp

		output:
		file "*_sub.fastq.gz" into ch_fastqc_mp

		script:
		
		read_out = read.getName().replaceAll(".fastq.gz", "_sub.fastq")
		
		"""
		seqtk sample -s100 ${read} 10000 > ${read_out}
		gzip ${read_out}
		"""
	}

	// The paired end libraries
	process sub_sample_pe{
		tag "${read}"
		publishDir 	path: "${params.pe_dir_sub}/raw_reads"

		input:
		file read from ch_seqtk_pe

		output:
		file "*_sub.fastq.gz" into ch_fastqc_pe

		script:
		
		read_out = read.getName().replaceAll(".fastq.gz", "_sub.fastq")
		
		"""
		seqtk sample -s100 ${read} 10000 > ${read_out}
		gzip ${read_out}
		"""
	}
	// Once seqtk is complete then will need to do the gz compression
	// Then will need to put this gz compressed file into the ch_fastqc_pr channel for the fastqc
}else{
	ch_fastqc_pe = Channel.fromPath("${params.mp_dir}/raw_reads/*.fastq.gz")
	ch_fastqc_mp = Channel.fromPath("${params.mp_dir}/raw_reads/*.fastq.gz")
}

	// Now to the pre_trim fastqc for each of the subsampled fastq.gz files
	process fastqc_pre_trim_mp{
		tag "${read}"
		publishDir 	path: "${params.mp_dir_sub}/fastqc_pre_trim", mode: 'copy',
			saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

		input:
		file read from ch_fastqc_mp

		output:
    	file "*_fastqc.{zip,html}" into fastqc_pre_trim_results_mp

		script:
		"""
		fastqc $read
		"""
	}

	// Now to the pre_trim fastqc for each of the subsampled fastq.gz files
	process fastqc_pre_trim_pe{
		tag "${read}"
		publishDir 	path: "${params.pe_dir_sub}/fastqc_pre_trim", mode: 'copy',
			saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

		input:
		file read from ch_fastqc_pe

		output:
    	file "*_fastqc.{zip,html}" into fastqc_pre_trim_results_pe

		script:
		"""
		fastqc $read
		"""
	}
// results.subscribe { println "value: $it" }
