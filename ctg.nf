#!/usr/bin/env nextflow
params.basedir = "/home/humebc/projects/st_genome"

params.mp_dir = "${params.basedir}/mate_pair_libs"
params.pe_dir = "${params.basedir}/paired_end_reads"
// The sub sampled directory if subsampling will be used
params.mp_dir_sub = "${params.basedir}/mate_pair_libs_sub"
params.pe_dir_sub = "${params.basedir}/paired_end_reads_sub"
// We don't want to have a set of processess for sub_sampled vs not subsampled
// So we will set a working output base dir here
if (params.sub_sample){
	params.mp_wkd = "${params.basedir}/mate_pair_libs_sub"
	params.pe_wkd = "${params.basedir}/paired_end_reads_sub"
}else{
	params.mp_wkd = "${params.basedir}/mate_pair_libs"
	params.pe_wkd = "${params.basedir}/paired_end_reads"
}

// If subsampling
if (params.sub_sample){
	// The full size raw files that will be used as input to the seqtk
	// Channel.fromFilePairs("${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz").set{ch_seqtk_mp}
	Channel.fromFilePairs(["${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz", "${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz"]).set{ch_seqtk_input}
	
	
	// The directory for the sub_sampled raw reads to be output
	params.mp_dir_sub_raw_reads = "${params.mp_wkd}/raw_reads" 	
	params.pe_dir_sub_raw_reads = "${params.pe_wkd}/raw_reads" 	
	
	// Now do the subsampling using seqtk
	// The mate pair libraries
	process sub_sample{
		tag "${base}"
		publishDir = [
			[path: "${params.mp_dir_sub_raw_reads}", mode: 'copy', overwrite: 'true', pattern: "M_18*.fastq.gz"],
			[path: "${params.pe_dir_sub_raw_reads}", mode: 'copy', overwrite: 'true', pattern: "M_17*.fastq.gz"]
		]

		input:
		tuple val(base), file(reads) from ch_seqtk_input

		output:
		file "*_sub.fastq.gz" into ch_fastqc_input
		tuple val(base), file("*.fastq.gz") into ch_trim_input

		script:
		read_out_one = reads[0].getName().replaceAll(".fastq.gz", "_sub.fastq")
		read_out_two = reads[1].getName().replaceAll(".fastq.gz", "_sub.fastq")
		
		"""
		seqtk sample -s100 ${reads[0]} 10000 > ${read_out_one}
		gzip ${read_out_one}
		seqtk sample -s100 ${reads[1]} 10000 > ${read_out_two}
		gzip ${read_out_two}
		"""
	}
	// Once seqtk is complete then will need to do the gz compression
	// Then will need to put this gz compressed file into the ch_fastqc_pr channel for the fastqc
}else{
	Channel.fromPath(["${params.pe_dir}/raw_reads/*.fastq.gz", "${params.mp_dir}/raw_reads/*.fastq.gz"]).set{ch_fastqc_input}
	Channel.fromFilePairs(["${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz", "${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz"]).set{ch_trim_input}
}


// Now to the pre_trim fastqc for each of the subsampled fastq.gz files
process fastqc_pre_trim{
	tag "${read}"
	publishDir = [
			[path: "${params.mp_wkd}/fastqc_pre_trim", mode: 'copy', overwrite: 'true', saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}, pattern: "M_18*"],
			[path: "${params.pe_wkd}/fastqc_pre_trim", mode: 'copy', overwrite: 'true', saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}, pattern: "M_17*"]
		]

	input:
	file read from ch_fastqc_input.flatten()

	output:
	file "*_fastqc.{zip,html}" into ch_fastqc_pre_trim_output

	script:
	"""
	fastqc $read
	"""
}

// Now trim the reads using trimmomatic
process trim_reads_mp{
	tag "${reads[0].getName()}"
	publishDir = [
			[path: "${params.mp_wkd}/trimmed", mode: 'copy', overwrite: 'true', pattern: "M_18*"],
			[path: "${params.pe_wkd}/trimmed", mode: 'copy', overwrite: 'true', pattern: "M_17*"]
		]

	input:
	tuple val(base), file(reads) from ch_trim_input
	
	output:
	// Output that will be used for the post_trim fastqc
	// It is a flast list of all of the trimmed files
	file "*.fq.gz" into ch_fastqc_post_trim_input
	// Output that will be used for the error_correction
	// This is a list of tuples that are the 1P and 2P output files only
	tuple file("*1P.fq.gz"), file("*2P.fq.gz") into ch_r_correct_input

	script:
	outbase = reads[0].getName().replaceAll("_sub", "").replaceAll('_R1.fastq.gz', '.trimmed.fq.gz')
	"""
	trimmomatic PE -threads ${params.trimmomatic_threads} -basein ${reads[0]} \\
		-baseout $outbase \\
		ILLUMINACLIP:${params.tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
		LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11
	"""
}

// // results.subscribe { println "value: $it" }
// Now to the pre_trim fastqc for each of the subsampled fastq.gz files
process fastqc_post_trim{
	tag "${read}"
	publishDir = [
			[path: "${params.mp_wkd}/fastqc_post_trim", mode: 'copy', overwrite: 'true', saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}, pattern: "M_18*"],
			[path: "${params.pe_wkd}/fastqc_post_trim", mode: 'copy', overwrite: 'true', saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}, pattern: "M_17*"]
		]

	input:
	file read from ch_fastqc_post_trim_input.flatten()

	output:
	file "*_fastqc.{zip,html}" into ch_fastqc_post_trim_output

	script:
	"""
	fastqc $read
	"""
}

// // Error correction of the read pairs with rcorrector
// // Attemp having two read in channels one which we get the val from and one which we get a 
// // tuple from. 
// // Then make a list that is then put ina channel
// process attempt_list{
// 	tag "${reads[0]}"
// }

// ch_fastqc_post_trim_results_pe.subscribe { println "ch_fastqc_post_trim_results_pe file: $it" }
// ch_fastqc_post_trim_results_mp.subscribe { println "ch_fastqc_post_trim_results_mp file: $it" }