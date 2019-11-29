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
	Channel.fromFilePairs("${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz").set{ch_seqtk_mp}
	Channel.fromFilePairs("${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz").set{ch_seqtk_pe}
	
	// The directory for the sub_sampled raw reads to be output
	params.mp_dir_sub_raw_reads = "${params.mp_wkd}/raw_reads" 	
	params.pe_dir_sub_raw_reads = "${params.pe_wkd}/raw_reads" 	
	
	// Now do the subsampling using seqtk
	// The mate pair libraries
	process sub_sample_mp{
		tag "${base}"
		publishDir 	path: "${params.mp_wkd}/raw_reads"

		input:
		tuple val(base), file(reads) from ch_seqtk_mp

		output:
		file "*_sub.fastq.gz" into ch_fastqc_mp
		tuple val(base), file("*.fastq.gz") into ch_trim_mp

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

	// The paired end libraries
	process sub_sample_pe{
		tag "${base}"
		publishDir 	path: "${params.pe_wkd}/raw_reads"

		input:
		tuple val(base), file(reads) from ch_seqtk_pe

		output:
		file "*_sub.fastq.gz" into ch_fastqc_pe
		tuple val(base), file("*.fastq.gz") into ch_trim_pe

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
	Channel.fromPath("${params.pe_dir}/raw_reads/*.fastq.gz").set{ch_fastqc_pe}
	Channel.fromPath("${params.mp_dir}/raw_reads/*.fastq.gz").set{ch_fastqc_mp}
	Channel.fromFilePairs("${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz").set{ch_trim_pe}
	Channel.fromFilePairs("${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz").set{ch_trim_mp}
}


// Now to the pre_trim fastqc for each of the subsampled fastq.gz files
process fastqc_pre_trim_mp{
	tag "${read}"
	publishDir 	path: "${params.mp_wkd}/fastqc_pre_trim", mode: 'copy',
		saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

	input:
	file read from ch_fastqc_mp.flatten()

	output:
	file "*_fastqc.{zip,html}" into ch_fastqc_pre_trim_results_mp

	script:
	"""
	fastqc $read
	"""
}

// Now to the pre_trim fastqc for each of the subsampled fastq.gz files
process fastqc_pre_trim_pe{
	tag "${read}"
	publishDir 	path: "${params.pe_wkd}/fastqc_pre_trim", mode: 'copy',
		saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

	input:
	file read from ch_fastqc_pe.flatten()

	output:
	file "*_fastqc.{zip,html}" into ch_fastqc_pre_trim_results_pe

	script:
	"""
	fastqc $read
	"""
}

// //Now trim the reads using trimmomatic
// // At this point we will have a set of mp/pe processes for if we are working with sub sample
// // vs is we are not. This is beacause I wasn't able to work out how to get the 
// // output of the subsampling back into the same tuple format that is created using
// // the fromFilePairs read in. After this trimming, the output files should be the same
// // irrespective of whether we are working with subsampled or not. Therefore we will
// // not need to have two sets of processess as we move forwards.

// process trim_reads_mp{
// 	tag "${reads[0].getName()}"
// 	publishDir 	path: "${params.mp_wkd}/trimmed", mode: 'copy'

// 	input:
// 	tuple val(base), file(reads) from ch_trim_mp
	
// 	output:
// 	// Output that will be used for the post_trim fastqc
// 	// It is a flast list of all of the trimmed files
// 	file "*.fq.gz" into ch_fastqc_post_trim_input_mp
// 	// Output that will be used for the error_correction
// 	// This is a list of tuples that are the 1P and 2P output files only
// 	tuple file("*1P.fq.gz"), file("*2P.fq.gz") into ch_r_correct_input_mp

// 	script:
// 	outbase = reads[0].getName().replaceAll("_sub", "").replaceAll('_R1.fastq.gz', '.trimmed.fq.gz')
// 	"""
// 	trimmomatic PE -threads ${params.trimmomatic_threads} -basein ${reads[0]} \\
// 		-baseout $outbase \\
// 		ILLUMINACLIP:${params.tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
// 		LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11
// 	"""
// }

// process trim_reads_pe{
// 	tag "${fwd_rreads[0].getName()ead}"
// 	publishDir 	path: "${params.pe_wkd}/trimmed", mode: 'copy'

// 	input:
// 	tuple val(base), file(reads) from ch_trim_pe
	
// 	output:
// 	file "*.fq.gz" into ch_fastqc_post_trim_input_pe
// 	tuple file("*1P.fq.gz"), file("*2P.fq.gz") into ch_r_correct_input_pe

// 	script:
// 	//At this point we can get rid of the _sub part of the name if it exists
// 	// This way the file names will be the same from here on
// 	// whether we are working with subsampled or not
// 	outbase = fwd_read.getName().replaceAll("_sub", "").replaceAll('_R1.fastq.gz', '.trimmed.fq.gz')
// 	"""
// 	trimmomatic PE -threads ${params.trimmomatic_threads} -basein $fwd_read \\
// 		-baseout $outbase \\
// 		ILLUMINACLIP:${params.tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
// 		LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11
// 	"""
// }


ch_trim_pe.subscribe { println "ch_trim_pe file: $it" }
ch_trim_mp.subscribe { println "ch_trim_mp file: $it" }


// // // results.subscribe { println "value: $it" }
// // Now to the pre_trim fastqc for each of the subsampled fastq.gz files
// process fastqc_post_trim_mp{
// 	tag "${read}"
// 	publishDir 	path: "${params.mp_wkd}/fastqc_post_trim", mode: 'copy',
// 		saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

// 	input:
// 	file read from ch_fastqc_post_trim_input_mp.flatten()

// 	output:
// 	file "*_fastqc.{zip,html}" into ch_fastqc_post_trim_results_mp

// 	script:
// 	"""
// 	fastqc $read
// 	"""
// }

// // Now to the pre_trim fastqc for each of the subsampled fastq.gz files
// process fastqc_post_trim_pe{
// 	tag "${read}"
// 	publishDir 	path: "${params.pe_wkd}/fastqc_post_trim", mode: 'copy',
// 		saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

// 	input:
// 	file read from ch_fastqc_post_trim_input_pe.flatten()

// 	output:
// 	file "*_fastqc.{zip,html}" into ch_fastqc_post_trim_results_pe

// 	script:
// 	"""
// 	fastqc $read
// 	"""
// }

// // Error correction of the read pairs with rcorrector
// // Attemp having two read in channels one which we get the val from and one which we get a 
// // tuple from. 
// // Then make a list that is then put ina channel
// process attempt_list{
// 	tag "${reads[0]}"
// }

// ch_fastqc_post_trim_results_pe.subscribe { println "ch_fastqc_post_trim_results_pe file: $it" }
// ch_fastqc_post_trim_results_mp.subscribe { println "ch_fastqc_post_trim_results_mp file: $it" }