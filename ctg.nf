#!/usr/bin/env nextflow

// The git and base directory
params.basedir = "/home/humebc/projects/st_genome"

// The directories for each the mate pair libraries and the paired end libraries
params.mp_dir = "${params.basedir}/mate_pair_libs"
params.pe_dir = "${params.basedir}/paired_end_reads"

// The sub sampled directory if subsampling will be used
params.mp_dir_sub = "${params.basedir}/mate_pair_libs_sub"
params.pe_dir_sub = "${params.basedir}/paired_end_reads_sub"


// We don't want to have a set of processess for sub_sampled vs not subsampled
// So we will set a working output base dir here so that the same processes can be used
if (params.sub_sample){
	params.mp_wkd = "${params.basedir}/mate_pair_libs_sub"
	params.pe_wkd = "${params.basedir}/paired_end_reads_sub"
}else{
	params.mp_wkd = "${params.basedir}/mate_pair_libs"
	params.pe_wkd = "${params.basedir}/paired_end_reads"
}

// We will enable subsampling for development purposes.
// We will subsample all of the sequencing reads down to 1 000 000 reads. 
if (params.sub_sample){
	// The full size raw files that will be used as input to the seqtk
	// We will put all of the sequencing files into a single channel
	Channel.fromFilePairs(["${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz", "${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz"]).set{ch_sub_sample_seqtk_input}
	// Channel.fromFilePairs(["${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz"]).set{ch_sub_sample_seqtk_input}
	
	// The directory for the sub_sampled raw reads to be output
	params.mp_dir_sub_raw_reads = "${params.mp_wkd}/raw_reads" 	
	params.pe_dir_sub_raw_reads = "${params.pe_wkd}/raw_reads" 	
	
	// Now do the subsampling using seqtk
	// Once seqtk is complete then will need to do the gz compression
	// Then will need to put this gz compressed file into the ch_fastqc_pr channel for the fastqc
	// and the ch_trim_input for the trimming
	process sub_sample_seqtk{
		tag "${base}"
		conda "envs/general_conda_env.yaml"
		// publishDir = [
		// 	[path: "${params.mp_dir_sub_raw_reads}", mode: 'copy', overwrite: 'true', pattern: "M_18*.fastq.gz"],
		// 	[path: "${params.pe_dir_sub_raw_reads}", mode: 'copy', overwrite: 'true', pattern: "M_17*.fastq.gz"]
		// ]

		input:
		tuple val(base), file(reads) from ch_sub_sample_seqtk_input

		output:
		file "*_sub.fastq.gz" into ch_fastqc_pre_trim_input
		tuple val(base), file("*.fastq.gz") into ch_trim_reads_input

		script:
		read_out_one = reads[0].getName().replaceAll(".fastq.gz", "_sub.fastq")
		read_out_two = reads[1].getName().replaceAll(".fastq.gz", "_sub.fastq")
		
		"""
		seqtk sample -s100 ${reads[0]} 1000000 > ${read_out_one}
		gzip ${read_out_one}
		seqtk sample -s100 ${reads[1]} 1000000 > ${read_out_two}
		gzip ${read_out_two}
		"""
	}
	
}else{
	Channel.fromPath(["${params.pe_dir}/raw_reads/*.fastq.gz", "${params.mp_dir}/raw_reads/*.fastq.gz"]).set{ch_fastqc_pre_trim_input}
	Channel.fromFilePairs(["${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz", "${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz"]).set{ch_trim_reads_input}
}


// Now trim the reads using trimmomatic
process trim_reads{
	tag "${reads[0].getName()}"
	conda "envs/general_conda_env.yaml"
	
	// I don't think we need to publish these
	// publishDir = [
	// 		[path: "${params.mp_wkd}/trimmed", mode: 'copy', overwrite: 'true', pattern: "M_18*"],
	// 		[path: "${params.pe_wkd}/trimmed", mode: 'copy', overwrite: 'true', pattern: "M_17*"]
	// 	]

	input:
	tuple val(base), file(reads) from ch_trim_reads_input
	
	output:
	// Output that will be used for the post_trim fastqc
	// It is a list of all of the trimmed files that will be flattened in the fastqz_post_trim process
	file "*.fq.gz" into ch_fastqc_post_trim_input
	// Output that will be used for the error_correction
	// This is a list of tuples that are the 1P and 2P output files only
	tuple file("*1P.fq.gz"), file("*2P.fq.gz") into ch_rcorrector_input

	script:
	// We will modify the name thusly, so that regardless of whether we started with subsampling
	// the names will be inagreement.
	outbase = reads[0].getName().replaceAll("_sub", "").replaceAll('_R1.fastq.gz', '.trimmed.fq.gz')
	"""
	trimmomatic PE -threads ${params.trimmomatic_threads} -basein ${reads[0]} \\
		-baseout $outbase \\
		ILLUMINACLIP:${params.tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
		LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11
	"""
}

if (params.make_fastqc){
	// Now to the pre_trim fastqc for each of the subsampled fastq.gz files
	process fastqc_pre_trim{
		tag "${read}"
		conda "envs/general_conda_env.yaml"
		publishDir = [
				[path: "${params.mp_wkd}/fastqc_pre_trim", mode: 'copy', overwrite: 'true', saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}, pattern: "M_18*"],
				[path: "${params.pe_wkd}/fastqc_pre_trim", mode: 'copy', overwrite: 'true', saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}, pattern: "M_17*"]
			]

		input:
		file read from ch_fastqc_pre_trim_input.flatten()

		output:
		file "*_fastqc.{zip,html}" into ch_fastqc_pre_trim_output

		script:
		"""
		fastqc $read
		"""
	}
	
	// Now to the pre_trim fastqc for each of the subsampled fastq.gz files
	process fastqc_post_trim{
		tag "${read}"
		conda "envs/general_conda_env.yaml"
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
}


// Now error correction
// The output here will be a little difficult because we need to separate the mp and pe
// reads into seperate channels. We may need to write a mapping function for this
process rcorrector{
	tag "${trimmed_read_one}"
	conda "envs/general_conda_env.yaml"
	cpus params.rcorrector_threads

	input:
	tuple file(trimmed_read_one), file(trimmed_read_two) from ch_rcorrector_input

	output:
	tuple file("*1P.cor.fq.gz"), file("*2P.cor.fq.gz") into ch_mp_out, ch_pe_out

	script:
	"""
	run_rcorrector.pl -1 $trimmed_read_one -2 $trimmed_read_two -od . -t ${task.cpus}
	"""
}

// ch_mp_out.toList().view()


// Each of the ouput channels from rcorrector has both the paired end
// and the mate pairs file sets in them.
// We want to remove the seq types that are not required respectively 
// from each channel. To do this we will write a map function
// to parse through the tuples and remove those tuples that don't contain the 
// respective M_18 or M_17 string
// Whoo! This works!
ch_mp_out.toList().flatMap{
		// Create the new list to return
		List output_list = new ArrayList();
		// for each tuple in the list
		for (i=0; i<(it.size()); i++){
			println("the get name function gives us ${it[i][0].getName()}")
			if (it[i][0].getName().contains("M_18")){
				// Then this is an mp tuple
				// add it to the output_list
				println("Tuple ${it[i]} contains M_18 so is an mp item. Adding to the output list")
				output_list.add(it[i])
				println("Output list now contains ${output_list}")
			}else{
				// Then this should be a pe tuple
				println("Tuple ${it[i]} does not contain M_18 so is a pe item.")
			}
		}
		println("We have finished iterating through the channel and are ready to output")
		println("output list looks like this:")
		println("${output_list}")
		return output_list
    }.view()


// // After error correction we will want to split the channels up by paired end and mate pair


// // The paired end reads will go into BBMerge, to merge as many as possible
// process bbMerge{
// 	tag "$rcorrected_read_one"
// 	conda "envs/general_conda_env.yaml"

// 	input:
// 	tuple file(rcorrected_read_one_pe), file(rcorrected_read_two_pe) from ch_pe_out

// 	output:
// 	tuple file(".merged.fastq.gz"), file("_1.merged.fastq.gz"), file("_2.merged.fastq.gz") into ch_bbmerge_out
// 	file "*ihist.txt" into ch_bbmerge_ihist_out

// 	script:
// 	out_merged_name = rcorrected_read_one_pe.getName().replaceAll("1P.cor.fq.gz", ".merged.fq.gz")
// 	one_unmerged_name = rcorrected_read_one_pe.getName().replaceAll(".fq.gz", ".unmerged.fq.gz")
// 	two_unmerged_name = rcorrected_read_two_pe.getName().replaceAll(".fq.gz", ".unmerged.fq.gz")
// 	ihist_name = rcorrected_read_one_pe.getName().replaceAll("1P.cor.fq.gz", ".ihist.txt")
// 	"""
// 	bbmerge-auto.sh in1=$rcorrected_read_one_pe in2=$rcorrected_read_two_pe out=$out_merged_name outu1=$one_unmerged_name outu2=$two_unmerged_name ihist=$ihist_name ecct extend2=20 iterations=5
// 	"""

// }

// The mate pair reads will go into NXtrim to characterise them

// process nxtrim{
// 	tag "$rcorrected_read_one"
// 	conda "envs/general_conda_env.yaml"

// 	input:
// 	tuple file(rcorrected_read_one_mp), file(rcorrected_read_two_mp) from ch_mp_out

// 	output:
// 	tuple file("*mp.fq.gz"), file("*pe.fq.gz"), file("*se.fq.gz"), file("*unknown.fq.gz") into ch_nxtrim_output

// 	script:
// 	sample_name_out = rcorrected_read_one_mp.getName().replaceAll("1P.cor.fq.gz", "")
// 	"""
// 	nxtrim -1 $rcorrected_read_one_mp -2 $rcorrected_read_two_mp -O $sample_name_out
// 	"""
// }

// So that we don't lose the ability to use the resume option in nextflow, we will start another conda environment here
// We'll call it stage_two.yaml.
// Eventually we will want to merge all of the environments together if possible.

// The liu et al 2018 paper that we were following used clc workbench to do an initial assembly to get accurate estimate
// of the insert lengths. However, this is commercial software and I'd rather stay clear.
// For their main assembly, they use both SPADES and ALLPATHS-LG. We could use one of these to make an initial assembly to check insert sizes.

// We also need to take into account the fact that we still have coral reads in our data. Let's dump those out before we do our first assembly.
// It looks like bbmap will be a good way to map the various reads to our p lutea genome.

