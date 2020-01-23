#!/usr/bin/env nextflow

/* The cacheing system that nextflow ships with by default is great when it works but it is causing a lot
of grief when it fails and this seems to happen quite often. Sometimes when looking at the dumped log files
it is not even possible to tell what the different hashes are referring to. As such, we are going to build
in a sort of caching to this nextflow doc. We will search to see if the files we need are already availale
and then we will work the first point at which the files are not available. In this way we won't have to 
worry about the cacheing mechanism failing. */

// The git and base directory
params.basedir = "/home/humebc/projects/st_genome"

// Directory with all of the python scripts
params.bin_dir = "${workflow.launchDir}/bin"

// The directories for each the mate pair libraries and the paired end libraries
params.mp_dir = "${params.basedir}/mate_pair_libs"
params.pe_dir = "${params.basedir}/paired_end_reads"

// The sub sampled directory if subsampling will be used
params.mp_dir_sub = "${params.basedir}/mate_pair_libs_sub"
params.pe_dir_sub = "${params.basedir}/paired_end_reads_sub"


// We don't want to have a set of processess for sub_sampled vs not subsampled
// So we will set a working output base dir here so that the same processes can be used
// Because subsampling takes so much time I am going to invest in a system here to check to see if the sub
// sampling has already happened.
// We will place subsampled sets of sequnces in their respective sub folders
// i.e. mate_pair_libs_sub and then within a foler names as the depth of the
// subsampling eg. mate_pair_libs_sub/10000
// This way we can simply check to see if the files already exist
// If they do then we can start a channel by reading in the files from these locations
// Else we will create them from scratch using the sub_sample_seqtk process below.
if (params.sub_sample){
	println("params.sub_sample is set to true")
	println("checking to see if the subsampled sequnce files already exist")
	params.mp_wkd = "${params.basedir}/mate_pair_libs_sub"
	params.pe_wkd = "${params.basedir}/paired_end_reads_sub"
	// Check to see if the directory already exists
	mp_dir = file("${params.mp_wkd}/raw_reads/${params.sub_sample_depth}")
	pe_dir = file("${params.pe_wkd}/raw_reads/${params.sub_sample_depth}")
	if ( mp_dir.exists() && pe_dir.exists() ){
		// Then the directories that will hold the subsampled sequences already exist
		// Check to see that they contain the correct number of folders
		if ( mp_dir.list().size() == 12 && pe_dir.list().size() == 4){
			// Then all the files are already here and so we can work from this dir
			params.sub_sample_from_scratch = false
			println("Successfully found subsampled seqs. We will not subsample from scratch")
		}else{
			params.sub_sample_from_scratch = true
			println("Unsuccessful in finding subsampled seqs at a depth of ${params.sub_sample_depth}. We will subsample from scratch")
		}
	}else{
		params.sub_sample_from_scratch = true
		println("Unsuccessful in finding subsampled seqs at a depth of ${params.sub_sample_depth}. We will subsample from scratch")
	}
}else{
	params.mp_wkd = "${params.basedir}/mate_pair_libs"
	params.pe_wkd = "${params.basedir}/paired_end_reads"
}

// We will enable subsampling for development purposes.
// We will subsample all of the sequencing reads down to the value that is set in the nextflow.config file reads. 

// We will use a different set of output directories for publishing depending on whether
// we are subsampling or not
if (params.sub_sample){
	// The directory for the sub_sampled raw reads to be output
	params.mp_dir_sub_raw_reads = "${params.mp_wkd}/raw_reads/${params.sub_sample_depth}" 	
	params.pe_dir_sub_raw_reads = "${params.pe_wkd}/raw_reads/${params.sub_sample_depth}"
	
	// The directory for the sub_sampled rimmed reads to be output
	params.mp_trimmed_reads_output_dir = "${params.mp_wkd}/trimmed_reads/${params.sub_sample_depth}" 	
	params.pe_trimmed_reads_output_dir = "${params.pe_wkd}/trimmed_reads/${params.sub_sample_depth}"

	// NB that the error correction will probably fail if the subsampled reads are low as there
	// is a required read depth required by tadpole.
	// The directory for the corrected reads to be output
	params.mp_corrected_reads_output_dir = "${params.mp_wkd}/corrected_reads/${params.sub_sample_depth}" 	
	params.pe_corrected_reads_output_dir = "${params.pe_wkd}/corrected_reads/${params.sub_sample_depth}"

	// The directory for the reads that are unmapped to lutea
	params.mp_unmapped_to_lutea_reads_output_dir = "${params.mp_wkd}/unmapped_to_lutea_reads/${params.sub_sample_depth}"
	params.pe_unmapped_to_lutea_reads_output_dir = "${params.pe_wkd}/unmapped_to_lutea_reads/${params.sub_sample_depth}"

	// The directory for the reads that have gone through nxtrim
	params.mp_nxtrim_output_dir = "${params.mp_wkd}/nxtrim_output/${params.sub_sample_depth}"

	// The directory for the discovar assembly output
	discovar_publish_path = "${params.pe_wkd}/discovardenovo_assembly/${params.sub_sample_depth}"
	
}else{
	// The directory for the trimmed reads to be output
	// When we are not subsampling then we don't need to have a separate directory
	// we just output the files directly into the trimmed read dir.
	params.mp_trimmed_reads_output_dir = "${params.mp_wkd}/trimmed_reads" 	
	params.pe_trimmed_reads_output_dir = "${params.pe_wkd}/trimmed_reads"

	// The directory for the corrected reads to be output
	params.mp_corrected_reads_output_dir = "${params.mp_wkd}/corrected_reads" 	
	params.pe_corrected_reads_output_dir = "${params.pe_wkd}/corrected_reads"

	// The directory for the reads that are unmapped to lutea
	params.mp_unmapped_to_lutea_reads_output_dir = "${params.mp_wkd}/unmapped_to_lutea_reads"
	params.pe_unmapped_to_lutea_reads_output_dir = "${params.pe_wkd}/unmapped_to_lutea_reads"

	// The directory for the reads that are unmapped to lutea
	params.mp_nxtrim_output_dir = "${params.mp_wkd}/nxtrim_output"

	// The directory for the discovar assembly output
	discovar_publish_path = "${params.pe_wkd}/discovardenovo_assembly"

	// The directory where we will store the size output information that will be required
	// by the assemblers
	mp_insert_size_info_output_dir = "${params.mp_wkd}/insert_size_info"
	pe_insert_size_info_output_dir = "${params.pe_wkd}/insert_size_info"
}


if (params.sub_sample){	
	if (params.sub_sample_from_scratch){
		println("subsampling")
		// The full size raw files that will be used as input to the seqtk
		// We will put all of the sequencing files into a single channel
		Channel.fromFilePairs(["${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz", "${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz"]).set{ch_sub_sample_seqtk_input}
		
		// Now do the subsampling using seqtk
		// Once seqtk is complete then will need to do the gz compression
		// Then will need to put this gz compressed file into the ch_fastqc_pr channel for the fastqc
		// and the ch_trim_input for the trimming
		process sub_sample_seqtk{
			cache 'lenient'
			tag "${base}"
			conda "envs/stage_one.yaml"
			publishDir = [
				[path: "${params.mp_dir_sub_raw_reads}", mode: 'copy', overwrite: 'true', pattern: "M_18*.fastq.gz"],
				[path: "${params.pe_dir_sub_raw_reads}", mode: 'copy', overwrite: 'true', pattern: "M_17*.fastq.gz"]
			]

			input:
			tuple val(base), file(reads) from ch_sub_sample_seqtk_input

			output:
			file "*_sub.fastq.gz" into ch_fastqc_pre_trim_input
			tuple val(base), file("*.fastq.gz") into ch_trim_reads_input

			script:
			read_out_one = reads[0].getName().replaceAll(".fastq.gz", "_sub.fastq")
			read_out_two = reads[1].getName().replaceAll(".fastq.gz", "_sub.fastq")
			
			"""
			seqtk sample -s100 ${reads[0]} ${params.sub_sample_depth} > ${read_out_one}
			gzip ${read_out_one}
			seqtk sample -s100 ${reads[1]} ${params.sub_sample_depth} > ${read_out_two}
			gzip ${read_out_two}
			"""
		}
	}else{
		// Then we generate the required channels directly from the already existing subsample directory
		// The directory for the sub_sampled raw reads to be output
		println("populating channels from already existant subsampled reads")
		Channel.fromPath(["${params.pe_dir_sub_raw_reads}/*.fastq.gz", "${params.mp_dir_sub_raw_reads}/*.fastq.gz"]).set{ch_fastqc_pre_trim_input}
		Channel.fromFilePairs(["${params.pe_dir_sub_raw_reads}/*{R1_sub,R2_sub}.fastq.gz", "${params.mp_dir_sub_raw_reads}/*{R1_sub,R2_sub}.fastq.gz"]).into{ch_trim_reads_input; ch_test}
	}
	
}else{
	println("populating channels from already existant raw reads")
	Channel.fromPath(["${params.pe_dir}/raw_reads/*.fastq.gz", "${params.mp_dir}/raw_reads/*.fastq.gz"]).set{ch_fastqc_pre_trim_input}
	Channel.fromFilePairs(["${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz", "${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz"]).set{ch_trim_reads_input}
}


// Now trim the reads using trimmomatic
// TODO check to see if the trimmed reads already exist
// If so, then scip the trimming and read in the trimmed reads directly from their containing directories.
mp_dir_as_file = file(params.mp_trimmed_reads_output_dir)
pe_dir_as_file = file(params.pe_trimmed_reads_output_dir)
if ( mp_dir_as_file.exists() && pe_dir_as_file.exists() ){
		// Then the directories that will hold the trimmed reads  already exist
		// Check to see that they contain the correct number of reads
		if ( mp_dir_as_file.list().size() == 24 && pe_dir_as_file.list().size() == 8){
			// Then all the files are already here and so we can work from this dir
			println("Successfully found trimmed seqs")
			params.trim_from_scratch = false
		}else{
			println("Could not find trimmed seqs")
			params.trim_from_scratch = true
		}
	}else{
		println("Could not find trimmed seqs")
		params.trim_from_scratch = true
	}

if (params.trim_from_scratch){
	println("Starting trimming")

	// Then we need to do the trimming and work from here.
	process trim_reads{
	cache 'lenient'
	tag "${reads[0].getName()}"
	conda "envs/stage_one.yaml"
	
	// These represent a considerable investment in time and so we should publish
	publishDir = [
			[path: params.mp_trimmed_reads_output_dir, mode: 'copy', overwrite: 'true', pattern: "M_18*"],
			[path: params.pe_trimmed_reads_output_dir, mode: 'copy', overwrite: 'true', pattern: "M_17*"]
		]

	input:
	tuple val(base), file(reads) from ch_trim_reads_input
	
	output:
	// Output that will be used for the post_trim fastqc
	// It is a list of all of the trimmed files that will be flattened in the fastqz_post_trim process
	file "*.fq.gz" into ch_fastqc_post_trim_input
	// Output that will be used for the error_correction
	// This is a list of tuples that are the 1P and 2P output files only
	tuple file("*1P.fq.gz"), file("*2P.fq.gz") into ch_tadpole_input

	script:
	// We will modify the name so that regardless of whether we started with subsampling
	// the names will be inagreement.
	paired_out_one = reads[0].getName().replaceAll("_sub", "").replaceAll('_R1.fastq.gz', '.trimmed.1P.fq.gz')
	unpaired_out_one = reads[0].getName().replaceAll("_sub", "").replaceAll('_R1.fastq.gz', '.trimmed.1U.fq.gz')
	paired_out_two = reads[1].getName().replaceAll("_sub", "").replaceAll('_R2.fastq.gz', '.trimmed.2P.fq.gz')
	unpaired_out_two = reads[1].getName().replaceAll("_sub", "").replaceAll('_R2.fastq.gz', '.trimmed.2U.fq.gz')
	"""
	trimmomatic PE -threads ${params.trimmomatic_threads} ${reads[0]} ${reads[1]} \\
		$paired_out_one $unpaired_out_one $paired_out_two $unpaired_out_two\\
		ILLUMINACLIP:${params.tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
		LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11
	"""
	}	
}else{
	println("Populating channels from already trimmed files")
	// Then the files have already been trimmed and we can create our channels directly from the directories

	// First we need to read in all of the files that are in each of the directories, join the lists
	// and then pass these into the ch_fastqc_post_trim_input list
	Channel.fromPath(["${params.mp_trimmed_reads_output_dir}/*.fq.gz", "${params.pe_trimmed_reads_output_dir}/*.fq.gz"]).set{ch_fastqc_post_trim_input}

	// Next we need to read in the paired files from each of the directories
	// We can use the fromFilePairs method but we will then need to get rid of the val in each of the tuples.
	// We can do this by using the map function as below
	Channel.fromFilePairs(["${params.mp_trimmed_reads_output_dir}/*{1P,2P}.fq.gz", "${params.pe_trimmed_reads_output_dir}/*{1P,2P}.fq.gz"]).map{it[1]}.set{ch_tadpole_input}
}


if (params.make_fastqc){
	// Now to the pre_trim fastqc for each of the subsampled fastq.gz files
	process fastqc_pre_trim{
		tag "${read}"
		conda "envs/stage_one.yaml"
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
		conda "envs/stage_one.yaml"
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
// The paper we were following was using quake to do the error correcting but
// this was a complete pain in the arse to get working and not available
// via conda so we will go with the bbtools offering and try to use tadpole
// TODO We will eventually want to check to see if the error correction has already been done
// and work from the containing directories if they have.
mp_dir_as_file = file(params.mp_corrected_reads_output_dir)
pe_dir_as_file = file(params.pe_corrected_reads_output_dir)
if ( mp_dir_as_file.exists() && pe_dir_as_file.exists() ){
	if ( mp_dir_as_file.list().size() == 12 && pe_dir_as_file.list().size() == 4){
			// Then all the files are already here and so we can work from this dir
			println("Successfully found error corrected seqs")
			params.correct_from_scratch = false
		}else{
			println("Could not find error corrected seqs")
			params.correct_from_scratch = true
		}
}else{
	println("Could not find error corrected seqs")
			params.correct_from_scratch = true
}
if (params.correct_from_scratch){
	// The output here will be a little difficult because we need to separate the mp and pe
	// reads into seperate channels. We may need to write a mapping function for this
	process tadpole{
		cache 'lenient'
		tag "${trimmed_read_one}"
		conda "envs/stage_one.yaml"
		// todpole should try to use all cores by default so I have set the tadpole_treads to be 24
		// so that nextflow will only try to run one at a time.
		// I haven't explicitly limited this in the running of tadpole.
		cpus params.tadpole_threads

		// These represent a considerable investment in time and so we should publish
		publishDir = [
				[path: params.mp_corrected_reads_output_dir, overwrite: 'true', pattern: "M_18*"],
				[path: params.pe_corrected_reads_output_dir, overwrite: 'true', pattern: "M_17*"]
			]

		input:
		tuple file(trimmed_read_one), file(trimmed_read_two) from ch_tadpole_input

		output:
		// We will put the pe files into one output for mapping against lutea
		// We will put the mp files into one output for going into nxtrim
		tuple file("*M_17*1P.cor.fastq.gz") optional true, file("*M_17*2P.cor.fastq.gz") optional true into ch_pe_map_against_lutea_input
		tuple file("*M_18*1P.cor.fastq.gz") optional true, file("*M_18*2P.cor.fastq.gz") optional true into ch_nxtrim_input

		script:
		output_one = trimmed_read_one.getName().replaceAll(".fq.gz", ".cor.fastq.gz")
		output_two = trimmed_read_two.getName().replaceAll(".fq.gz", ".cor.fastq.gz")
		
		"""
		tadpole.sh in=$trimmed_read_one in2=$trimmed_read_two out=$output_one out2=$output_two mode=correct
		"""
	}
}else{
	println("Populating channels from already corrected files")
	// Then the files have already been error corrected and we can create our channels directly from the directories
	// The pe files can go straight into the map against lutea but the mp files need to go through nxtrim first
	Channel.fromFilePairs("${params.pe_corrected_reads_output_dir}/*{1P,2P}.cor.fastq.gz").map{it[1]}.set{ch_pe_map_against_lutea_input}
	Channel.fromFilePairs("${params.mp_corrected_reads_output_dir}/*{1P,2P}.cor.fastq.gz").map{it[1]}.set{ch_nxtrim_input}

}


// The mate pair reads will go into NXtrim to characterise them, change orientation to innie and remove the adapters.
process nxtrim{
	cache 'lenient'
	tag "$mp_corrected_fastq_gz_one"
	conda "envs/stage_one.yaml"

	publishDir path: params.mp_nxtrim_output_dir

	input:
	tuple file(mp_corrected_fastq_gz_one) , file(mp_corrected_fastq_gz_two) from ch_nxtrim_input

	output:
	tuple file("*mp.fastq.gz"), file("*pe.fastq.gz"), file("*se.fastq.gz"), file("*unknown.fastq.gz") into ch_mp_map_against_lutea_input

	script:
	
	sample_name_out = mp_corrected_fastq_gz_one.getName().replaceAll(".fastq.gz", "")
	"""
	nxtrim -1 $mp_corrected_fastq_gz_one -2 $mp_corrected_fastq_gz_two -O $sample_name_out
	"""
}

// We can do a bbmap against lutea and a bwa map against lutea
// We have done a comparison of a small number of seq files to check for differences between bbmap and bwa.
// For the mp file M_18_1922_HUME-ST-35-45_AD002_Lane1.trimmed.1P.cor.mp.fastq.gz
// the number of proper pairs seem very comparable, as do the estimated insert sizes and standard deviations.
// In terms of the sequences that are contained by each there is also a large agreement (80%). 
// I.e. if there were 100 unique proper pairs found in total between bowth mapping types, 80 of them
// were found in both mapping type results.
// The results for the se reads in single read mode i.e. using interleaved=false with bbmap, were poor in agreement
// Only a very small proportion were matched using bwa wheras many more were matched using bbmap even with minratio set to 0.9
// When working with the pe files, we again saw very good agreement between the propoer pairs that were mapped approx 47% mapped
// And we saw that it was basically the same sequnces that were mapped. Again about 80%.
// The insert sizes and standard deviations were a bit higher with bbmap compared to bwa.
// bbmap mean was 353 and stdv 1174 vs for bwa that was mean 248 stdev 515.
// I think over all we should probably be working with bwa here rather than bbmap.
// Here are a few useful notes for running this part of the analysis:
// To run the bbmap mapping we were running:
// bbmap.sh ref=../../plutea_reference_genome/plut_v1.1.genome.fa in=M_18_1922_HUME-ST-35-45_AD002_Lane1.trimmed.1P.cor.pe.fastq.gz out=M_18_1922_HUME-ST-35-45_AD002_Lane1.trimmed.1P.cor.pe.bbmap.sam
// When we were doing it with the se, we also added interleaved=true and minratio=0.9
// For bwa mapping we did:
// bwa mem -p -t 26 ../../plutea_reference_genome/plut_v1.1.genome.fa M_18_1922_HUME-ST-35-45_AD002_Lane1.trimmed.1P.cor.pe.fastq.gz > M_18_1922_HUME-ST-35-45_AD002_Lane1.trimmed.1P.cor.pe.bwa.sam
// To convert sam to bam --> samtools view -S -b <sam file> > <bam file>
// To create a text file that was only the proper pairs:
// samtools view -f 2 M_18_1922_HUME-ST-35-45_AD002_Lane1.trimmed.1P.cor.pe.bwa.bam > proper_pairs.bwa.txt
// Then to get the agreement between sequences we ran both of the proper_pairs files as input to the compare_bams.py in bin.
// Then to get the insert mean and stdev we used an awk command:
//  awk -F '\t' '{if ($9 > 0) {tot+=$9; totsqr+= $9^2; cnt++;} } END{printf "mean insert size of %d; stdv of %d \n", (tot/cnt), (sqrt(totsqr/cnt-(tot/cnt)^2))}' proper_pairs.bbmap.txt
// Obiously we didn't do any of the paired testing when we were running the se files.

// For the output of the mapping we will want a couple of different threads
// For the M_17 reads (pe) we will want them going into the discovar assmbly.
// To do that we will want the unmapped reads in fastq.gz format.
// Then we'll want to re map these same reads to the discovar assembly along with the nxtrim unmapped reads
// critically we are going to want to ensure that the discovar assembly is created first.
// before we try to map the reads to it. 
// so it is important that that we don't simply have a channel that goes into the discovar assemply with just the pe
// reads and one that goes into the discovar mapping with both the mp and the pe reads. Rather we need to have one
// channel that contains the pe reads that will go into the discovar assembly making and a second channel containing
// jus tthe mp reads. Then as an ouput from the discovar creation we will re output the pe reads and use this channel
// joined with the mp channel from the lutea mapping joined together to map against the lutea. The other option
// would be to have an output from the discovar that was a file object for the actual disovarassembled .fa but 
// I don't know how easy it would be to pass this around and if we'd have to end up having multiple copies of it that would take
// up loads of space. So for the time being lets try the first approach.
// First thing to do is to map against lutea using bwa, discard the properpairs and then convertback to fastq.gz

process bwa_map_against_lutea{
	cache 'lenient'
	tag "$fastq_gz_to_map"
	conda "envs/stage_one.yaml"
	cpus params.bbmap_lutea_threads
	publishDir = [
		[path: params.mp_unmapped_to_lutea_reads_output_dir, overwrite: 'true', pattern: "M_18*"],
		[path: params.pe_unmapped_to_lutea_reads_output_dir, overwrite: 'true', pattern: "M_17*"]
	]

	// Here we will not work in tuples but rather work seq file by seq file
	// So we will need to mix the two channels that are coming in and then flatten them

	input:
	// The mp files will be coming in as individual files
	// but we will have the pe files coming in as tuples.
	file (fastq_gz_to_map) from ch_mp_map_against_lutea_input.flatten().mix(ch_pe_map_against_lutea_input)

	// One channel for the pe outputs to go into the discovar assembly
	// One channel for the mp outputs that will go into disovar mapping
	output:
	file "*M_17*bwa_lutea.unmapped.fastq.gz" optional true into ch_discovar_assembly_creation_input
	file "*M_18*bwa_lutea.unmapped.fastq.gz" optional true into ch_mp_discovar_mapping_input

	script:
	// Do mapping (fastq.gz --> sam)
	// Conver to bam (sam --> bam)
	// Discard non proper pairs or mapped reads (bam --> bam)
	// Convert bam to fastq (bam --> fastq)
	// compress fastq (fastq --> fastq.gz)
	// Clean up unused files

	// This conditional will check to see if we are working with a tuple (i.e. the pe)
	// or not (i.e. the mp)
	if (fastq_gz_to_map.size() > 1){
		// Then this is a tuple and so must be the pe reads
		read_one = fastq_gz_to_map[0]
		read_two = fastq_gz_to_map[1]
		out_sam_file_name = read_one.getName().replaceAll("fastq.gz", "bwa_lutea.sam")
	}else{
		read_one = fastq_gz_to_map
		out_sam_file_name = read_one.getName().replaceAll("fastq.gz", "bwa_lutea.sam")
	}
	out_bam_file_name = out_sam_file_name.replaceAll(".sam", ".bam")
	out_bam_unmapped_file_name = out_bam_file_name.replaceAll(".bam", ".unmapped.bam")
	out_fastq_unmapped_file_name = out_bam_unmapped_file_name.replaceAll(".bam", ".fastq")
	"""
	if [[ $read_one == *".se."* ]]; then
		bwa mem -t 24 ${params.lutea_ref_genome_path} $read_one > $out_sam_file_name
		samtools view -S -b $out_sam_file_name > $out_bam_file_name
		samtools view -b -f 4 $out_bam_file_name > $out_bam_unmapped_file_name
	elif [[ $read_one == *"M_17"* ]]; then
		bwa mem -t 24 ${params.lutea_ref_genome_path} $read_one $read_two > $out_sam_file_name
		samtools view -S -b $out_sam_file_name > $out_bam_file_name
		samtools view -b -F 2 $out_bam_file_name > $out_bam_unmapped_file_name
	else
		bwa mem -p -t 24 ${params.lutea_ref_genome_path} $read_one > $out_sam_file_name
		samtools view -S -b $out_sam_file_name > $out_bam_file_name
		samtools view -b -F 2 $out_bam_file_name > $out_bam_unmapped_file_name
	fi
	samtools fastq $out_bam_unmapped_file_name > $out_fastq_unmapped_file_name
	gzip $out_fastq_unmapped_file_name
	rm $out_sam_file_name $out_bam_file_name $out_bam_unmapped_file_name
	"""

}

// ch_discovar_assembly_creation_input.collect().view()

// // Here we have the unmapped fastq.gz files for both the mp and pe reads.
// // Here we will do a discovar assembly using on the pe reads. We will then use this assembly to map
// // reads to to get the ideas of insert lengths that we will eventually use for the LG_allpaths assembly

// I managed to invalidate the cache by accidentally modifying one of the output files.
// As such I will check to see if the discovar assembly file has already been created and skip this if it has.
// Check to see if the a.lines.fasta file already exists. If it does, then have a basic process that
// collects the inputs from ch_discovar_assembly_creation_input and puts them into ch_pe_discovar_mapping_input
// NB remember that index files must be created from the a.lines.fasta file for bwa mem to use.
// run bwa index a.lines.fasta and then point bwa mem to the a.lines.fasta file. It will find the index files.
discovar_denovo_fasta_file_path = file("${discovar_publish_path}/output/a.final/a.lines.fasta")
if ( discovar_denovo_fasta_file_path.exists() ){
	ch_discovar_assembly_creation_input.collect().set{ch_pe_discovar_mapping_input}
}else{
	process do_discovar_assembly{
	cache 'lenient'
	tag "discovar_assembly"
	conda "envs/stage_one.yaml"
	publishDir path: discovar_publish_path

	input:
	tuple file(unmapped_pe_fastq_gz_one), file(unmapped_pe_fastq_gz_two) from ch_discovar_assembly_creation_input.collect()

	output:
	// We will keep the output a.lines.fasta to do our insert map against
	tuple file("**/stats"), file("**/a.lines.fasta") into ch_discovar_out
	tuple file(unmapped_pe_fastq_gz_one), file(unmapped_pe_fastq_gz_two) into ch_pe_discovar_mapping_input

	script:
		"""
		mkdir output
		DiscovarDeNovo READS=${unmapped_pe_fastq_gz_one},${unmapped_pe_fastq_gz_two} OUT_DIR=\${PWD}/output/ MAX_MEM_GB=450
		"""
	}
}

// Once we have the discovar assembly we can map the pe reads and mp reads against it to get
// our insert lengths.
// We will not map the se reads as these do not have inserts
process bwa_map_against_discovar{
	cache 'lenient'
	tag "$fastq_gz_to_map"
	conda "envs/stage_one.yaml"
	cpus params.bbmap_lutea_threads
	publishDir = [
		[path: mp_insert_size_info_output_dir, overwrite: 'true', pattern: "M_18*"],
		[path: pe_insert_size_info_output_dir, overwrite: 'true', pattern: "M_17*"]
	]

	input:
	// All of the files coming in will now be single files that are interleaved pairs except for the 
	// se files that should be treated as single pairs when mapping.
	file fastq_gz_to_map from ch_mp_discovar_mapping_input.mix(ch_pe_discovar_mapping_input.flatten())


	// We will use all of the size info together and we may as well store it
	// with the file that it is related to. So, let's output to a single channel
	// That can be put into the prep process for the allpathsLG assembly
	// NB there is little point in doing a mapping for the se reads as these don't have insert lengths
	// Instead we will simply output a dud size_info file
	output:
	tuple file(fastq_gz_to_map), file("*bwa_discovar.mapped.size_info.txt") into ch_all_paths_input
	
	script:
	// Do mapping (fastq.gz --> sam)
	// Conver to bam (sam --> bam)
	// Put proper pairs into .txt for awk (bam --> txt)
	// Use awk to get the average insert size and the stdev
	// Clean up unused files
	discovar_assembly_fasta_path = discovar_denovo_fasta_file_path
	out_sam_file_name = fastq_gz_to_map.getName().replaceAll("fastq.gz", "bwa_discovar.sam")
	out_bam_file_name = out_sam_file_name.replaceAll(".sam", ".bam")
	out_mapped_txt_file_name = out_bam_file_name.replaceAll(".bam", ".mapped.txt")
	size_info_output_file_name = out_mapped_txt_file_name.replaceAll(".txt", ".size_info.txt")
	"""
	if [[ $fastq_gz_to_map == *".se."* ]]; then
		echo No size info for se reads > $size_info_output_file_name
	else
		bwa mem -p -t 24 ${discovar_assembly_fasta_path} $fastq_gz_to_map > $out_sam_file_name
		samtools view -S -b $out_sam_file_name > $out_bam_file_name
		samtools view -f 2 $out_bam_file_name > $out_mapped_txt_file_name
		awk -F '\\t' '{if (\$9 > 0) {tot+=\$9; totsqr+= \$9^2; cnt++;} } END{printf "mean insert size of %d; stdv of %d \\n", (tot/cnt), (sqrt(totsqr/cnt-(tot/cnt)^2))}' $out_mapped_txt_file_name > $size_info_output_file_name
		rm $out_sam_file_name $out_bam_file_name $out_mapped_txt_file_name
	fi
	"""
}