#!/usr/bin/env nextflow

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
		}else{
			params.sub_sample_from_scratch = true
		}
	}else{
		params.sub_sample_from_scratch = true
	}
}else{
	params.mp_wkd = "${params.basedir}/mate_pair_libs"
	params.pe_wkd = "${params.basedir}/paired_end_reads"
}

// We will enable subsampling for development purposes.
// We will subsample all of the sequencing reads down to 1 000 000 reads. 

if (params.sub_sample){
				

	// The directory for the sub_sampled raw reads to be output
	params.mp_dir_sub_raw_reads = "${params.mp_wkd}/raw_reads/${params.sub_sample_depth}" 	
	params.pe_dir_sub_raw_reads = "${params.pe_wkd}/raw_reads/${params.sub_sample_depth}"
	println("params.sub_sample_from_scratch is ${params.sub_sample_from_scratch}")
	if (params.sub_sample_from_scratch){
		// The full size raw files that will be used as input to the seqtk
		// We will put all of the sequencing files into a single channel
		Channel.fromFilePairs(["${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz", "${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz"]).set{ch_sub_sample_seqtk_input}
		
		 	
		
		// Now do the subsampling using seqtk
		// Once seqtk is complete then will need to do the gz compression
		// Then will need to put this gz compressed file into the ch_fastqc_pr channel for the fastqc
		// and the ch_trim_input for the trimming
		process sub_sample_seqtk{
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
		// Then we generate the required chanels directly from the already existing subsample directory
		// The directory for the sub_sampled raw reads to be output
		println("WE FOT HERE")
		Channel.fromPath(["${params.pe_dir_sub_raw_reads}/*.fastq.gz", "${params.mp_dir_sub_raw_reads}/*.fastq.gz"]).set{ch_fastqc_pre_trim_input}
		Channel.fromFilePairs(["${params.pe_dir_sub_raw_reads}/*{R1_sub,R2_sub}.fastq.gz", "${params.mp_dir_sub_raw_reads}/*{R1_sub,R2_sub}.fastq.gz"]).into{ch_trim_reads_input; ch_test}
	}
	
}else{
	Channel.fromPath(["${params.pe_dir}/raw_reads/*.fastq.gz", "${params.mp_dir}/raw_reads/*.fastq.gz"]).set{ch_fastqc_pre_trim_input}
	Channel.fromFilePairs(["${params.pe_dir}/raw_reads/*{R1,R2}.fastq.gz", "${params.mp_dir}/raw_reads/*{R1,R2}.fastq.gz"]).set{ch_trim_reads_input}
}


// Now trim the reads using trimmomatic
process trim_reads{
	tag "${reads[0].getName()}"
	conda "envs/stage_one.yaml"
	
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
	paired_out_one = reads[0].getName().replaceAll("_sub", "").replaceAll('_R1.fastq.gz', '.trimmed.1P.fq.gz')
	unpaired_out_one = reads[0].getName().replaceAll("_sub", "").replaceAll('_R1.fastq.gz', '.trimmed.1U.fq.gz')
	paired_out_two = reads[1].getName().replaceAll("_sub", "").replaceAll('_R2.fastq.gz', '.trimmed.2P.fq.gz')
	unpaired_out_two = reads[1].getName().replaceAll("_sub", "").replaceAll('_R2.fastq.gz', '.trimmed.2U.fq.gz')
	"""
	trimmomatic PE -threads ${params.trimmomatic_threads} ${reads[0]} ${reads[0]} \\
		$paired_out_one $unpaired_out_one $paired_out_two $unpaired_out_two\\
		ILLUMINACLIP:${params.tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
		LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11
	"""
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
// The output here will be a little difficult because we need to separate the mp and pe
// reads into seperate channels. We may need to write a mapping function for this
process rcorrector{
	tag "${trimmed_read_one}"
	conda "envs/stage_one.yaml"
	cpus params.rcorrector_threads

	input:
	tuple file(trimmed_read_one), file(trimmed_read_two) from ch_rcorrector_input

	output:
	tuple file("*1P.cor.fq.gz"), file("*2P.cor.fq.gz") into ch_bbmerge_input, ch_nxtrim_input, ch_map_pe_against_lutea_unmerged_input

	script:
	"""
	run_rcorrector.pl -1 $trimmed_read_one -2 $trimmed_read_two -od . -t ${task.cpus}
	"""
}

// ch_mp_out.toList().view()


// // After error correction we will want to split the channels up by paired end and mate pair
// We will map diferently according to whether we are doing paired end or mate pair




// TODO it actually makes more sense to do the mapping after the bbmerge and nxtrim treatement
// SO that we can get an idea of insert length.

// // The paired end reads will go into BBMerge, to merge as many as possible
// process bbMerge{
// 	tag "$error_corrected_fastq_gz_one"
// 	conda "envs/stage_one.yaml"
// 	// I have set the CPUs here but will not limit the threads in the actual
// 	// command. This should be a good balance between still making use of as many threads
// 	// as possible but not having innefficient load distribution across the threads
// 	// i.e. spending most of the cpu time in the kernel.
// 	cpus params.bbmerge_threads
// 	input:
// 	// Importantly we only want the paried end reads here so the M_17 files
// 	tuple file(error_corrected_fastq_gz_one), file(error_corrected_fastq_gz_two) from ch_bbmerge_input.toList().flatMap{
// 		// Create the new list to return
// 		List output_list = new ArrayList();
// 		// for each tuple in the list
// 		for (i=0; i<(it.size()); i++){
// 			if (it[i][0].getName().contains("M_17")){
// 				output_list.add(it[i])
// 			}
// 		}
// 		return output_list
//     }

// 	output:
// 	tuple file("*.merged.fastq.gz"), file("*.unmerged.fastq.gz") into ch_pe_bbmap_filter_lutea_input, ch_pe_bbmap_gc_investigation_input, ch_pe_bbmap_gc_investigation_overall_input
// 	file "*ihist.txt" into ch_bbmerge_ihist_out

// 	script:
// 	out_merged_name = error_corrected_fastq_gz_one.getName().replaceAll(".fq.gz", ".merged.fastq.gz")
// 	unmerged_name = error_corrected_fastq_gz_one.getName().replaceAll(".fq.gz", ".unmerged.fastq.gz")
// 	ihist_name = error_corrected_fastq_gz_one.getName().replaceAll(".fq.gz", ".ihist.txt")
// 	"""
// 	bbmerge-auto.sh in1=$error_corrected_fastq_gz_one in2=$error_corrected_fastq_gz_two out=$out_merged_name outu=$unmerged_name ihist=$ihist_name ecct extend2=20 iterations=5
// 	"""
// }

// // The mate pair reads will go into NXtrim to characterise them

// process nxtrim{
// 	tag "$error_corrected_fastq_gz_one"
// 	conda "envs/stage_one.yaml"

// 	input:
// 	// Importantly we only want the mate pair reads here so the M_18 files
// 	tuple file(error_corrected_fastq_gz_one), file(error_corrected_fastq_gz_two) from ch_nxtrim_input.toList().flatMap{
// 		// Create the new list to return
// 		List output_list = new ArrayList();
// 		// for each tuple in the list
// 		for (i=0; i<(it.size()); i++){
// 			if (it[i][0].getName().contains("M_18")){
// 				output_list.add(it[i])
// 			}
// 		}
// 		return output_list
//     }	

// 	output:
// 	tuple file("*mp.fastq.gz"), file("*pe.fastq.gz"), file("*se.fastq.gz"), file("*unknown.fastq.gz") into ch_mp_bbmap_filter_lutea_input, ch_mp_bbmap_gc_investigation_input, ch_mp_bbmap_gc_investigation_overall_input

// 	script:
// 	sample_name_out = error_corrected_fastq_gz_one.getName().replaceAll(".fq.gz", "")
// 	// We need to have paired reads here rather than the single interleaved fastq that we get
// 	// from the bbmap. We will use the deinterleave_fastq.sh in bin to unpair the reads
// 	"""
// 	nxtrim -1 $error_corrected_fastq_gz_one -2 $error_corrected_fastq_gz_two -O $sample_name_out
// 	"""
// }

if (params.do_bbmap_lutea_gc_investigations){
	// bbmap is really useful. You can set the minratio score which is the ratio of th actual best score
	// over the best possible score. It basically a match quality cutoff.
	// I want to try something quite ambitious which is to compare the percentage of mapping
	// to the GC content of the mapped and unmapped reads.
	// In theory this will tell us whether the mappings are really host or whether they are
	// zooxs reads being mapped to the host genome.
	// To do this, I will use the each method of next flow to run instances of bbmap for every
	// library and for each library I will run 5% quality thresholds sets.
	// We will then be able to plot up the matches vs the GC content change.
	// To minimise space we can delete the output files once they have been made.
	// The other infomation we need gets output to stdout so we will collect that too.
	// Once the bbmap is complete we will run stats.sh to generate the gc content of the mapped and unmapped
	// We will then need to collect all of this for plotting with Python.
	// This can be a figure in the paper.

	// The command for getting average length and stdev: awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f\n",n,m,sqrt(sq/n-m*m));}' *.fastq
	// NB we were having some problems with the new line character in the above awk command so we removed it.
	// The results are very interesting. They show that the pe reads should probably have a higher
	// minratio than the mp reads. For the pe a ratio of 0.85 looks like it would be good
	// For the mp much lower i.e. 0.6 looks like it would be appropriate.
	quality_values = [1, 0.95, 0.90, 0.85, 0.8, 0.75, 0.7, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30]

	//TODO the gc values that are being output are really quite crap.
	// I will calculate them by hand using some more awesome awk and then we can replot
	process bbmap_gc_investigation{
		tag "${fastq_gz_to_map}_$qual_val"
		conda "envs/stage_one.yaml"
		// I have set the CPUs here but will not limit the threads in the actual
		// command. This should be a good balance between still making use of as many threads
		// as possible but not having innefficient load distribution across the threads
		// i.e. spending most of the cpu time in the kernel.
		cpus params.bbmap_gc_investigation_threads
		publishDir path: "gc_lutea_maping"

		input:
		file fastq_gz_to_map from ch_pe_bbmap_gc_investigation_input.toList().mix(ch_mp_bbmap_gc_investigation_input.toList()).flatten()
		each qual_val from quality_values

		output:
		tuple file("*.log.txt"), file("*.gchist.txt"), file("*.unmapped.awk_gc.txt"), file("*.mapped.awk_gc.txt") into ch_mapped_unmapped_output

		script:
		// This is very useful: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
		seq_sample_basename = fastq_gz_to_map.getName().replaceAll(".fastq.gz", "")
		"""
		# First run bbmap.sh to produce the mapped and unmapped fastqs
		# By defualt, a pair that have one read mapped and one unmapped will be put into the mapped file (i.e. it is conservative). We will not change this default.
		# Because we have already run nxtrim on the mp libraries it has already changed RF orientation to FR so we don't have to worry about the rcs=f
		# Because we have no merged some of the pe reads, this makes some of the reads longer than the maximum size of 600 that bbmap.sh allows.
		# We can use the option maxlen=600 to break up the longer sequences into chuncks. It seems that the chunks are put pack together in the out put as the
		# average length of the mapped and unmapped are almost identical to the input sequences. As are the standard deviations.
		# We should not that for the single (i.e. non paired) reads (e.g. the merged reads) the log file does not contain insert data length so we will have to work this out on the command line.
		# I will use a cool awk command to calculate this and put it into len_info.txt files. This will give us average and stdev of the seq lengths. (Acutally, we will only need this info when we do the actual mapping)
		
		bbmap.sh in=$fastq_gz_to_map ref=${params.lutea_ref_genome_path} outm=${seq_sample_basename}.mapped.fastq outu=${seq_sample_basename}.unmapped.fastq nodisk printunmappedcount minratio=$qual_val k=14 gchist=${seq_sample_basename}.gchist.txt maxsites=1 gcbins=1000 maxlen=600
		# then remove the large datafiles
		awk '(NR%4==2) {N1+=length(\$0);gsub(/[AT]/,"");N2+=length(\$0);}END{print N2/N1;}' ${seq_sample_basename}.unmapped.fastq > ${seq_sample_basename}.${qual_val}.unmapped.awk_gc.txt
		awk '(NR%4==2) {N1+=length(\$0);gsub(/[AT]/,"");N2+=length(\$0);}END{print N2/N1;}' ${seq_sample_basename}.mapped.fastq > ${seq_sample_basename}.${qual_val}.mapped.awk_gc.txt
		rm ${seq_sample_basename}.unmapped.fastq
		rm ${seq_sample_basename}.mapped.fastq
		# The output that we need is not written to stdout but rather to the .command.log file output by nextflow
		# We will attempt to rename this and extract it to the outputs
		cp .command.log ${seq_sample_basename}.${qual_val}.log.txt
		"""
	}

	// The above process still does not do an awk calculation of the overall gc value of the inital reads that go into bbmap. We will only need to do this once per
	// library. As I have already done the above process I will calculate this in a seperate process
	process get_overall_gc_of_libs{
		tag "${fastq_gz_to_map}"
		conda "envs/stage_one.yaml"
		// I have set the CPUs here but will not limit the threads in the actual
		// command. This should be a good balance between still making use of as many threads
		// as possible but not having innefficient load distribution across the threads
		// i.e. spending most of the cpu time in the kernel.
		cpus params.bbmap_gc_investigation_threads
		publishDir path: "gc_lutea_mapping_investigation"

		input:
		file fastq_gz_to_map from ch_pe_bbmap_gc_investigation_overall_input.toList().mix(ch_mp_bbmap_gc_investigation_overall_input.toList()).flatten()
		
		output:
		file "*.all.awk_gc.txt" into ch_get_overall_gc_of_libs_output

		script:
		// This is very useful: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
		seq_sample_basename = fastq_gz_to_map.getName().replaceAll(".fastq.gz", "")
		"""
		# First run bbmap.sh to produce the mapped and unmapped fastqs
		# By defualt, a pair that have one read mapped and one unmapped will be put into the mapped file (i.e. it is conservative). We will not change this default.
		# Because we have already run nxtrim on the mp libraries it has already changed RF orientation to FR so we don't have to worry about the rcs=f
		# Because we have no merged some of the pe reads, this makes some of the reads longer than the maximum size of 600 that bbmap.sh allows.
		# We can use the option maxlen=600 to break up the longer sequences into chuncks. It seems that the chunks are put pack together in the out put as the
		# average length of the mapped and unmapped are almost identical to the input sequences. As are the standard deviations.
		# We should not that for the single (i.e. non paired) reads (e.g. the merged reads) the log file does not contain insert data length so we will have to work this out on the command line.
		# I will use a cool awk command to calculate this and put it into len_info.txt files. This will give us average and stdev of the seq lengths. (Acutally, we will only need this info when we do the actual mapping)
		
		gzip -df $fastq_gz_to_map
		awk '(NR%4==2) {N1+=length(\$0);gsub(/[AT]/,"");N2+=length(\$0);}END{print N2/N1;}' ${seq_sample_basename}.fastq > ${seq_sample_basename}.all.awk_gc.txt
		rm ${seq_sample_basename}.fastq
		"""
	}

	// In theory we should just be able to flatten this.
	// process bbmap_gc_supp_figure{
	// 	tag "${gc_content_hist_mapped.getName().replaceAll(".gc_content_hist_mapped.txt", "")}"
	// 	conda "envs/standard_python.yaml"
	// 	publishDir path: "figures"

	// 	input:
	// 	// tuple file(log.txt), file(gchist.txt), file(gc_content_hist_unmapped.txt), file(gc_content_hist_mapped.txt) from ch_mapped_unmapped_output.collect()
	// 	file bbmap_output_files from ch_mapped_unmapped_output.collect()

	// 	output:
	// 	tuple file("gc_to_mapping_pct.svg"), file("gc_to_mapping_pct.svg") into bbmap_gc_supp_figure_output

	// 	"""
	// 	python3 ${params.bin_dir}/summarise_gc_content_mapping_figure.py
	// 	"""
	// }
}


// ch_pe_bbmap_gc_investigation_input.toList().mix(ch_mp_bbmap_gc_investigation_input.toList()).flatten().view()

// // mapp the error corrected reads to the lutea genome using the predetermined minratio scores.
// // These were determined from the bbmap_gc_investigation process above
// // Above we were running the gc mapping investigation process once for every type and every library.
// // From, this we have predetermined minratios for each of the types of library to perform the mapping with.
// // Here we will also run once per type and library as we don't need to keep the sets togheter.
// // I.e. it doesn't matter if we don't have the pe, mp, se and unknowns in tuple,
// // we can just run them all seperately and derive which library they are from in the next
// // python script for making the data inputs to the AllPathsLG assembler.
// process bbmap_filter_lutea{
// 	tag "${fastq_gz_to_map}"
// 	conda "envs/stage_one.yaml"
// 	// I have set the CPUs here but will not limit the threads in the actual
// 	// command. This should be a good balance between still making use of as many threads
// 	// as possible but not having innefficient load distribution across the threads
// 	// i.e. spending most of the cpu time in the kernel.
// 	cpus params.bbmap_lutea_threads
// 	publishDir path: "gc_lutea_mapping"
// 	input:
// 	file fastq_gz_to_map from ch_pe_bbmap_filter_lutea_input.toList().mix(ch_mp_bbmap_filter_lutea_input.toList()).flatten()

// 	output:
// 	//From this output we will need the unmapped seq file so that we can get the name from that
// 	// We will also want to have either:
// 	// a) the mean length of the fragment and the standard deviation of this (for the non-matepair reads), or
// 	// b) the mean insert length and the standard deviation of this.
// 	// We will get the insert length and standard deviation hopefully from the log file output
// 	// Acutally only the mean length was available from the log file.
// 	// So rather we will output the ihist file for the matepair mapping.
// 	// This will allow us to calcuate mean and stdevs.
// 	// We will generate the men length and stdev info using the awk command shown below and be sure to collect
// 	// this output.
// 	// We will not output the mapped reads.
// 	// I think it will actually be worth having a look at the insert lengths for all of the paired sequences so I will add this in below.
// 	tuple file("*.unmapped.fastq"), file("*.unmapped.len_info.txt"), file("*.ihist.txt") optional true into ch_bbmap_filter_lutea_output

// 	script:
// 	// This is very useful: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
// 	seq_sample_basename = fastq_gz_to_map.getName().replaceAll(".fastq.gz", "")
// 	"""
	
// 	if [[ $seq_sample_basename == *".merged"* ]]; then
// 	# Then we need to run with the maxlen=600 because these are longer due to being merged
//   		bbmap.sh in=$fastq_gz_to_map ref=${params.lutea_ref_genome_path} outu=${seq_sample_basename}.unmapped.fastq nodisk printunmappedcount minratio=0.6 k=14 gchist=${seq_sample_basename}.gchist.txt maxsites=1 gcbins=1000 maxlen=600
// 	elif [[ $seq_sample_basename == *".mp"* ]] || [[ $seq_sample_basename == *"unknown"* ]]; then
// 		# Then these are matepair and we need to know the mean and stdev of the insert size.
// 		# We get these out of ihist.
// 		bbmap.sh in=$fastq_gz_to_map ref=${params.lutea_ref_genome_path} outu=${seq_sample_basename}.unmapped.fastq nodisk printunmappedcount minratio=0.6 k=14 gchist=${seq_sample_basename}.gchist.txt maxsites=1 gcbins=1000 ihist=${seq_sample_basename}.mapped.ihist.txt
// 	elif [[ $seq_sample_basename == *".se"* ]]; then
// 		# Then this is the single end and we will not require the ihist
// 		bbmap.sh in=$fastq_gz_to_map ref=${params.lutea_ref_genome_path} outu=${seq_sample_basename}.unmapped.fastq nodisk printunmappedcount minratio=0.6 k=14 gchist=${seq_sample_basename}.gchist.txt maxsites=1 gcbins=1000
// 	else
// 		# This for all other sequencing types that will be paired. We will require the ihist file for all of these.
// 		bbmap.sh in=$fastq_gz_to_map ref=${params.lutea_ref_genome_path} outu=${seq_sample_basename}.unmapped.fastq nodisk printunmappedcount minratio=0.6 k=14 gchist=${seq_sample_basename}.gchist.txt maxsites=1 gcbins=1000 ihist=${seq_sample_basename}.mapped.ihist.txt
// 	fi
// 	awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length(\$0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f",n,m,sqrt(sq/n-m*m));}' ${seq_sample_basename}.unmapped.fastq > ${seq_sample_basename}.unmapped.len_info.txt
// 	"""
// }

// TODO we have written a .py script for making the input data files for all pathslg and 
// for making histogras of the insertlengths that uses the outputs from bbmap_filter_lutea.
// We just need to make a process that uses it.

// TODO I think its probably worth putting in a discovar assembly using only the unmapped short reads and then doing the bbmap to this to get the insert lengths
// We can then compare these to the insert length estimates using the lutea genome.

// We will do an assembly using only the paired end sequences so that we can get a better idea of insert
// length distributions. We will then be able to compare the the distributions to those that we got from the lutea mapping
// Discovar takes in 250 length paired reads. As such, for the input to this we will need to work from before the bbmerge of the pe sequeneces
// We will then need to run these against the lutea genome and then we can use this paired unmapped fastq as input to discovar

process map_pe_against_lutea_unmerged{
	tag "$error_corrected_fastq_gz_one"
	conda "envs/stage_one.yaml"
	publishDir path: "gc_lutea_mapping_unmapped_bam"
	input:
	tuple file(error_corrected_fastq_gz_one), file(error_corrected_fastq_gz_two) from ch_map_pe_against_lutea_unmerged_input.toList().flatMap{
			// Create the new list to return
			List output_list = new ArrayList();
			// for each tuple in the list
			for (i=0; i<(it.size()); i++){
				if (it[i][0].getName().contains("M_17")){
					output_list.add(it[i])
				}
			}
			return output_list
		}
	
	output:
	file "*unmapped.bam" into ch_discovar_assembly_input
	
	script:
	seq_sample_basename = error_corrected_fastq_gz_one.getName().replaceAll(".fq.gz", "")
	"""
	bbmap.sh in1=$error_corrected_fastq_gz_one in2=$error_corrected_fastq_gz_two ref=${params.lutea_ref_genome_path} outu=${seq_sample_basename}.unmapped.bam nodisk printunmappedcount minratio=0.6 k=14 maxsites=1
	"""
}

// process fastq_to_bam{
// 	tag "$unmapped_fastq_gz"
// 	conda "envs/stage_one.yaml"
// 	publishDir path: "gc_lutea_mapping_unmapped_bam"
// 	input:
// 	file unmapped_fastq_gz from ch_fastq_to_bam_input

// 	output:
// 	file "*.unaligned.bam" into ch_discovar_input

// 	script:
// 	sam_output = unmapped_fastq_gz.getName().replaceAll(".fastq.gz", ".unaligned.bam")
// 	sample_name = unmapped_fastq_gz.getName().replaceAll(".fastq.gz", "")
// 	"""
// 	picard FastqToSam F1=$unmapped_fastq_gz O=$sam_output SM=$sample_name
// 	"""
// }

if (params.sub_sample){
	discovar_publish_path = "${params.pe_wkd}/raw_reads/${params.sub_sample_depth}_discovar_assembly"
}else{
	discovar_publish_path = "${params.pe_wkd}/discovardenovo_assembly"
}

process discovar_assembly{
tag "Discovar_assembly"
conda "envs/stage_one.yaml"
publishDir path: discovar_publish_path
input: 
tuple file(unaligned_paired_bam_one), file(unaligned_paired_bam_two) from ch_discovar_assembly_input.toList()

output:
// We will keep the output a.lines.fasta to do our insert map against
tuple file("**/stats"), file("**/a.lines.fasta") into ch_discovar_out

script:
	"""
	mkdir output
	
	DiscovarDeNovo READS=${unaligned_paired_bam_one},${unaligned_paired_bam_two} OUT_DIR=\${PWD}/output/
	"""
}

// 
// ch_bbmerge_out.toList().view()
// ch_nxtrim_output.toList().view()


// So that we don't lose the ability to use the resume option in nextflow, we will start another conda environment here
// We'll call it stage_two.yaml.
// Eventually we will want to merge all of the environments together if possible.

// The liu et al 2018 paper that we were following used clc workbench to do an initial assembly to get accurate estimate
// of the insert lengths. However, this is commercial software and I'd rather stay clear.
// For their main assembly, they use both SPADES and ALLPATHS-LG. We could use one of these to make an initial assembly to check insert sizes.

// We also need to take into account the fact that we still have coral reads in our data. Let's dump those out before we do our first assembly.
// It looks like bbmap will be a good way to map the various reads to our p lutea genome.

