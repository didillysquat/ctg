configfile: "config.yaml"

LANES = ["Lane1", "Lane2"]
DIRECTIONS = ["R1", "R2"]
MP_SAMPLES = ["M_18_1922_HUME-ST-35-45_AD002", "M_18_1923_HUME-ST-5-7_AD007", "M_18_1924_HUME-ST-8-11_AD019"]


# rule all:
#     input:
#         "mate_pair_libs/fastqc/fastqc_complete.txt",
#         "paired_end_reads/fastqc/fastqc_complete.txt"

# Create fastqc files for each of the sequencing files before we start trimming and correction
rule make_fastqc_mp_pre_trim:
    input:
        expand("mate_pair_libs/{mp_sample}_{lane}_{direction}_001.fastq.gz", mp_sample=MP_SAMPLES, lane=LANES, direction=DIRECTIONS)
    threads:12
    output:
        expand("mate_pair_libs/fastqc/{mp_sample}_{lane}_{direction}_001.fastqc.zip", mp_sample=MP_SAMPLES, lane=LANES, direction=DIRECTIONS),
        "mate_pair_libs/fastqc/fastqc_complete.txt"
    shell:
        "parallel -j {threads} fastqc -o mate_pair_libs/fastqc ::: {input}; "
        "touch mate_pair_libs/fastqc/fastqc_complete.txt"


rule make_fastqc_pe_pre_trim:
    input:
        expand("paired_end_reads/M_17_426_1_AD08_{lane}_{direction}_001.fastq.gz", lane=LANES, direction=DIRECTIONS)
    threads:4
    output:
        expand("paired_end_reads/fastqc/M_17_426_1_AD08_{lane}_{direction}_001.fastq.gz", lane=LANES, direction=DIRECTIONS),
        "paired_end_reads/fastqc/fastqc_complete.txt"
    shell:
        "parallel -j {threads} fastqc -o paired_end_reads/fastqc ::: {input}; "
        "touch paired_end_reads/fastqc/fastqc_complete.txt"

# Use trimmomatic to remove 3' adapters and and to remove first few bp from beginning
# snakemake --use-conda --cores 24 mate_pair_libs/trimmed/M_18_1922_HUME-ST-35-45_AD002_Lane1.trimmed_1P.fq.gz
# mate_pair_libs/trimmed/M_18_1922_HUME-ST-35-45_AD002_Lane2.trimmed_1P.fq.gz
# mate_pair_libs/trimmed/M_18_1923_HUME-ST-5-7_AD007_Lane1.trimmed_1P.fq.gz
# mate_pair_libs/trimmed/M_18_1923_HUME-ST-5-7_AD007_Lane2.trimmed_1P.fq.gz
# mate_pair_libs/trimmed/M_18_1924_HUME-ST-8-11_AD019_Lane1.trimmed_1P.fq.gz
# mate_pair_libs/trimmed/M_18_1924_HUME-ST-8-11_AD019_Lane2.trimmed_1P.fq.gz
rule trim_mp:
    input:
        "mate_pair_libs/{mp_sample}_R1_001.fastq.gz"
    output:
        "mate_pair_libs/trimmed/{mp_sample}.trimmed_1P.fq.gz",
        "mate_pair_libs/trimmed/{mp_sample}.trimmed_2P.fq.gz"
    threads: 6
    conda:
        "envs/trim_correct.yaml"
    shell:
        "trimmomatic PE -threads {threads} -basein {input} -baseout mate_pair_libs/trimmed/{wildcards.mp_sample}.trimmed.fq.gz "
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11"


rule trim_pe:
    input:
        "paired_end_reads/{pe_sample}_R1_001.fastq.gz"
    output:
        "paired_end_reads/trimmed/{pe_sample}.trimmed_1P.fq.gz",
        "paired_end_reads/trimmed/{pe_sample}.trimmed_2P.fq.gz"
    threads: 12
    conda:
        "envs/trim_correct.yaml"
    shell:
        "trimmomatic PE -threads {threads} -basein {input} -baseout paired_end_reads/trimmed/{wildcards.pe_sample}.trimmed.fq.gz "
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11"


DIRECTIONS_TRIM = ["1P", "2P"]
# Create fastqc files after trimming
rule make_fastqc_mp_post_trim:
    input:
        expand("mate_pair_libs/trimmed/{mp_sample}_{lane}.trimmed_{direction}.fq.gz", mp_sample=MP_SAMPLES, lane=LANES, direction=DIRECTIONS_TRIM)
    threads:6
    output:
        expand("mate_pair_libs/fastqc_post_trim/{mp_sample}_{lane}.trimmed_{direction}.fastqc.zip", mp_sample=MP_SAMPLES, lane=LANES, direction=DIRECTIONS_TRIM),
        expand("mate_pair_libs/fastqc_post_trim/{mp_sample}_{lane}.trimmed_{direction}.fastqc.html", mp_sample=MP_SAMPLES, lane=LANES, direction=DIRECTIONS_TRIM),
        "mate_pair_libs/fastqc_post_trim/fastqc_complete.txt"
    shell:
        "parallel -j {threads} fastqc -o mate_pair_libs/fastqc_post_trim ::: {input}; "
        "touch mate_pair_libs/fastqc_post_trim/fastqc_complete.txt"


rule make_fastqc_pe_post_trim:
    input:
        expand("paired_end_reads/trimmed/M_17_426_1_AD08_{lane}.trimmed_{direction}.fq.gz", lane=LANES, direction=DIRECTIONS_TRIM)
    threads:4
    output:
        expand("paired_end_reads/fastqc_post_trim/M_17_426_1_AD08_{lane}.trimmed_{direction}.fastqc.html", lane=LANES, direction=DIRECTIONS_TRIM),
        expand("paired_end_reads/fastqc_post_trim/M_17_426_1_AD08_{lane}.trimmed_{direction}.fastqc.zip", lane=LANES, direction=DIRECTIONS_TRIM),
        "paired_end_reads/fastqc/fastqc_complete.txt"
    shell:
        "parallel -j {threads} fastqc -o paired_end_reads/fastqc_post_trim ::: {input}; "
        "touch paired_end_reads/fastqc_post_trim/fastqc_complete.txt"

# Use rcorrect to do kmer based error correction of the sequences
# snakemake -np --cores 24 err_corrected/b_minutum/SRR17933{20..23}.trimmed_{1,2}P.cor.fq.gz err_corrected/b_psygmophilum/SRR17933{24..27}.trimmed_{1,2}P.cor.fq.gz
rule err_correct:
	input:
		fwd = "trimmed/{species}/{sra}.trimmed_1P.fq.gz",
		rev = "trimmed/{species}/{sra}.trimmed_2P.fq.gz"
	output:
		"err_corrected/{species}/{sra}.trimmed_1P.cor.fq.gz",
		"err_corrected/{species}/{sra}.trimmed_2P.cor.fq.gz"
	conda:
		"envs/strain_deg.yaml"
	threads: 6
	shell:
		"run_rcorrector.pl -1 {input.fwd} -2 {input.rev} -od err_corrected/{wildcards.species}/ -t {threads}"