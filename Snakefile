configfile: "config.yaml"

LANES = ["Lane1", "Lane2"]
DIRECTIONS = ["R1", "R2"]
MP_SAMPLES = ["M_18_1922_HUME-ST-35-45_AD002", "M_18_1923_HUME-ST-5-7_AD007", "M_18_1924_HUME-ST-8-11_AD019"]

# Create fastqc files for each of the sequencing files
rule make_fastqc_pe:
    input:
        "paired_end_reads/M_17_426_1_AD08_{lane}_{direction}_001.fastq.gz"
    output:
        "paired_end_reads/meta/fastqc/M_17_426_1_AD08_{lane}_{direction}_001.fastqc.zip"
    shell:
        "fastqc -o paired_end_reads/meta/fastqc {input}"

# rule make_fastqc_mp:
#     input:
#         expand("mate_pair_libs/{mp_sample}_{lane}_{direction}_001.fastq.gz", mp_sample=MP_SAMPLES, lane=LANES, direction=DIRECTIONS)
#     threads:1
#     output:
#         expand("mate_pair_libs/meta/fastqc/{mp_sample}_{lane}_{direction}_001.fastqc.zip", mp_sample=MP_SAMPLES, lane=LANES, direction=DIRECTIONS)
#     shell:
#         "fastqc -o mate_pair_libs/meta/fastqc {input}"

# rule make_fastqc_mp:
#     input:
#         "mate_pair_libs/{mp_sample}_{lane}_{direction}_001.fastq.gz"
#     threads:1
#     output:
#         "mate_pair_libs/meta/fastqc/{mp_sample}_{lane}_{direction}_001.fastqc.zip"
#     shell:
#         "parallelfastqc -o mate_pair_libs/meta/fastqc {input}"

rule make_fastqc_mp:
    input:
        expand("mate_pair_libs/{mp_sample}_{lane}_{direction}_001.fastq.gz", mp_sample=MP_SAMPLES, lane=LANES, direction=DIRECTIONS)
    threads:20
    output:
        expand("mate_pair_libs/fastqc/{mp_sample}_{lane}_{direction}_001.fastqc.zip", mp_sample=MP_SAMPLES, lane=LANES, direction=DIRECTIONS),
        "mate_pair_libs/fastqc/fastqc_complete.txt"
    shell:
        "parallel -j {threads} fastqc -o mate_pair_libs/fastqc ::: {input}; "
        "touch mate_pair_libs/fastqc/fastqc_complete.txt"

def make_all_fastqcs(wildcards):
    return [
        "mate_pair_libs/meta/fastqc/M_18_1922_HUME-ST-35-45_AD002_Lane1_R1_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1923_HUME-ST-5-7_AD007_Lane1_R1_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1924_HUME-ST-8-11_AD019_Lane1_R1_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1922_HUME-ST-35-45_AD002_Lane1_R2_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1923_HUME-ST-5-7_AD007_Lane1_R2_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1924_HUME-ST-8-11_AD019_Lane1_R2_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1922_HUME-ST-35-45_AD002_Lane2_R1_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1923_HUME-ST-5-7_AD007_Lane2_R1_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1924_HUME-ST-8-11_AD019_Lane2_R1_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1922_HUME-ST-35-45_AD002_Lane2_R2_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1923_HUME-ST-5-7_AD007_Lane2_R2_001.fastq.gz",
        "mate_pair_libs/meta/fastqc/M_18_1924_HUME-ST-8-11_AD019_Lane2_R2_001.fastq.gz",
        "paired_end_reads/meta/fastqc/M_17_426_1_AD08_Lane1_R1_001.fastq.gz",
        "paired_end_reads/meta/fastqc/M_17_426_1_AD08_Lane1_R2_001.fastq.gz",
        "paired_end_reads/meta/fastqc/M_17_426_1_AD08_Lane2_R1_001.fastq.gz",
        "paired_end_reads/meta/fastqc/M_17_426_1_AD08_Lane2_R2_001.fastq.gz"
    ]

rule make_all_fastqcs:
    input:
        make_all_fastqcs
    output:
        pe_sum="paired_end_reads/meta/fastqc/all_complete.txt",
        mp_sum="mate_pair_libs/meta/fastqc/all_complete.txt"
    shell:
        "touch paired_end_reads/meta/fastqc/all_complete.txt mate_pair_libs/meta/fastqc/all_complete.txt"

