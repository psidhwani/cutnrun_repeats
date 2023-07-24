# === PATH OF SCRIPTS ===
# Paths are always the first obstacle to science, so let's start by
# dealing with them properly

#ml biology and ml bedtools prior to running
import os
shell.prefix("PATH=$PATH:"+config["scripts_path"]+"; ")

workdir: config["workdir"]

# === CONFIGURE PATH OF INPUT FILES ===
# SAMPLES = config["samples"]
# for key, value in SAMPLES.items():
#     SAMPLES[key]["PE1"]=srcdir(os.path.join(config["data_root"],value["PE1"]))
#     SAMPLES[key]["PE2"]=srcdir(os.path.join(config["data_root"],value["PE2"]))
# SAMPLESID=config.get("run_samples",list(SAMPLES.keys()))

#config['samples'].keys()

mymacs=list(config['macs'].keys())

OUTPUTS=[expand('{sampleID}/peaks/{aligner_config}/{rgP}/tracks/{macs}_{mtype}.uniq.bw', sampleID=config['samples'].keys(), aligner_config='chm13v2_unmasked',rgP='allPairs',macs=config['macs'].keys(), mtype=['FE','logLR'])]

wildcard_constraints:
    sampleID="[^/]+",
    macs="[^/]+",
    chip="[^/]+",
    bowtie_config="[^/]+",
    aligner_config="[^/]+",
    rgP="allPairs|unmergedPairs_pear",
    rgS="readthrough|mergedPairs_pear",#"([^/]+)/([^/]+)(/D|/R|)",
    fq="reads",
    mates_pairing_mode="[^/]+",
    track_type="[^/]+"


### =========== RULES START ===========================
rule all:
    input: OUTPUTS

# === LINK INPUT FILE INTO THE PIPELINE RUNNING FOLDER ===
# a symlink to the raw data is created in the pipeline running folders {sampleID}. If the "nreads" option is set (not -1), then instead of a symlink, the first "nreads" of the raw data are COPIED in the {sampleID} folder.
rule prepare_input:
    input:
        PE1=lambda wildcards: config['samples'][wildcards.sampleID][wildcards.chip]["PE1"],
        PE2=lambda wildcards: config['samples'][wildcards.sampleID][wildcards.chip]["PE2"]
    output:
        PE1='{sampleID}/{chip}/deduped/reads.1.fastq.gz',
        PE2='{sampleID}/{chip}/deduped/reads.2.fastq.gz'
    run:
        nlines=4*config["nreads"]
        if nlines<0:
            print("Running pipeline on ALL reads")
            shell("mkdir -p {wildcards.sampleID}")
            # if ("{input.PE1}".endswith(".gz")):
            shell("ln -s $(readlink -e {input.PE1}) {output.PE1}")
            # else:
            #     shell("gzip -c $(readlink -e {input.PE1})> {output.PE1}")
            # if ("{input.PE2}".endswith(".gz")):
            shell("ln -s $(readlink -e {input.PE2}) {output.PE2}")
            # else:
            #     shell("gzip -c $(readlink -e {input.PE2})> {output.PE2}")
        else:
            print("Running pipeline on subset of " + str(nlines) + " reads")
            shell("mkdir -p {wildcards.sampleID}")
            if ("{input.PE1}".endswith(".gz")):
                shell("set +o pipefail; pwd; zcat $(readlink -f {input.PE1}) | head -n " + str(nlines) + " | gzip > {output.PE1}")
                shell("set +o pipefail; pwd; zcat $(readlink -f {input.PE2}) | head -n " + str(nlines) + " | gzip > {output.PE2}")
            else:
                shell("set +o pipefail; pwd; cat $(readlink -f {input.PE1}) | head -n " + str(nlines) + " | gzip > {output.PE1}")
                shell("set +o pipefail; pwd; cat $(readlink -f {input.PE2}) | head -n " + str(nlines) + " | gzip > {output.PE2}")


#rule deduplicate:
 #   input:
  #      PE1='{sampleID}/{chip}/raw/reads.1.fastq',
   #     PE2='{sampleID}/{chip}/raw/reads.2.fastq'
    #output:
     #   PE1='{sampleID}/{chip}/deduped/reads.1.fastq.gz',
      #  PE2='{sampleID}/{chip}/deduped/reads.2.fastq.gz'
    #benchmark: 'benchmarks/{sampleID}_preprocessing_2_deduplicate.tsv'
    #log: '{sampleID}/{chip}/deduped/deduping.log'
    #params:
     #   xmx=config["max_mem"]
    #threads: config["max_threads"]
    #shell: """
#clumpify.sh t={threads} in1={input.PE1} in2={input.PE2} out1={output.PE1} out2={output.PE2} dedupe=t subs=0 reorder=f overwrite=t -Xmx{params.xmx} 2>{log}
#"""

# === REMOVE ADAPTERS AND TRIM OUT LOW QUALITY BASES===
# trims adapters
rule trim:
    input: 
        PE1='{sampleID}/{chip}/deduped/reads.1.fastq.gz',
        PE2='{sampleID}/{chip}/deduped/reads.2.fastq.gz'
    output: 
        PE1='{sampleID}/{chip}/deduped.trimmed/long/reads.1.fastq.gz',
        PE2='{sampleID}/{chip}/deduped.trimmed/long/reads.2.fastq.gz',
        readthrough='{sampleID}/{chip}/deduped.trimmed/readthrough/reads.fastq.gz'
    #benchmark: 'benchmarks/{sampleID}_preprocessing_1_trim.tsv'
    # log: '{sampleID}/chimeras/deduped.trimmed/trimming.log'
    params:
        trimmomatic_options=config["trimmomatic_options"],
        trimmomatic_filepath=srcdir(config["trimmomatic_filepath"]),
        trimfasta=lambda wildcards: config["samples"][wildcards.sampleID][wildcards.chip]["clip"],
        xmx=config["max_mem"]
    threads: config["max_threads"]
    shell: """
java -XX:Par{rgP}elGCThreads={threads} -Xmx{params.xmx} -jar {params.trimmomatic_filepath} PE -threads {threads} -phred33 {input.PE1} {input.PE2} {output.PE1} {output.readthrough}_U1.fastq.gz {output.PE2} {output.readthrough}_U2.fastq.gz ILLUMINACLIP:{params.trimfasta}:{params.trimmomatic_options} && gzip -c <(gunzip -c {output.readthrough}_U1.fastq.gz) <(gunzip -c {output.readthrough}_U2.fastq.gz) > {output.readthrough} && rm {output.readthrough}_U1.fastq.gz && rm {output.readthrough}_U2.fastq.gz
"""

rule trim_nodrop:
    input: 
        PE1='{sampleID}/{chip}/deduped/reads.1.fastq.gz',
        PE2='{sampleID}/{chip}/deduped/reads.2.fastq.gz'
    output: 
        PE1='{sampleID}/{chip}/deduped.trimmed/allPairs/reads.1.fastq.gz',
        PE2='{sampleID}/{chip}/deduped.trimmed/allPairs/reads.2.fastq.gz'
    #benchmark: 'benchmarks/{sampleID}_preprocessing_1_trim.tsv'
    # log: '{sampleID}/chimeras/deduped.trimmed/trimming.log'
    params:
        trimmomatic_options=config["trimmomatic_options_nodrop"],
        trimmomatic_filepath=srcdir(config["trimmomatic_filepath"]),
        trimfasta=lambda wildcards: config["samples"][wildcards.sampleID][wildcards.chip]["clip"],
        xmx=config["max_mem"]
    threads: config["max_threads"]
    shell: """
java -XX:ParallelGCThreads={threads} -Xmx{params.xmx} -jar {params.trimmomatic_filepath} PE -threads {threads} -phred33 {input.PE1} {input.PE2} {output.PE1} "$(dirname {output.PE1})"/tmp_U1.fastq.gz {output.PE2} "$(dirname {output.PE1})"/tmp_U2.fastq.gz ILLUMINACLIP:{params.trimfasta}:{params.trimmomatic_options}
"""

rule pair_reads:
    input:
        PE1='{sampleID}/{chip}/deduped.trimmed/{rgP}/reads.1.fastq.gz',
        PE2='{sampleID}/{chip}/deduped.trimmed/{rgP}/reads.2.fastq.gz'
    output:'{sampleID}/{chip}/deduped.trimmed/{rgP}/reads.fastq.gz'
    threads: config["max_threads"]
    shell: """
/home/groups/astraigh/software/bbmap/bbmerge-auto.sh in1={input.PE1} in2={input.PE2} out={output} rem k=62 extend2=50 ecct
"""
rule align_bowtie:
    input: '{sampleID}/{chip}/deduped.trimmed/{rgP}/reads.fastq.gz'
    params:
        bowtie_config=lambda wildcards: config["bowtie_configurations"][wildcards.aligner_config]['S']
    output:
        bam='{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/reads.bam',
        stats='{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/reads.stats.txt'
    # benchmark: '{sampleID}/benchmarks/align_dna.tsv'
    threads: config["max_threads"]
    # group: "align"
    shell: """
bowtie2 -p{threads} {params.bowtie_config} -U {input} 2>{output.stats} | samtools view -bh -F2308 | samtools sort > {output.bam}
"""

rule align_bowtie_ecoli:
    input: '{sampleID}/{chip}/deduped.trimmed/{rgP}/reads.fastq.gz'
    output:
        bam_ecoli='{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/ecolireads.bam',
        stats_ecoli= '{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/ecolireads.stats.txt'
    threads: config["max_threads"]
    shell: """
bowtie2 --very-sensitive -x /oak/stanford/groups/astraigh/charseq2.0/genomes/ecoli/Ecoli -U {input} 2>{output.stats_ecoli} | samtools view -bh >{output.bam_ecoli}
"""

rule calculate_norm_factor:
    input: '{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/ecolireads.stats.txt'
    output: '{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/normfactor.txt'
    threads: config["max_threads"]
    shell: """
gawk '{{print $1}}' FPAT='[0-9]+' {input} | xargs | awk '{{print ($2/200)/($4 + $5)}}' > {output}
"""

rule scale_data:
    input:
        bam='{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/reads.bam',
        norm_factor='{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/normfactor.txt',
         
    output:'{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/reads.norm.bdg'
    threads: config["max_threads"]
    shell: """
bedtools genomecov -bg -scale `awk '{{print $1}}'{input.norm_factor}` -ibam {input.bam} > {output}
"""
# rule make_genomecov:
#     input:
#         bam='{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/{fq}.PE.bam'
#     output:
#         bg='{sampleID}/{chip}/bigwigs/{aligner_config}/{rgP}/genomecov.bg'
#     shell:"""    
# bedtools genomecov -bg -ibam {input.bam}> {output.bedgraph}
# """

# rule make_binned_bw:
#     input:
#         bg='{sampleID}/{chip}/bigwigs/{aligner_config}/{rgP}/genomecov.bg'
#     output:
#         raw_bg='{sampleID}/{chip}/bigwigs/{aligner_config}/{rgP}/binned_{binsize}.bg',
#         normalized_bg='{sampleID}/{chip}/bigwigs/{aligner_config}/{rgP}/binned_{binsize}.bg'

#     shell:"""    
# bedtools map -c 4 -null 0 -a /oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/genome.{binsize}.bed -b {input.bg} > {output.rawbg};
# awk -F $'\t' 'BEGIN{OFS=FS}{h=($5==0 ? $4*100 : $4/$5*100); print $1, $2, $3, h}' <(paste {output.raw_bg} <(cut -f4,4 /oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/dpnII.{binsize}.bed)) | sort -k1,1 -k2,2n > {output.normalized_bg}
# """


# def get_peaks_file(wildards):
#     y="summits.bed"
#     # if ("broad" in config["macs"][wildcards.macs]['options']):
#     #     y="broadPeaks"
#     return y

rule macs_PE_narrow:
    input:
        bamT=lambda wildcards: '{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/reads.norm.bdg'.format(sampleID=wildcards.sampleID, aligner_config=wildcards.aligner_config, rgP=wildcards.rgP, chip=config['macs'][wildcards.macs]['treatment']),
        bamC=lambda wildcards: '{sampleID}/{chip}/bowtie/{aligner_config}/{rgP}/reads.norm.bdg'.format(sampleID=wildcards.sampleID, aligner_config=wildcards.aligner_config, rgP=wildcards.rgP, chip=config['macs'][wildcards.macs]['control'])
    params:
        macs_options=lambda wildcards: config["macs"][wildcards.macs]['options'],
        peaks_file="summits.bed" #lambda wildcards: get_peaks_file(wildards)
    output:
        #peaks='{sampleID}/peaks/{aligner_config}/{rgP}/{macs}_peaks.peak',
        spmr_c='{sampleID}/peaks/{aligner_config}/{rgP}/{macs}_control_lambda.bdg',
        spmr_t='{sampleID}/peaks/{aligner_config}/{rgP}/{macs}_treat_pileup.bdg'
    # benchmark: '{sampleID}/benchmarks/align_dna.tsv'
    threads: config["max_threads"]
    # group: "align"
    shell: """
/home/groups/astraigh/miniconda3_py38_4.9.2/envs/macs2/bin/macs2 callpeak -B --SPMR -t {input.bamT} -c {input.bamC} -n {wildcards.macs} --outdir $(dirname {output.spmr_c}) -f BED {params.macs_options}
"""
#source activate macs; 
#ln -s $(readlink -e $(dirname {output.peaks}/{wildcards.macs}_{params.peaks_file})) {output.peaks} 

rule generate_FE_tracks:
    input:
        spmr_c='{sampleID}/peaks/{aligner_config}/{rgP}/{macs}_control_lambda.bdg',
        spmr_t='{sampleID}/peaks/{aligner_config}/{rgP}/{macs}_treat_pileup.bdg'
    output:
        FE='{sampleID}/peaks/{aligner_config}/{rgP}/tracks/{macs}_FE.bdg',
        logLR='{sampleID}/peaks/{aligner_config}/{rgP}/tracks/{macs}_logLR.bdg'
    params: 
        chromosomes=config['chromosomes']
    threads: 1
    # group: "align"
    shell: """
/home/groups/astraigh/miniconda3_py38_4.9.2/envs/macs2/bin/macs2 bdgcmp -t {input.spmr_t} -c {input.spmr_c} -o {output.FE} -m FE ; /home/groups/astraigh/miniconda3_py38_4.9.2/envs/macs2/bin/macs2 bdgcmp -t {input.spmr_t} -c {input.spmr_c} -o {output.logLR} -m logLR -p 0.00001
"""

rule sort_bdg:
    input:'{mytrack}.bdg'
    output: '{mytrack}.sorted.bdg'
    threads: 1
    shell: """
export LC_COLLATE=C; sort -k1,1 -k2,2n {input} >{output}
"""

rule intersect_uniq_100mers:
    input:'{mytrack}.sorted.bdg'
    output:'{mytrack}.uniq.sorted.bdg'
    threads: 1
    shell: """
bedtools intersect -a {input} -b /oak/stanford/groups/astraigh/Pragya/T2T/uniq_kmers/chm13v2.0_100mers.bed > {output}
"""

rule generate_bigwigs:
    input: '{mytrack}.uniq.sorted.bdg'
    output: '{mytrack}.uniq.bw'
    # benchmark: '{sampleID}/benchmarks/align_dna.tsv'
    params: 
        chromosomes=config['chromosomes']
    threads: 1
    # group: "align"
    shell: """
bedGraphToBigWig {input} {params.chromosomes} {output}
"""
