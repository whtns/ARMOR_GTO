configfile: "config.yaml"
# configfile: "config/samples.yaml"

def get_trimmed_reads(wildcards):
  """Get trimmed reads of given sample-unit."""
  # paired-end sample
  return expand([outputdir + "FASTQtrimmed/{sample}_R1_val_1.fq.gz", outputdir + "FASTQtrimmed/{sample}_R2_val_2.fq.gz"], **wildcards)

rule bwa_mem_hg38:
    input:
        reads = get_trimmed_reads
    output:
        bwa_bam = outputdir + "bwa/hg38/{sample}.sorted.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=config["genome"],
        sort="samtools",
        sort_order="coordinate"
    threads: config["ncores"]
    wrapper:
        "0.27.1/bio/bwa/mem"
        
## Index bam files
rule bamindexbwa:
	input:
		bwa_bam = outputdir + "bwa/hg38/{sample}.sorted.bam"
	output:
		bam_index = outputdir + "bwa/hg38/{sample}.sorted.bam.bai"
	log:
		outputdir + "logs/samtools_index_{sample}.log"
	benchmark:
		outputdir + "benchmarks/samtools_index_{sample}.txt"
	conda:
		"../envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bwa_bam}"



rule read_counter:
	input:
	  bwa_bam = outputdir + "bwa/hg38/{sample}.sorted.bam",
	  bam_index = outputdir + "bwa/hg38/{sample}.sorted.bam.bai"
	output:
		outputdir + "ichorCNA/readDepth/{sample}.bin" + str(config["binSize"])+".wig"		
	params:
		readCounter=config["readCounterScript"],
		binSize=config["binSize"],
		qual="20",
		chrs=config["chrs"]
	resources:
		mem=4
	log:
		"logs/ichorCNA/readDepth/{sample}.bin" + str(config["binSize"])+".log"
	shell:
		"{params.readCounter} {input.bwa_bam} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"
		

rule ichorCNA:
	input:
		tum=outputdir + "ichorCNA/readDepth/{sample}.bin" + str(config["binSize"]) + ".wig",
		script = "scripts/run_ichorcna.R"
		#norm=lambda wildcards: "output/ichorCNA/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
	output:
		#corrDepth=outputdir + "ichorCNA/{sample}/{sample}.correctedDepth.txt",
		#param=outputdir + "ichorCNA/{sample}/{sample}.params.txt",
		cna=outputdir + "ichorCNA/{sample}/{sample}.cna.seg",
		#segTxt=outputdir + "ichorCNA/{sample}/{sample}.seg.txt",
		#seg=outputdir + "ichorCNA/{sample}/{sample}.seg",
		#rdata=outputdir + "ichorCNA/{sample}/{sample}.RData"
	params:
		outDir=outputdir + "ichorCNA/{sample}/",
	resources:
		mem=4
	log:
		"logs/ichorCNA/{sample}.log"
	shell:
		'''{Rbin} CMD BATCH --no-restore --no-save "--args samplewig='{input.tum}' outDir='{params.outDir}'" {input.script} {log}'''


