
def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{sample}-{unit}.bam",
                  sample=wildcards.sample,
                  unit=units.loc[wildcards.sample].unit)
                  
def get_call_variants_params(wildcards, input):
    return (config["params"]["gatk"]["HaplotypeCaller"])

def get_fai():
    return config["genome"] + ".fai"

rule picard_readgroups:
    input:
        bam = outputdir + "HISAT2/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        rg_bam = outputdir + "Picard/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    log:
        "logs/gatk/Picard/{sample}.log"
    params:
        extra=config["params"]["gatk"]["HaplotypeCaller"]
    shell:
        """
        picard AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.rg_bam} \
            -RGID 4 \
            -RGLB lib1 \
            -RGPL ILLUMINA \
            -RGPU unit1 \
            -RGSM 20
        """

## Index bam files
rule bamindexpicard:
	input:
		bam = outputdir + "Picard/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		bai = outputdir + "Picard/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
	log:
		outputdir + "logs/samtools_index_{sample}.log"
	benchmark:
		outputdir + "benchmarks/samtools_index_{sample}.txt"
	conda:
		"envs/environment.yaml"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"

rule call_variants:
    input:
        bam = outputdir + "Picard/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        bai = outputdir + "Picard/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
        ref=config["genome"],
        known=config["known-variants"],
        regions=[]
    output:
        gvcf="called/{sample}.g.vcf.gz"
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra=config["params"]["gatk"]["HaplotypeCaller"]
    shell:
        "echo {input.bai}; "
        "gatk HaplotypeCaller {params.extra} --reference {input.ref} --input {input.bam} --output {output.gvcf} -ERC GVCF"

rule combine_calls:
    input:
        ref=config["genome"],
        gvcfs=expand("called/{sample}.g.vcf.gz", sample = samples.names.values.tolist())
    output:
        gvcf="called/all.g.vcf.gz"
    log:
        "logs/gatk/combinegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["genome"],
        gvcf="called/all.g.vcf.gz"
    output:
        vcf=temp("genotyped/all.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotypegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below
        vcfs=lambda w: expand("genotyped/all.vcf.gz", contig=get_contigs()),
    output:
        vcf=outputdir + "genotyped/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    wrapper:
        "0.40.2/bio/picard/mergevcfs"
