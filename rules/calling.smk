
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
            -RGSM {wildcards.sample}
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
		"../envs/environment.yaml"
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

        
rule genomics_db_import:
    input:
        gvcfs=expand("called/{sample}.g.vcf.gz", sample = samples.names.values.tolist())
    output:
        db=directory("called/db"),
    log:
        "logs/gatk/genomicsdb.log"
    params:
        intervals="1",
        db_action="create", # optional
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=2048
    wrapper:
        "v0.75.0/bio/gatk/genomicsdbimport"
        
rule genotype_gvcfs:
    input:
        genomicsdb="called/db",  # combined gvcf over multiple samples
        ref=config["genome"]
    output:
        vcf="genotyped/all.vcf.gz"
    log:
        "logs/gatk/genotypegvcfs.log"
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    resources:
        mem_mb=1024
    wrapper:
        "v0.75.0/bio/gatk/genotypegvcfs"


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
        
rule download_vep_plugins:
    output:
        directory("resources/vep/plugins")
    params:
        release=100
    wrapper:
        "v0.75.0/bio/vep/plugins"
        
rule annotate_variants:
    input:
        calls=outputdir + "genotyped/all.vcf.gz",  # .vcf, .vcf.gz or .bcf
        cache="/dataVolume/.vep/homo_sapiens/102_GRCh37", # can be omitted if fasta and gff are specified
        plugins="resources/vep/plugins",
        # optionally add reference genome fasta
        # fasta="genome.fasta",
        # fai="genome.fasta.fai", # fasta index
        # gff="annotation.gff",
        # csi="annotation.gff.csi", # tabix index
    output:
        calls = outputdir + "genotyped/variants.annotated.bcf", # .vcf, .vcf.gz or .bcf
        stats = outputdir + "genotyped/variants.html"
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool"],
        extra="--everything"  # optional: extra arguments
    log:
        "logs/vep/annotate.log"
    threads: 4
    wrapper:
        "v0.75.0/bio/vep/annotate"
