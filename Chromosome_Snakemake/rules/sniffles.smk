def get_T0_fasta(wildcards):
    return dictReads[wildcards.sample]["T0_fasta"]
def get_T1_long(wildcards):
    return dictReads[wildcards.sample]["T1_long"]


rule minimap2:
    """Run nucdiff."""
    input:
        get_T0_fasta,
        get_T1_long
    output:
        os.path.join(SNIFFLES, "{sample}.sam")
    conda:
        os.path.join('..', 'envs','minimap2.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        minimap2 -ax map-ont {input[0]} {input[1]} > {output[0]}
        """

rule sam_to_bam:
    """Run sam to bam."""
    input:
        os.path.join(SNIFFLES, "{sample}.sam")
    output:
        os.path.join(SNIFFLES, "{sample}.bam"),
        os.path.join(SNIFFLES, "sorted_{sample}.bam")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools view -S -b {input[0]} >  {output[0]}
        samtools sort {output[0]}  -o {output[1]}
        """

rule bam_index:
    """index bam."""
    input:
        os.path.join(SNIFFLES, "sorted_{sample}.bam")
    output:
        os.path.join(SNIFFLES, "sorted_{sample}.bam.bai")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools index {input[0]}     
        """

rule sniffles:
    """sniffles."""
    input:
        os.path.join(SNIFFLES, "sorted_{sample}.bam"),
        os.path.join(SNIFFLES, "sorted_{sample}.bam.bai")
    output:
        os.path.join(SNIFFLES, "{sample}_sniffles.vcf")
    conda:
        os.path.join('..', 'envs','sniffles.yaml')
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        sniffles -i {input[0]} -v {output[0]}
        """




rule aggr_sniffles:
    """Aggregate."""
    input:
        expand(os.path.join(SNIFFLES, "{sample}_sniffles.vcf"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_sniffles.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """

