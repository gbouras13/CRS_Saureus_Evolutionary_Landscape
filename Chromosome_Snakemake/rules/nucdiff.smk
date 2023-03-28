def get_T0_fasta(wildcards):
    return dictReads[wildcards.sample]["T0_fasta"]
def get_T1_fasta(wildcards):
    return dictReads[wildcards.sample]["T1_fasta"]


rule nucdiff:
    """Run nucdiff."""
    input:
        get_T0_fasta,
        get_T1_fasta
    output:
        os.path.join(NUCDIFF, "{sample}", "results", "{sample}_ref_struct.gff")
    conda:
        os.path.join('..', 'envs','nucdiff.yaml')
    params:
        os.path.join(NUCDIFF, "{sample}")
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        nucdiff {input[0]} {input[1]} {params[0]} {wildcards.sample}
        """

rule aggr_nucdiff:
    """Aggregate."""
    input:
        expand(os.path.join(NUCDIFF, "{sample}","results", "{sample}_ref_struct.gff"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_nucdiff.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """
