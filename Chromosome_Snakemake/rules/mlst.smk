
def get_T0_fasta(wildcards):
    return dictReads[wildcards.sample]["T0_fasta"]
def get_T1_fasta(wildcards):
    return dictReads[wildcards.sample]["T1_fasta"]

rule mlst:
    """Run mlst."""
    input:
        get_T0_fasta,
        get_T1_fasta
    output:
        os.path.join(MLST,"{sample}_T0.csv"),
        os.path.join(MLST,"{sample}_T1.csv")
    conda:
        os.path.join('..', 'envs','mlst.yaml')
    resources:
        mem_mb=4000,
        time=10
    threads:
        1
    shell:
        """
        mlst --scheme saureus --nopath  --csv {input[0]} > {output[0]}
        mlst --scheme saureus --nopath  --csv {input[1]} > {output[1]}
        """

rule aggr_mlst:
    """Aggregate."""
    input:
        expand(os.path.join(MLST,"{sample}_T0.csv"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_mlst.txt")
    resources:
        mem_mb=4000,
        time=2
    threads:
        1
    shell:
        """
        touch {output[0]}
        """

