
rule isescan:
    """Run isescan."""
    input:
        get_T0_fasta,
        get_T1_fasta
    output:
        os.path.join(ISESCAN, "{sample}_T0", "{sample}.touch"),
        os.path.join(ISESCAN, "{sample}_T1", "{sample}.touch")
    conda:
        os.path.join('..', 'envs','isescan.yaml')
    params:
        os.path.join(ISESCAN, "{sample}_T0"),
        os.path.join(ISESCAN, "{sample}_T1")
    threads:
        8
    resources:
        mem_mb=16000,
        time=60
    shell:
        """
        isescan.py --seqfile {input[0]} --output {params[0]} --nthread {threads}
        isescan.py --seqfile {input[1]} --output {params[1]} --nthread {threads}
        touch {output[0]}
        touch {output[1]}
        """



rule aggr_isescan:
    """Aggregate."""
    input:
        expand(os.path.join(ISESCAN, "{sample}_T0", "{sample}.touch"), sample = SAMPLES),
        expand(os.path.join(ISESCAN, "{sample}_T1", "{sample}.touch"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_isescan.txt")
    threads:
        1
    resources:
        mem_mb=1000,
        time=2
    shell:
        """
        touch {output[0]}
        """

