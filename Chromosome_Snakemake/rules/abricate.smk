

rule abricate_card:
    """Run card."""
    input:
        get_T0_fasta,
        get_T1_fasta
    output:
        os.path.join(ABRICATE,"{sample}_T0_card.tsv"),
        os.path.join(ABRICATE,"{sample}_T1_card.tsv")
    conda:
        os.path.join('..', 'envs','abricate.yaml')
    threads:
        1
    resources:
        mem_mb=4000,
        time=5
    shell:
        """
        abricate {input[0]} --db card > {output[0]}
        abricate {input[1]} --db card > {output[1]}
        """

rule abricate_vfdb:
    """Run vfdb."""
    input:
        get_T0_fasta,
        get_T1_fasta
    output:
        os.path.join(ABRICATE,"{sample}_T0_vfdb.tsv"),
        os.path.join(ABRICATE,"{sample}_T1_vfdb.tsv")
    conda:
        os.path.join('..', 'envs','abricate.yaml')
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=5
    shell:
        """
        abricate {input[0]} --db vfdb > {output[0]}
        abricate {input[1]} --db vfdb > {output[1]}
        """

rule aggr_abricate:
    """Aggregate."""
    input:
        expand(os.path.join(ABRICATE,"{sample}_T0_vfdb.tsv"), sample = SAMPLES),
        expand(os.path.join(ABRICATE,"{sample}_T0_card.tsv"), sample = SAMPLES),
    output:
        os.path.join(FLAGS, "aggr_abricate.txt")
    resources:
        mem_mb=SmallJobMem,
        time=2
    threads:
        1
    shell:
        """
        touch {output[0]}
        """

