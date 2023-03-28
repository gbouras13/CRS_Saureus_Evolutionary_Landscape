

rule abricate_card:
    """Run card."""
    input:
        os.path.join(INPUT, "{plasmid}.fasta")
    output:
        os.path.join(ABRICATE,"{plasmid}_card.tsv")
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
        """
        

rule aggr_abricate:
    """Aggregate."""
    input:
        expand(os.path.join(ABRICATE,"{plasmid}_card.tsv"), plasmid = PLASMIDS),
    output:
        os.path.join(LOGS, "aggr_abricate.txt")
    resources:
        mem_mb=SmallJobMem,
        time=2
    threads:
        1
    shell:
        """
        touch {output[0]}
        """

