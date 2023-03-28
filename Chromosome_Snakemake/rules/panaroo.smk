rule panaroo:
    """Run panaroo."""
    input:
        get_all_T0_gffs,
        get_all_T1_gffs
    output:
        os.path.join(PANAROO, "gene_presence_absence.Rtab"),
        os.path.join(PANAROO, "core_gene_alignment.aln")
    conda:
        os.path.join('..', 'envs','panaroo.yaml')
    params:
        PANAROO
    threads:
        32
    resources:
        mem_mb=50000,
        time=1000
    shell:
        """
        panaroo -i {input[0]} {input[1]} -o {params[0]} -a core --mode strict --core_threshold 0.98 -t {threads} 
        """

rule aggr_panaroo:
    """Aggregate."""
    input:
        os.path.join(PANAROO, "gene_presence_absence.Rtab"),
        os.path.join(PANAROO, "core_gene_alignment.aln")
    output:
        os.path.join(FLAGS, "aggr_panaroo.txt")
    threads:
        1
    resources:
        mem_mb=1000,
        time=2
    shell:
        """
        touch {output[0]}
        """

