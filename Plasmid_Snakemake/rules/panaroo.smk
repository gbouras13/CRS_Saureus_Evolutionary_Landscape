rule panaroo:
    """Run panaroo."""
    input:
        gffs = expand(os.path.join(PLASMID_GFFS,"{plasmid}.gff") , plasmid = PLASMIDS)
    output:
        os.path.join(PANAROO, "gene_presence_absence.Rtab")
    conda:
        os.path.join('..', 'envs','panaroo.yaml')
    params:
        out_dir = PANAROO
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        panaroo -i {input.gffs} -o {params.out_dir}   -t {threads}  --clean-mode sensitive 
        """

# C285.5 is causing issues so remove it 

rule aggr_panaroo:
    """Aggregate."""
    input:
        os.path.join(PANAROO, "gene_presence_absence.Rtab")
    output:
        os.path.join(LOGS, "aggr_panaroo.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem
    shell:
        """
        touch {output[0]}
        """

