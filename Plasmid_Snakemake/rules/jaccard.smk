rule calc_jaccard:
    input:
        rtab = os.path.join(PANAROO, "gene_presence_absence.Rtab")
    output:
        matrix = os.path.join(RESULTS,"jaccard_matrix.csv")
    threads:
        BigJobCpu
    conda:
        os.path.join('..', 'envs','scripts.yaml')
    resources:
        mem_mb=BigJobMem
    script:
        '../scripts/calc_jaccard.py'


rule aggr_jaccard:
    input:
        os.path.join(RESULTS,"jaccard_matrix.csv"),
    output:
        os.path.join(LOGS, "aggr_jaccard.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """

