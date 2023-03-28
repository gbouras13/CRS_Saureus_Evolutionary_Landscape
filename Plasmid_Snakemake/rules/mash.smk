rule run_mash:
    input:
        query = os.path.join(INPUT, "{plasmid}.fasta"),
        refs = expand( os.path.join(INPUT,"{plasmid}.fasta"), plasmid = PLASMIDS)
    output:
        mash = os.path.join(MASH,"{plasmid}_mash.txt")
    threads:
        BigJobCpu
    conda:
        os.path.join('..', 'envs','mash.yaml')
    resources:
        mem_mb=BigJobMem
    shell:
        """
        mash dist  {input.query} {input.refs} > {output.mash}
        """

rule mash_matrix:
    input:
        mashes = expand( os.path.join(MASH,"{plasmid}_mash.txt"), plasmid = PLASMIDS)
    output:
        matrix = os.path.join(RESULTS,"mash_matrix.csv")
    threads:
        BigJobCpu
    params:
        TMP,
        PLASDB
    conda:
        os.path.join('..', 'envs','scripts.yaml')
    resources:
        mem_mb=BigJobMem
    script:
        '../scripts/mash_matrix.py'

rule aggr_mash:
    input:
        os.path.join(RESULTS,"mash_matrix.csv"),
    output:
        os.path.join(LOGS, "aggr_mash.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """
