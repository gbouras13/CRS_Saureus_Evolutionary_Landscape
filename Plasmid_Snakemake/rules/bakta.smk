rule concat_plasmids:
    input:
        list = expand(os.path.join(INPUT, "{plasmid}.fasta") , plasmid = PLASMIDS)
    output:
        train = os.path.join(BAKTA, "all_plasmids.fa")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        cat {input.list} > {output.train}
        """

rule progidal_train:
    input:
        os.path.join(BAKTA, "all_plasmids.fa")
    output:
        os.path.join(BAKTA, "train.trn")
    threads:
        1
    conda:
        os.path.join('..', 'envs','prodigal.yaml')
    resources:
        mem_mb=BigJobMem
    shell:
        """
        prodigal -i {input[0]} -t {output[0]}
        """



### https://github.com/tseemann/prokka/issues/203 
### metagenome mode giving weird results for prokka which was breaking panaroo
### so train my own prodigal file instead after concatenating the plasmids
### edit line 714 of prokka 


rule bakta:
    """Run bakta."""
    input:
        os.path.join(INPUT, "{plasmid}.fasta"),
        os.path.join(BAKTA, "train.trn")
    output:
        os.path.join(BAKTA,"{plasmid}","{plasmid}.gff3"),
        os.path.join(BAKTA,"{plasmid}","{plasmid}.ffn"),
        os.path.join(BAKTA,"{plasmid}","{plasmid}.gbff")
    conda:
        os.path.join('..', 'envs','bakta.yaml')
    params:
        BAKTA_DB,
        os.path.join(BAKTA, "{plasmid}")
    threads:
        8
    resources:
        mem_mb=4000
    shell:
        """
        bakta --db {params[0]} --verbose --output {params[1]} --prefix {wildcards.plasmid} --locus-tag {wildcards.plasmid} --prodigal-tf {input[1]} --threads {threads} {input[0]} 
        """

rule move_gff:
    input:
        os.path.join(BAKTA,"{plasmid}","{plasmid}.gff3")
    output:
        os.path.join(PLASMID_GFFS,"{plasmid}.gff")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        cp {input[0]} {output[0]} 
        """
    
rule aggr_bakta:
    input:
        expand(os.path.join(PLASMID_GFFS,"{plasmid}.gff") , plasmid = PLASMIDS)
    output:
        os.path.join(LOGS, "aggr_bakta.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """
