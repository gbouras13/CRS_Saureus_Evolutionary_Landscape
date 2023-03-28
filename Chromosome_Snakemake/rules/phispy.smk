rule phispy:
    """Run phispy."""
    input:
        get_T0_gbk,
        get_T1_gbk
    output:
        os.path.join(PHISPY,"{sample}_T0", "prophage_coordinates.tsv"),
        os.path.join(PHISPY,"{sample}_T0", "phage.fasta"),
        os.path.join(PHISPY,"{sample}_T1", "prophage_coordinates.tsv"),
        os.path.join(PHISPY,"{sample}_T1", "phage.fasta")
    conda:
        os.path.join('..', 'envs','phispy.yaml')
    params:
        os.path.join(PHISPY,"{sample}_T0"),
        os.path.join(PHISPY,"{sample}_T1")
    threads:
        8
    resources:
        mem_mb=16000,
        time=60
    shell:
        """
        phispy {input[0]} --output_choice 512 -o {params[0]} --phage_genes 0
        phispy {input[1]} --output_choice 512 -o {params[1]} --phage_genes 0
        """


rule aggr_phispy:
    """Aggregate."""
    input:
        expand(os.path.join(PHISPY,"{sample}_T0", "prophage_coordinates.tsv"), sample = SAMPLES),
        expand(os.path.join(PHISPY,"{sample}_T1", "phage.fasta"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_phispy.txt")
    threads:
        1
    resources:
        mem_mb=1000,
        time=1
    shell:
        """
        touch {output[0]}
        """

