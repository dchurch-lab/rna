rule multiqc:
    input:
        # FastQC reports
        expand(os.path.join(config['FASTQC'], '{sample}{tmp}_fastqc.zip'),
               sample = SAMPLE_IDS,
               tmp = ["_1", "_2"] if SEQUENCING_TYPE == "PE" else [""]),

        # Mark Duplicates metrcis
        expand(os.path.join(config['ALIGNMENT_QUAL'],
                            'DuplicationMetrics_{sample}.txt'),
               sample = SAMPLE_IDS),

        # Alignment Stats
        expand(os.path.join(config['ALIGNMENT_QUAL'],
                            'flagstat_{sample}.txt'),
               sample = SAMPLE_IDS),

        # Salmon quanitification
        expand(os.path.join(config['QUANT'], '{sample}', '{result}'),
               sample = SAMPLE_IDS,
               result = ['lib_format_counts.json']),

        # STAR alignment
        expand(os.path.join(config['ALIGNMENT'], '{sample}',
               'Log.final.out'),
               sample = SAMPLE_IDS)

    output:
        "multiqc_report.html"
    threads:
        1
    conda:
        config['CONDA_QUALITY']
    log:
        os.path.join(config["LOGS"], 'multiqc.log')
    message:
        '\n######################### Multiqc ##########################\n'
        'Running multiqc on all intermediate\n'
        'quality checks.\n'
        '############################################################'
    shell:
        """
        export LC_ALL=en_US.utf-8
        export LANG=en_US.utf-8
        multiqc --force --ignore .snakemake/ -d -dd 1 . &> {log}
        """
