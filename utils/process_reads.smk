ruleorder: trimming_paired > trimming_single
ruleorder: fastqc_paired > fastqc_single

rule trimming_paired:
    input:
        reads_input
    output:
        os.path.join(config['TRIMMED'], '{sample}_1_val_1.fq.gz'),
        os.path.join(config['TRIMMED'], '{sample}_2_val_2.fq.gz'),
        os.path.join(config['TRIMMED'],
                     '{sample}_1.fastq.gz_trimming_report.txt'),
        os.path.join(config['TRIMMED'],
                     '{sample}_2.fastq.gz_trimming_report.txt')
    params:
        qc_outdir = config["FASTQC"],
        outdir = config['TRIMMED']
    threads:
        1
    conda:
        config['CONDA_QUALITY']
    log:
        os.path.join(config["LOGS"], 'trim_galore', '{sample}.log')
    message:
        '\n######################## Trimming PE #######################\n'
        'Running trim_galore on: {wildcards.sample}\n'
        '############################################################'
    shell:
        'trim_galore --fastqc_args "--outdir {params.qc_outdir}"'
        ' --gzip -o {params.outdir} --paired {input} &> {log}'

rule trimming_single:
    input:
        reads_input
    output:
        os.path.join(config['TRIMMED'], '{sample}_trimmed.fq.gz'),
        os.path.join(config['TRIMMED'],
                     '{sample}.fastq.gz_trimming_report.txt')
    params:
        qc_outdir = config["FASTQC"],
        outdir = config['TRIMMED']
    threads:
        1
    conda:
        config['CONDA_QUALITY']
    log:
        os.path.join(config["LOGS"], 'trim_galore', '{sample}.log')
    message:
        '\n####################### Trimming SE ########################\n'
        'Running trim_galore on: {input}\n'
        '############################################################'
    shell:
        'trim_galore --fastqc_args "--outdir {params.qc_outdir}"'
        ' --gzip -o {params.outdir} {input} &> {log}'

rule fastqc_paired:
    input:
        reads_input
    output:
        os.path.join(config['FASTQC'], '{sample}_1_fastqc.zip'),
        os.path.join(config['FASTQC'], '{sample}_2_fastqc.zip')
    params:
        outdir = config["FASTQC"]
    threads:
        1
    conda:
        config["CONDA_QUALITY"]
    log:
        os.path.join(config["LOGS"], "fastqc", "{sample}.log")
    shell:
        "fastqc {input[0]} {input[1]} -o {params.outdir} &> {log}"

rule fastqc_single:
    input:
        reads_input
    output:
        os.path.join(config['FASTQC'], '{sample}_fastqc.zip')
    params:
        outdir = config["FASTQC"]
    threads:
        1
    conda:
        config["CONDA_QUALITY"]
    log:
        os.path.join(config["LOGS"], "fastqc", "{sample}.log")
    shell:
        "fastqc {input} -o {params.outdir} &> {log}"
