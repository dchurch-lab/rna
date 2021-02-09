rule mark_duplicates:
    input:
        sats = os.path.join(config['ALIGNMENT'], '{sample}',
                            'Log.final.out'),
        success = os.path.join(config['ALIGNMENT'], '{sample}',
                               'sort.success')
    output:
        bam = os.path.join(config['ALIGNMENT'], '{sample}',
                           '{sample}_sorted_md.bam'),
        index = os.path.join(config['ALIGNMENT'], '{sample}',
                             '{sample}_sorted_md.bam.bai'),
        flagstat_metrics = os.path.join(config['ALIGNMENT_QUAL'],
                                        'flagstat_{sample}.txt'),
        md_metrics = os.path.join(config['ALIGNMENT_QUAL'],
                                 'DuplicationMetrics_{sample}.txt')
    params:
        tmp_dir = os.path.join(config['ALIGNMENT'],
                                '{sample}', "md_tmp") + '/',
        sorted_bam = os.path.join(config['ALIGNMENT'], '{sample}',
                                  'Aligned.sortedByCoord.out.bam')
    threads:
        1
    conda:
        config['CONDA_ALIGNMENT']
    log:
        os.path.join(config['LOGS'], 'MarkDuplicates', '{sample}.log')
    message:
        '\n##################### Marking Duplicates ###################\n'
        'Running Picard MarkDuplicates followed by indexing\n'
        '{output.bam}\n'
        '{output.index}\n'
        'Log file: {log}\n'
        '############################################################'
    shell:
        'picard -Xmx4g -Xms4g MarkDuplicates'
        ' I={params.sorted_bam} O={output.bam} ASSUME_SORTED=true'
        ' METRICS_FILE={output.md_metrics}'
        ' TMP_DIR={params.tmp_dir} 2> {log};'
        ' samtools index {output.bam};'
        ' samtools flagstat {output.bam} > {output.flagstat_metrics};'
        ' rm -f {params.sorted_bam}*;'
        ' rm -fr {params.tmp_dir}'

rule create_star_index:
    input:
        genome = config['GENOME'],
        annotations = config['ANNOTATIONS']
    output:
        expand(os.path.join(config['STAR_INDEX'], '{index_files}'),
               index_files = ['SA', 'SAindex', 'Genome', 'genomeParameters.txt',
                              'chrNameLength.txt', 'chrName.txt', 'chrStart.txt',
                              'exonGeTrInfo.tab', 'exonInfo.tab',
                              'geneInfo.tab', 'transcriptInfo.tab',
                              'sjdbInfo.txt', 'sjdbList.fromGTF.out.tab',
                              'sjdbList.out.tab'])
    params:
        index = config['STAR_INDEX'],
        overhang = 100
    threads:
        16
    conda:
        config['CONDA_ALIGNMENT']
    log:
        os.path.join(config["LOGS"], 'star_index.log')
    shell:
        'STAR --runThreadN {threads} --runMode genomeGenerate'
        ' --genomeDir {params.index}'
        ' --genomeFastaFiles {input.genome}'
        ' --sjdbGTFfile {input.annotations}'
        ' --sjdbOverhang $(({params.overhang} - 1)) &> {log}'

rule alignment:
    input:
        expand(os.path.join(config['STAR_INDEX'], '{index_files}'),
               index_files = ['SA', 'SAindex', 'Genome', 'genomeParameters.txt',
                              'chrNameLength.txt', 'chrName.txt', 'chrStart.txt',
                              'exonGeTrInfo.tab', 'exonInfo.tab',
                              'geneInfo.tab', 'transcriptInfo.tab',
                              'sjdbInfo.txt', 'sjdbList.fromGTF.out.tab',
                              'sjdbList.out.tab']),
        reads = trimming_output
    output:
        sats = os.path.join(config['ALIGNMENT'], '{sample}',
                            'Log.final.out'),
        success = os.path.join(config['ALIGNMENT'], '{sample}',
                               'sort.success')
    params:
        star_index = config['STAR_INDEX'],
        # STAR output
        unsorted_bam = os.path.join(config['ALIGNMENT'], '{sample}',
                                    'Aligned.out.bam'),
        out = os.path.join(config['ALIGNMENT'], '{sample}') + '/',
        bam = os.path.join(config['ALIGNMENT'], '{sample}',
                                   'Aligned.sortedByCoord.out.bam'),
        index = os.path.join(config['ALIGNMENT'], '{sample}',
                             'Aligned.sortedByCoord.out.bam.bai')
    threads:
        16
    conda:
        config['CONDA_ALIGNMENT']
    log:
        os.path.join(config['LOGS'], 'STAR', '{sample}.log')
    message:
        '\n######################### Mapping ##########################\n'
        'Running STAR to produce:\n'
        '{params.bam}\n'
        'Log file: {log}\n'
        '############################################################'
    shell:
        'STAR --runThreadN {threads} --genomeDir {params.star_index}'
        ' --readFilesIn {input.reads} '
        ' --outFileNamePrefix {params.out} --outSAMtype BAM Unsorted'
        ' --outStd Log --readFilesCommand "gunzip -c" > {log};'
        ' samtools sort -o {params.bam} -@ {threads} -T {params.out}'
        ' {params.unsorted_bam};'
        ' rm -f {params.unsorted_bam};'
        ' samtools index {params.bam};'
        ' touch {output.success}'
