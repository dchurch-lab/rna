#ruleorder: quantify_paired > quantify_single

rule create_salmon_index:
    input:
        gentrome = config['GENTROME'],
        decoys = config['DECOYS']
    output:
        os.path.join(config['SALMON_INDEX'], 'versionInfo.json')
    params:
        index = config['SALMON_INDEX'],
        kmer = 31
    threads:
        8
    conda:
        config['CONDA_ALIGNMENT']
    log:
        os.path.join(config["LOGS"], 'salmon_index.log')
    shell:
        'salmon index -t {input.gentrome} -i {params.index}'
        ' -k {params.kmer} --decoys {input.decoys} --threads {threads}'
        ' --gencode &> {log}'

rule quantify:
    input:
        os.path.join(config['SALMON_INDEX'], 'versionInfo.json'),
        reads = trimming_output
    output:
        os.path.join(config['QUANT'], '{sample}', 'quant.sf'),
        os.path.join(config['QUANT'], '{sample}', 'lib_format_counts.json')
    params:
        index = config['SALMON_INDEX'],
        out_dir = os.path.join(config['QUANT'], '{sample}'),
        extra = ''
    threads:
        8
    conda:
        config['CONDA_ALIGNMENT']
    log:
        os.path.join(config["LOGS"], 'salmon', '{sample}.log')
    message:
        '\n######################## Quantifing #########################\n'
        'Quantifying transcript abundance for:\n'
        '{input[0]}\n'
        'Log file: {log}\n'
        '############################################################'
    shell:
        """
        input_reads=({input.reads})
        if [ "${{#input_reads[@]}}" -eq 2 ]; then
            salmon quant \
             -i {params.index} \
             -l A --no-version-check \
             -1 ${{input_reads[0]}} -2 ${{input_reads[1]}} \
             -p {threads} \
             --validateMappings \
             -o {params.out_dir} \
             --gcBias {params.extra} &> {log}
        elif [ "${{#input_reads[@]}}" -eq 1 ]; then
            salmon quant \
             -i {params.index} \
             -l A --no-version-check \
             -r ${{input_reads[0]}} \
             -p {threads} \
             --validateMappings \
             -o {params.out_dir} \
             --gcBias {params.extra} &> {log}
        else
            echo "Wrong number of fastq files!"
            exit 1
        fi
        """
