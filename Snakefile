"""
Author: Katarzyna Kedzierska
Affiliation: WCHG, Unieversity of Oxford
Aim: A simple Snakemake workflow for processing RNA-seq data.
Date: July 2019
Run: snakemake --use-conda -j NUM_OF_THREADS

"""

from glob import glob
import re
import os

configfile: "config.yaml"
SEQUENCING_TYPE = config["SEQUENCING_TYPE"].upper()

SAMPLES_DICT = dict()

with open(config['SAMPLE_SHEET'], "r+") as fil:
    next(fil) # skip the header
    for lin in fil.readlines():
        row = lin.strip("\n").split("\t")
        sample_id = row[0]
        sample_name = row[1]
        if sample_name in SAMPLES_DICT.keys():
            SAMPLES_DICT[sample_name].append(sample_id)
        else:
            SAMPLES_DICT[sample_name] = [sample_id]

SAMPLES = list(SAMPLES_DICT.keys())
SAMPLE_IDS = [sample_id for sample in SAMPLES_DICT.values()
                for sample_id in sample]

include: "utils/helper.smk"

rule all:
    input:
        # Salmon quanitification
        expand(os.path.join(config['QUANT'], '{sample}', '{result}'),
               sample = SAMPLE_IDS,
               result = ['lib_format_counts.json', 'quant.sf']),

        # STAR alignment
        expand(os.path.join(config['ALIGNMENT'], '{sample}', '{sample}_{result}'),
               sample = SAMPLE_IDS,
               result = ['sorted_md.bam',
                         'sorted_md.bam.bai']),

        # multiqc report
        "multiqc_report.html"

    message:
        '\n#################### RNA-seq pipeline #####################\n'
        'Running all necessary rules to produce complete output.\n'
        '############################################################'

include: "utils/process_reads.smk"
include: "utils/align.smk"
include: "utils/quantify.smk"
include: "utils/report.smk"
