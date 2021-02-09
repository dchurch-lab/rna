import os
def reads_related_input(wildcards, path, ext,
                        seq_type = SEQUENCING_TYPE,
                        strand_annot = config["STRAND"]):
    if type(ext) == dict:
        ext = ext[seq_type]
    if seq_type == "PE":
        return expand(os.path.join(path, "{sample}_{strand}" + ext),
                       strand = strand_annot,
                       sample = wildcards.sample)
    elif seq_type == "SE":
        return expand(os.path.join(path, "{sample}" + ext),
                      sample = wildcards.sample)

def reads_input(wildcards,
                path = config['READS'],
                ext = ".fastq.gz"):
    return reads_related_input(**locals())

def fastqc_output(wildcards,
                  path = config["FASTQC"],
                  ext = ".fastqc.zip"):
    return reads_related_input(**locals())

def trimming_output(wildcards, path = config['TRIMMED'],
                    ext = {"SE": "_trimmed.fq.gz", "PE": "_val_{strand}.fq.gz"}):
    return reads_related_input(**locals())

def trimming_qc_output(wildcards, path = config['TRIMMED'],
                       ext = ".fastq.gz_trimming_report.txt"):
    return reads_related_input(**locals())
