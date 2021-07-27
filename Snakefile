# Bulk RNASeq
# author yhao at cshl dot edu
import pandas as pd
import os
import glob
import numpy as np
from itertools import combinations
from snakemake.utils import validate, min_version
from snakemake.shell import shell

min_version("5.1.2")

configfile: "config.yaml"
#validate(config, schema="config.schema.yaml")

threads_max = config["params"]["threads"]

def get_fq(wildcards):
    reads=SAMPLES.loc[wildcards.id, ["R1","R2"]].dropna(how="all")
    in_dir=config['fastqc']['fastq_dir']
    if(not in_dir or in_dir == './' ):
        return reads
    else:
        return [in_dir+r for r in reads]

def get_CONTRAST(conditions):
    cond=np.unique(conditions)
    cond=sorted(cond, key=str.casefold)
    b = list()
    for i in combinations(cond,2):
        b.append(i[0]+"_vs_"+i[1])
    return b

onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")

try:
    SAMPLES = pd.read_table(config["samples"]).set_index("sample", drop=False)
except:
    SAMPLES = pd.read_csv(config["samples"]).set_index("sample", drop=False)
else:
    print("There is an error with sample file.")

IDS = SAMPLES['sample']
READS = [os.path.basename(x).split(".")[0] for x in SAMPLES['R1'].dropna().tolist() + SAMPLES['R2'].dropna().tolist()]
CONTRASTS=get_CONTRAST(SAMPLES['condition'])

rule all:
    input:
        expand("gsea/gsea_{contrast}.log", contrast=CONTRASTS)

rule fastqc:
    input:
        get_fq
    threads: threads_max
    output:
        html=temp("fastqc/{id}.html"),
        zip=temp("fastqc/{id}.zip")
    params:
        prefix=config['fastqc']['out_dir']
    run:
        if not os.path.exists("fastqc"):
            os.makedirs("fastqc")

        shell("fastqc -t 4 --outdir {params.prefix} {input}")
        shell("touch {output}")

rule index_genome:
    input:
        fa=config['align']['genome_fasta'],
        gtf=config['align']['gtf']
    output:
        "genome/genomeLog.out"
    threads: threads_max
    params:
        ref=config['align']['reference'],
        starpath=config['align']['star_path'],
        fout="genome/genome",
        overhang=config['align']['star_sjdboverhang']
    run:
        if {params.ref}.pop() is None or not list({params.ref})[0].strip():
            if {input.fa}.pop() is None or not list({input.fa})[0].strip() or {input.gtf}.pop() is None or not list({input.gtf})[0].strip() or {params.overhang}.pop() is None or not str(list({params.overhang})[0]).strip():
                sys.exit('error: missing genome index!')
            else:
                # do genome indexing
                shell("echo \"Indexing genome ... \"")
                shell("rsem-prepare-reference --star "
                "--star-path {params.starpath} "
                "-p {threads} "
                "--gtf {input.gtf} "
                "--star-sjdboverhang {params.overhang} "
                "{input.fa} {params.fout}")
        else:
            shell("touch genome/genomeLog.out")

rule align:
    input:
        sp=get_fq,
        gn="genome/genomeLog.out"
    threads: threads_max
    output:
        result="mapping/{id}.genes.results",
        bam="mapping/{id}.STAR.genome.bam"
    params:
        fraglensd=config['align']['fragment_len_sd'],
        fraglenmean=config['align']['fragment_len_mean'],
        ref=config['align']['reference'],
        gtf=config['align']['gtf'],
        starpath=config['align']['star_path'],
        prefix="mapping/{id}"
    log:
        "mapping/{id}.align.log"
    run:
        if {params.ref}.pop() is None or not list({params.ref})[0].strip():
            this_ref="genome/genome"
        else:
            this_ref={params.ref}

        fc = len(unpack({input.sp})[0])
        if fc>1:
            #paired end
            shell("rsem-calculate-expression --star "
            "--star-path {params.starpath} "
            "-p {threads} "
            "--star-gzipped-read-file "
            "--star-output-genome-bam "
            "--paired-end "
            "{input.sp[0]} {input.sp[1]} "
            "{this_ref} {params.prefix}")
        elif fc == 1:
            #single end
            shell("rsem-calculate-expression --star "
            "--star-path {params.starpath} "
            "-p {threads} "
            "--star-gzipped-read-file "
            "--star-output-genome-bam "
            "--fragment-length-mean {params.fraglenmean} "
            "--fragment-length-sd {params.fraglensd} "
            "{input.sp[0]} {this_ref} {params.prefix}")
        else:
            print("Wrong number of input fastq. Exit.");exit(1)

rule get_exp_table:
    input:
        expand("mapping/{sample}.genes.results",sample=SAMPLES['sample'])
    output:
        expand("expression/{type}.txt", type=["expected_count","TPM","FPKM"])
    params:
        tpm="TPM", fpkm="FPKM", cnt="CNT"
    run:
        if not os.path.exists("expression"):
            os.makedirs("expression")
        if not os.path.exists("diffexp"):
            os.makedirs("diffexp")

        sel_cols=["expected_count", "TPM","FPKM"]
        all_files=glob.glob("mapping/*.genes.results")
        results = [pd.read_table(i).set_index("gene_id")[sel_cols] for i in all_files]

        for col in sel_cols:
            tab= pd.concat([ df[col] for df in results ],axis=1)
            new_cols=[ n.split(".genes.results")[0]+"_"+col for n in all_files]
            tab.columns = [n.split("/")[1] for n in new_cols]
            tab.to_csv('expression/{}.txt'.format(col), sep='\t', index=True)

rule deseq2:
    input:
        counts="expression/expected_count.txt"
    output:
        "diffexp/dds.rds","diffexp/diffexp.pdf",
         expand("diffexp/gsea_{contrast}.rnk", contrast=CONTRASTS)
    params:
        coldata=config['samples'],
        adjusted_pvalue=config['diffexp']['adjusted_pvalue'],
        rnk=config['diffexp']['rnk'],
        genome=config['diffexp']['genome'],
        gtf=config['align']['gtf']
    log:
        "diffexp/diffexp.log"
    script:
        "script/deseq2.R"

rule gsea:
    input:
        expand("diffexp/gsea_{contrast}.rnk", contrast=CONTRASTS)
    output:
        expand("gsea/gsea_{contrast}.log", contrast=CONTRASTS)
    params:
        jar=config['gsea']['exe_jar'],
        gmt=config['gsea']['gmt'],
        no_perm=config['gsea']['no_perm'],
        no_plot=config['gsea']['no_plot'],
        set_max=config['gsea']['set_max'],
        set_min=config['gsea']['set_min'],
        seed=config['gsea']['seed'],
    log:
        "gsea/gsea.log"
    script:
        "script/gsea_preranked.R"
