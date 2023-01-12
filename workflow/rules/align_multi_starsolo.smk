# -*- coding: utf-8 -*-

""" align_multi_starsolo

Performs multimapping alignment using STARsolo, allowing for a high number of
multimappers for downstream reassignment using stellarscope

Setup:

SAMPLES is a dictionary with the sample name as the key. The value is a dictionary
with two keys, 'r1' and 'r2', with the values being a list of the FASTQ files belonging
to the sample. For example:


SAMPLES = {
    'sample1': {
        'r1': ['run919.R1.fastq.gz', 'run542.R1.fastq.gz', 'run265.R1.fastq.gz',],
        'r2': ['run919.R2.fastq.gz', 'run542.R2.fastq.gz', 'run265.R2.fastq.gz',]
    },
    'sample2': {
        'r1': ['run623.R1.fastq.gz', 'run983.R1.fastq.gz',],
        'r2': ['run623.R2.fastq.gz', 'run983.R2.fastq.gz',]
    },
}

has two samples. Sample 1 has three PE runs and Sample 2 has two PE runs. Note that the
order of file names should be consistent in 'r1' and 'r2'. This dictionary should be 
straightforward to create from SRA projects by using the BioSample (SAMN*) or GEO (GSM*) 
sample ID as the sample name and SRA Run accessions (SRR*) as the run names.

There are three variables that need to be set in `config.yaml`:

  + `star_index`: is the path to the STAR index
  + `whitelist`: is the whitelist for the chemistry
  + 'align_multi_starsolo': common command-line arguments for STAR. Each key will be 
     prepended with '--' and values are converted to strings.

Example of config.yaml

```
star_index: 'refs/indexes/STAR_gdc38_gencode38'
whitelist:
    v3: "resources/whitelist/whitelist.10x.v3.txt"
align_multi_starsolo:
    soloType: "CB_UMI_Simple"
    soloCBstart: 1
    soloCBlen: 16
    soloUMIstart: 17
    soloUMIlen: 12
    outSAMattributes: "NH HI nM AS CR UR CB UB GX GN sS sQ sM"
    outSAMtype: "BAM SortedByCoordinate"
    clipAdapterType: "CellRanger4"
    outFilterScoreMin: 30
    soloCBmatchWLtype: "1MM_multi_Nbase_pseudocounts"
    soloUMIfiltering: "MultiGeneUMI_CR"
    soloUMIdedup: "1MM_CR"
    limitOutSJcollapsed: 5000000
    outFilterMultimapScoreRange: 5
    outFilterMultimapNmax: 500
```

"""


rule align_multi_starsolo:
    output:
        "results/align_multi_starsolo/{samp}/Aligned.sortedByCoord.out.bam",
        "results/align_multi_starsolo/{samp}/Solo.out/Gene/filtered/barcodes.tsv"
    input:
        lambda wc: SAMPLES[wc.samp]['r1'] + SAMPLES[wc.samp]['r2']
    params:
        r1 = lambda wc: ','.join(SAMPLES[wc.samp]['r1']),
        r2 = lambda wc: ','.join(SAMPLES[wc.samp]['r2']),
        common_args = lambda wc: ' '.join(f'--{k} {v}' for k,v in config['align_multi_starsolo'].items())
    threads: 6
    conda: "../envs/star.yaml"
    shell:
        '''
mkdir -p $(dirname {output[0]})
STAR\
 --runThreadN {threads}\
 --genomeDir {config[star_index]}\
 --readFilesIn {params.r2} {params.r1}\
 --readFilesCommand zcat\
 --soloCBwhitelist {config[whitelist][v3]}\
 {params.common_args}\
 {params.msr_arg}\
 {params.mnm_arg}\
 --outFileNamePrefix $(dirname {output[0]})/
        '''