#!/usr/bin/env nextflow

// please run this file with this command code: nextflow run codes.nf --reads1 /path/to/SRR26244988_1.fastq.gz --reads2 /path/to/SRR26244988_2.fastq.gz

nextflow.enable.dsl=2

// Define default values for the input parameters
params.reads1 = null
params.reads2 = null

// Check if input parameters are provided
if (!params.reads1 || !params.reads2) {
    exit 1, "Please provide input FASTQ file paths using '--reads1' and '--reads2' parameters"
}

process skesa {
    input:
    file reads_ch1
    file reads_ch2

    output:
    path "assembled.fasta", emit: assembled_fasta
    path "skesa_stdout.txt", emit: skesa_stdout

    script:
    """
    mkdir -p results
    skesa --fastq ${reads_ch1} --fastq ${reads_ch2} --contigs_out assembled.fasta > skesa_stdout.txt 2>&1
    """
}

process quality_assessment {

    publishDir "results/quality_assessment", mode: 'copy', pattern: "*.txt"

    input:
    path assembled_fasta

    output:
    file "${assembled_fasta.baseName}_quast.txt"

    script:
    """
    quast "${assembled_fasta}" -o ${assembled_fasta.baseName}_quast_output
    mv ${assembled_fasta.baseName}_quast_output/report.txt ${assembled_fasta.baseName}_quast.txt
    """
}


process mlst {
    input:
    path assembled_fasta

    output:
    path "mlst_results.txt"

    publishDir "results/mlst", mode: 'symlink'

    script:
    """
    mlst $assembled_fasta > mlst_results.txt
    """
}



workflow {
    reads_ch1 = file(params.reads1)
    reads_ch2 = file(params.reads2)
    

    main:
    skesa_results = skesa(reads_ch1, reads_ch2)
    quality_assessment(skesa_results.assembled_fasta)
    mlst_results = mlst(skesa_results.assembled_fasta)
    
    
}