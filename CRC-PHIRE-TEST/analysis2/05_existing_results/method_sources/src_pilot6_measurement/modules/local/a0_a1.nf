process PREPARE_GENOME {
    tag "${meta.id}"
    publishDir { "${params.outdir}/samples/${meta.id}/input" }, mode:'copy'
    input:
    tuple val(meta), path(fasta, stageAs: 'input.fna')
    path helper, stageAs: 'phageflow.py'
    output:
    tuple val(meta), path('prepared'), emit: prepared
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 phageflow.py prepare --genome-id '${meta.id}' --fasta '${fasta}' --outdir prepared
    """
}

process NORMALIZE_CALLS {
    tag "${meta.id}"
    publishDir { "${params.outdir}/samples/${meta.id}/normalized" }, mode:'copy'
    input:
    tuple val(meta), path(prepared), path(genomad)
    path helper, stageAs: 'phageflow.py'
    output:
    tuple val(meta), path('normalized'), emit: normalized
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 phageflow.py normalize --genome-id '${meta.id}' --fasta '${prepared}/genome.fna' --mapping '${prepared}/contig_map.tsv' --genomad-dir '${genomad}' --outdir normalized
    """
}

process ANNOTATE_QUALITY {
    tag "${meta.id}"
    publishDir { "${params.outdir}/samples/${meta.id}/final" }, mode:'copy'
    input:
    tuple val(meta), path(normalized), path(checkv)
    path helper, stageAs: 'phageflow.py'
    output:
    tuple val(meta), path('final'), emit: finished
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 phageflow.py quality --candidates '${normalized}/caller_hits.tsv' --checkv '${checkv}' --outdir final
    """
}

process EMPTY_QUALITY {
    tag "${meta.id}:zero_candidates"
    publishDir { "${params.outdir}/samples/${meta.id}/final" }, mode:'copy'
    input:
    tuple val(meta), path(normalized)
    path helper, stageAs: 'phageflow.py'
    output:
    tuple val(meta), path('final'), emit: finished
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 phageflow.py quality --candidates '${normalized}/caller_hits.tsv' --outdir final
    """
}

process COLLECT_REPORT {
    // Path.rglob() does not descend into child directory symlinks. Real copies
    // also prevent basename collisions between multiple samples' final dirs.
    stageInMode 'copy'
    publishDir "${params.outdir}", mode:'copy', saveAs: { filename -> filename.replaceFirst('report/','') }
    input:
    path sample_dirs, stageAs: 'samples??/*'
    path prepared_dirs, stageAs: 'prepared??/*'
    path helper, stageAs: 'phageflow.py'
    output:
    path 'report/*', emit: report
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 phageflow.py aggregate --inputs . --outdir report
    """
}
