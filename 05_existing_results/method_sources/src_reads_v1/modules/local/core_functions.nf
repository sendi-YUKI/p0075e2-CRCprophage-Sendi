process CORE_FUNCTIONS {
    tag "${meta.id}"
    publishDir "${params.outdir}/samples/${meta.id}", mode: 'copy'
    cpus 6
    memory '12 GB'
    maxForks 1
    container 'sha256:5bcc494406a8137d2103870fa7d3e71e46ec357c3b9af6b2e6ff770090b603f4'
    input:
    tuple val(meta), path(genes_dir)
    path master
    path functions_db, stageAs: 'functions_db'
    path helper
    output:
    tuple val(meta), path('functions'), emit: functions
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 ${helper} run --genes-dir ${genes_dir} --master ${master} \
      --database-root ${functions_db} --outdir functions --cpus ${task.cpus} ${params.conjscan_models_dir ? "--conjscan-models-dir '${params.conjscan_models_dir}'" : ""} ${params.element_pfam_db ? "--pfam-db '${params.element_pfam_db}'" : ""}
    """
}

process CORE_FUNCTION_SUMMARY {
    publishDir "${params.outdir}/functions", mode: 'copy'
    stageInMode 'copy'
    cpus 1
    memory '2 GB'
    container 'sha256:5bcc494406a8137d2103870fa7d3e71e46ec357c3b9af6b2e6ff770090b603f4'
    input:
    path votu_members
    path(function_dirs), stageAs: 'samples??/*'
    path helper
    output:
    path 'function_summary', emit: summary
    script:
    def dirs = (function_dirs instanceof List ? function_dirs : [function_dirs]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 ${helper} summarize --votu-members ${votu_members} \
      --functions-dirs ${dirs} --outdir function_summary
    """
}

process CORE_ELEMENT_COLLECT {
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/element_classification", mode: 'copy'
    stageInMode 'copy'
    input:
    path master
    path function_dirs, stageAs:'samples??/*'
    path helper
    output:
    path 'elements', emit: elements
    script:
    def dirs = (function_dirs instanceof List ? function_dirs : [function_dirs]).collect { "'${it}'" }.join(' ')
    """
    python3 ${helper} collect-elements --master '${master}' --functions-dirs ${dirs} --outdir elements
    """
}
