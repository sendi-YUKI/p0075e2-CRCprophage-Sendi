process CORE_KOFAM_DB_VALIDATE {
    cpus 1
    memory '2 GB'
    input:
    path database
    val expected_manifest_sha256
    path helper
    path common
    output:
    path 'kofam_db_check/database_verified.json', emit: receipt
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 ${helper} verify-db --database '${database}' --expected-manifest-sha256 '${expected_manifest_sha256}' --outdir kofam_db_check
    """
}

process CORE_KOFAM {
    tag "${meta.id}"
    publishDir "${params.outdir}/samples/${meta.id}", mode: 'copy'
    cpus 6
    memory '12 GB'
    maxForks 1
    input:
    tuple val(meta), path(genes)
    path database
    path database_receipt
    path helper
    path common
    path shim
    output:
    tuple val(meta), path('kofam'), emit: result
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 ${helper} run --genes-dir '${genes}' --database '${database}' \\
      --database-receipt '${database_receipt}' --shim '${shim}' --cpus ${task.cpus} --outdir kofam
    """
}

process CORE_KOFAM_SUMMARY {
    publishDir "${params.outdir}/functions", mode: 'copy'
    cpus 1
    memory '4 GB'
    input:
    path votu_members
    path(kofam_dirs), stageAs: 'kofam??/*'
    path(function_dirs), stageAs: 'functions??/*'
    path helper
    path common
    output:
    path 'kofam_summary', emit: summary
    script:
    def kd = (kofam_dirs instanceof List ? kofam_dirs : [kofam_dirs]).collect { "'${it}'" }.join(' ')
    def fd = (function_dirs instanceof List ? function_dirs : [function_dirs]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 ${helper} summarize --votu-members '${votu_members}' --kofam-dirs ${kd} \\
      --function-dirs ${fd} --outdir kofam_summary
    """
}
