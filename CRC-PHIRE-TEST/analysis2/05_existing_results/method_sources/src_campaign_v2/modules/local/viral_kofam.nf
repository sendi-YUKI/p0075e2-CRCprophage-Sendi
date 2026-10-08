process VIRAL_KOFAM {
    cpus 6
    memory '12 GB'
    maxForks 1
    publishDir "${params.outdir}/viral_kofam", mode:'copy'
    input:
    path genes
    path catalog
    path database
    path receipt
    path helpers
    output:
    path 'kofam', emit: result
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 core_kofam.py run --genes-dir '${genes}' --database '${database}' --database-receipt '${receipt}' \\
      --shim "\$PWD/kofam_hmmsearch_serial" --cpus ${task.cpus} --outdir kofam
    python3 viral_kofam.py --catalog '${catalog}' --genes '${genes}' --kofam kofam
    """
}
