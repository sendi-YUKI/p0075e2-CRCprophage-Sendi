process HOST_RECOMBINATION {
    cpus 6
    memory '12 GB'
    maxForks 1
    publishDir "${params.outdir}/host/recombination", mode:'copy'
    input:
    path config
    path assembly_manifest
    path source_fastas, stageAs:'source??/*'
    path helper
    output:
    path 'recombination', emit: result
    path 'recombination/masked', emit: masked
    path 'recombination/tree', emit: tree
    script:
    def fastaArgs = (source_fastas instanceof List ? source_fastas : [source_fastas]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 '${helper}' run --config '${config}' --assembly-manifest '${assembly_manifest}' \\
      --source-fastas ${fastaArgs} --cpus ${task.cpus} --outdir recombination_work
    cp -a recombination_work/deliverables recombination
    """
}
