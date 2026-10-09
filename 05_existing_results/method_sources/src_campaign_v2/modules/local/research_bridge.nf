process RESEARCH_V1_BRIDGE {
    cpus 6
    memory '8 GB'
    maxForks 1
    publishDir "${params.outdir}/research_bridge", mode:'copy'
    input:
    path research, stageAs:'research_input'
    path independent, stageAs:'independent_input'
    path source_manifest, stageAs:'source_manifest.json'
    path source_fastas, stageAs:'sources??/*'
    path helpers
    output:
    path 'catalog', emit: catalog
    script:
    def ff = (source_fastas instanceof List ? source_fastas : [source_fastas]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 research_bridge.py --research-catalog '${research}' --independent-catalog '${independent}' \\
      --source-manifest '${source_manifest}' --source-fastas ${ff} --cpus ${task.cpus} --outdir catalog
    """
}
