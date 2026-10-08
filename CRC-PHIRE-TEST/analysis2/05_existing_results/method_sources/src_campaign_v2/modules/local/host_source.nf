process HOST_SOURCE_ASSOCIATION {
    cpus 1
    memory '46 GB'
    maxForks 1
    publishDir "${params.outdir}/host_source", mode:'copy'
    input:
    path research
    path mask
    path tree
    path config
    path bundle
    path helper, stageAs:'host_source.py'
    path r_helper, stageAs:'host_source_models.R'
    path recombination_result
    path sensitivity_helper, stageAs:'host_sensitivity.R'
    output:
    path 'host_source', emit: result
    script:
    def recArgs = params.host_recombination_config ? "--recombination-result '${recombination_result}'" : ''
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 '${helper}' run --research '${research}' --mask '${mask}' --tree '${tree}' \\
      --config '${config}' --bundle '${bundle}' --r-script '${r_helper}' --sensitivity-script '${sensitivity_helper}' ${recArgs} --outdir host_source
    """
}
