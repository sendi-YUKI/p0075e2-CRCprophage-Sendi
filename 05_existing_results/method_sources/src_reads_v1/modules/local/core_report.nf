process CORE_REPORT {
    cpus 1
    memory '2 GB'
    container 'python:3.12.13-bookworm@sha256:3cd9086bdb30f7c9bc08a3fa621d9842e0d3f6f9291aeb4677e0547817c10b12'
    publishDir "${params.outdir}", mode:'copy', saveAs:{ name -> name.replaceFirst('core_report/','') }
    input:
    path a1, stageAs:'baseline'
    path catalog
    path qc
    path genes, stageAs:'genes??/*'
    path taxonomy, stageAs:'taxonomy??/*'
    path functions, stageAs:'functions??/*'
    path function_summary
    path caller_comparison
    path kofam_summary
    path(host_dirs), stageAs:'host??/*'
    path locus_outputs
    path research_outputs
    path host_source_outputs
    path helper, stageAs:'core_report.py'
    output:
    path 'core_report/*', emit: report
    script:
    def sourceArg = params.host_association_config ? "--host-source '${host_source_outputs}'" : ''
    def researchArg = params.catalog_policy ? "--locus '${locus_outputs}' --research-summary '${research_outputs}'" : ''
    def koArg = params.enable_kofam ? "--kofam-summary '${kofam_summary}'" : ''
    def hostArg = params.host_config ? '--host-dirs ' + (host_dirs instanceof List ? host_dirs : [host_dirs]).collect { "'${it}'" }.join(' ') : ''
    def geneArgs = (genes instanceof List ? genes : [genes]).collect { "'${it}'" }.join(' ')
    def taxArgs = (taxonomy instanceof List ? taxonomy : [taxonomy]).collect { "'${it}'" }.join(' ')
    def functionArgs = (functions instanceof List ? functions : [functions]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 core_report.py --a1 '${a1}' --catalog '${catalog}' --qc '${qc}' \\
      --genes ${geneArgs} --taxonomy ${taxArgs} --functions ${functionArgs} \\
      --function-summary '${function_summary}' --caller-comparison '${caller_comparison}' ${koArg} ${hostArg} ${researchArg} ${sourceArg} --outdir core_report
    """
}
