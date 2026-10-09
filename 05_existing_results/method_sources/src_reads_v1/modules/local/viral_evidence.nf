process VIRAL_BACPHLIP {
    cpus 6
    memory '12 GB'
    maxForks 1
    publishDir "${params.outdir}/evidence/${tool}", mode:'copy'
    input:
    val tool
    path catalog, stageAs:'catalog_input'
    path contract, stageAs:'viral_evidence_tools.json'
    path helper, stageAs:'viral_evidence.py'
    output:
    path 'evidence', emit: evidence
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    /usr/bin/python3 viral_evidence.py --tool '${tool}' --catalog catalog_input --contract viral_evidence_tools.json --outdir evidence --cpus ${task.cpus}
    """
}

process VIRAL_VIBRANT {
    cpus 6
    memory '12 GB'
    maxForks 1
    publishDir "${params.outdir}/evidence/${tool}", mode:'copy'
    input:
    val tool
    path catalog, stageAs:'catalog_input'
    path contract, stageAs:'viral_evidence_tools.json'
    path helpers
    path genes, stageAs:'genes_input'
    path canonical_genes, stageAs:'canonical??.tsv'
    output:
    path 'evidence', emit: evidence
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    /usr/bin/python3 viral_evidence.py --tool '${tool}' --catalog catalog_input --contract viral_evidence_tools.json --genes-dir genes_input --canonical-genes ${canonical_genes.join(' ')} --outdir evidence --cpus ${task.cpus}
    """
}

process VIRAL_EVIDENCE_SUMMARY {
    cpus 1
    memory '2 GB'
    maxForks 1
    publishDir "${params.outdir}", mode:'copy'
    input:
    path catalog, stageAs:'catalog_input'
    path evidence_dirs, stageAs:'tool??/*'
    val enabled_tools
    path helpers
    output:
    path 'evidence_summary', emit: summary
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 viral_evidence_summary.py --catalog catalog_input --evidence-dirs ${evidence_dirs.join(' ')} --enabled-tools '${enabled_tools}' --outdir evidence_summary
    """
}

process VIRAL_VCONTACT3 {
    cpus 6
    memory '12 GB'
    maxForks 1
    publishDir "${params.outdir}/evidence/${tool}", mode:'copy'
    input:
    val tool
    path catalog, stageAs:'catalog_input'
    path contract, stageAs:'viral_evidence_tools.json'
    path helper, stageAs:'viral_evidence.py'
    output:
    path 'evidence', emit: evidence
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    /usr/bin/python3 viral_evidence.py --tool '${tool}' --catalog catalog_input --contract viral_evidence_tools.json --outdir evidence --cpus ${task.cpus}
    """
}
