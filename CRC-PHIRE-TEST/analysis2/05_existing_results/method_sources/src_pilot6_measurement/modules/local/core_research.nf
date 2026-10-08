process CORE_RESEARCH_CATALOG {
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/catalog", mode:'copy'
    input:
    path master, stageAs:'baseline_master.tsv'
    path fasta, stageAs:'candidates.fna'
    path blast
    path qc
    path element_class, stageAs:'element_class.tsv'
    path helper, stageAs:'core_votu.py'
    path catalog_helper, stageAs:'catalog_policy.py'
    path common, stageAs:'mainline_common.py'
    path policy, stageAs:'catalog_policy.json'
    path evidence, stageAs:'catalog_evidence'
    output:
    path 'catalog', emit: catalog
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 core_votu.py --master '${master}' --fasta '${fasta}' --blast '${blast}/alignments.tsv' \\
      --qc '${qc}/genome_qc.tsv' --catalog-policy '${policy}' --catalog-evidence '${evidence}' --element-class '${element_class}' --require-element-evidence --outdir catalog
    """
}

process CORE_LOCUS_EVIDENCE {
    cpus 6
    memory '8 GB'
    publishDir "${params.outdir}/locus", mode:'copy'
    input:
    path prepared_dirs, stageAs:'prepared??/*'
    path catalog
    path policy, stageAs:'catalog_policy.json'
    path evidence, stageAs:'catalog_evidence'
    path helper, stageAs:'core_locus.py'
    path catalog_helper, stageAs:'catalog_policy.py'
    path votu_helper, stageAs:'core_votu.py'
    path common, stageAs:'mainline_common.py'
    output:
    path 'locus', emit: locus
    script:
    def pp = (prepared_dirs instanceof List ? prepared_dirs : [prepared_dirs]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 core_locus.py --catalog '${catalog}' --policy '${policy}' --evidence '${evidence}' \\
      --prepared-dirs ${pp} --cpus ${task.cpus} --outdir locus
    """
}

process CORE_RESEARCH_SUMMARY {
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/research", mode:'copy'
    input:
    path catalog
    path locus
    path qc
    path gene_dirs, stageAs:'genes??/*'
    path function_dirs, stageAs:'functions??/*'
    path function_summary
    path kofam_summary
    path whole_genome_evidence
    val whole_enabled
    path helper, stageAs:'core_research_summary.py'
    path catalog_helper, stageAs:'catalog_policy.py'
    path common, stageAs:'mainline_common.py'
    path study_schema, stageAs:'study_schema.json'
    path callability_helper, stageAs:'genome_callability.py'
    output:
    path 'research', emit: summary
    script:
    def genes = (gene_dirs instanceof List ? gene_dirs : [gene_dirs]).collect { "'${it}'" }.join(' ')
    def functions = (function_dirs instanceof List ? function_dirs : [function_dirs]).collect { "'${it}'" }.join(' ')
    def ko = params.enable_kofam ? "--kofam-summary '${kofam_summary}'" : ''
    def wg = whole_enabled ? "--whole-genome-evidence '${whole_genome_evidence}'" : ''
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 core_research_summary.py --catalog '${catalog}' --locus '${locus}' --qc '${qc}' \\
      --gene-dirs ${genes} --function-dirs ${functions} --function-summary '${function_summary}' \\
      ${ko} ${wg} --study-schema '${study_schema}' --outdir research
    """
}

process CORE_WHOLE_GENOME_EVIDENCE {
    cpus 6
    memory '8 GB'
    maxForks 1
    publishDir "${params.outdir}/whole_genome_evidence", mode:'copy'
    input:
    path prepared_dirs, stageAs:'prepared??/*'
    path catalog
    path config
    path helpers
    output:
    path 'whole_genome', emit: evidence
    script:
    def pp = (prepared_dirs instanceof List ? prepared_dirs : [prepared_dirs]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 whole_genome_evidence.py --catalog '${catalog}' --config '${config}' \\
      --prepared-dirs ${pp} --cpus ${task.cpus} --outdir whole_genome
    """
}
