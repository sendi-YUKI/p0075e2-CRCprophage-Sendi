process VIRAL_PREPARE {
    tag "${meta.id}"
    publishDir "${params.outdir}/assemblies/${meta.id}/input", mode:'copy'
    input:
    tuple val(meta), path(fasta, stageAs:'assembly.fna')
    path helpers
    output:
    tuple val(meta), path('prepared'), emit: prepared
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 viral_catalog.py prepare --assembly-id '${meta.id}' --fasta assembly.fna --expected-sha256 '${meta.sha256}' --outdir prepared
    """
}

process VIRAL_COLLECT {
    stageInMode 'copy'
    publishDir "${params.outdir}", mode:'copy'
    input:
    path sample_dirs, stageAs:'samples??/*'
    path prepared_dirs, stageAs:'prepared??/*'
    val manifest_b64
    path legacy_catalog, stageAs:'legacy_input'
    val use_legacy
    path legacy_sources, stageAs:'legacy_sources.json'
    val use_legacy_sources
    path legacy_fastas, stageAs:'legacy_source??.fna'
    path helpers
    output:
    path 'discovery', emit: catalog
    script:
    def legacyArgs = use_legacy ? '--legacy-catalog legacy_input' : ''
    def sourceArgs = use_legacy_sources ? "--legacy-source-manifest legacy_sources.json --legacy-source-fastas ${legacy_fastas.join(' ')}" : ''
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 -c "import base64;open('assemblies.json','wb').write(base64.b64decode('${manifest_b64}'))"
    python3 viral_catalog.py assemble --manifest assemblies.json --sample-dirs ${sample_dirs.join(' ')} --prepared-dirs ${prepared_dirs.join(' ')} --catalog-scope '${params.catalog_scope}' ${legacyArgs} ${sourceArgs} --outdir discovery
    """
}

process VIRAL_BLAST {
    cpus 6
    memory '12 GB'
    publishDir "${params.outdir}/alignment", mode:'copy'
    input:
    path catalog, stageAs:'catalog_input'
    output:
    path 'blast.tsv', emit: alignments
    path 'blast_version.txt', emit: version
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    blastn -version > blast_version.txt
    n=\$(grep -c '^>' catalog_input/viral_sequences.fna || true)
    if [ "\$n" -gt 1 ]; then
      makeblastdb -in catalog_input/viral_sequences.fna -dbtype nucl -out viraldb
      blastn -query catalog_input/viral_sequences.fna -db viraldb -out blast.tsv -outfmt '6 std qlen slen qseq sseq' -num_threads ${task.cpus} -evalue 1e-5 -max_target_seqs "\$n"
    else
      touch blast.tsv
    fi
    """
}

process VIRAL_CLUSTER {
    publishDir "${params.outdir}", mode:'copy'
    input:
    path catalog, stageAs:'catalog_input'
    path blast, stageAs:'blast.tsv'
    path helpers
    output:
    path 'catalog', emit: catalog
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 viral_catalog.py cluster --catalog catalog_input --blast blast.tsv --outdir catalog --ani-percent ${params.viral_ani_percent} --af-shorter-percent ${params.viral_af_percent}
    """
}

process VIRAL_GENES {
    cpus 1
    memory '3 GB'
    publishDir "${params.outdir}", mode:'copy'
    input:
    path catalog, stageAs:'catalog_input'
    path helpers
    path gene_tools, stageAs:'viral_gene_tools.json'
    output:
    path 'genes', emit: genes
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 viral_genes.py call --catalog catalog_input --tools-manifest viral_gene_tools.json --outdir genes
    """
}

process VIRAL_GENE_CROSSWALK {
    publishDir "${params.outdir}", mode:'copy'
    input:
    path catalog, stageAs:'catalog_input'
    path genes, stageAs:'genes_input'
    path canonical_genes, stageAs:'canonical??.tsv'
    path helpers
    output:
    path 'crosswalk', emit: crosswalk
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 viral_genes.py crosswalk --catalog catalog_input --genes-dir genes_input --canonical-genes ${canonical_genes.join(' ')} --outdir crosswalk
    """
}

process VIRAL_FUNCTIONS {
    cpus 6
    memory '12 GB'
    publishDir "${params.outdir}", mode:'copy'
    input:
    path catalog, stageAs:'catalog_input'
    path genes, stageAs:'genes_input'
    path database, stageAs:'functions_db'
    path tools_manifest, stageAs:'core_functions_tools.json'
    path helpers
    output:
    path 'functions', emit: functions
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 viral_functions.py run --catalog catalog_input --genes-dir genes_input --database-root functions_db --tools-manifest core_functions_tools.json --cpus ${task.cpus} --outdir functions
    """
}

process VIRAL_FUNCTION_SUMMARY {
    publishDir "${params.outdir}", mode:'copy'
    input:
    path catalog, stageAs:'catalog_input'
    path functions, stageAs:'functions_input'
    path helpers
    output:
    path 'summary', emit: summary
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 viral_functions.py summarize --catalog catalog_input --functions functions_input --outdir summary
    """
}
