process CORE_PHISPY {
    tag "${meta.id}"
    cpus 6
    memory '12 GB'
    maxForks 1
    container 'sha256:1ad724df938c3e52641769fcc642e729cd728c22bbb4a23d8f4f7782e24debb1'
    publishDir "${params.outdir}/callers/${meta.id}/phispy", mode:'copy'
    input:
    tuple val(meta), path(prepared), path(genbank)
    path helper, stageAs:'core_callers.py'
    output:
    tuple val(meta), path('phispy'), emit: result
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 core_callers.py run-phispy --genome-id '${meta.id}' \
        --fasta '${prepared}/genome.fna' --contig-map '${prepared}/contig_map.tsv' \
        --genbank '${genbank}' --threads ${task.cpus} --outdir phispy
    """
}

process CORE_PHAGEBOOST {
    tag "${meta.id}"
    cpus 6
    memory '12 GB'
    maxForks 1
    container 'sha256:1ad724df938c3e52641769fcc642e729cd728c22bbb4a23d8f4f7782e24debb1'
    publishDir "${params.outdir}/callers/${meta.id}/phageboost", mode:'copy'
    input:
    tuple val(meta), path(prepared)
    path helper, stageAs:'core_callers.py'
    output:
    tuple val(meta), path('phageboost'), emit: result
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 core_callers.py run-phageboost --genome-id '${meta.id}' \
        --fasta '${prepared}/genome.fna' --contig-map '${prepared}/contig_map.tsv' \
        --threads ${task.cpus} --outdir phageboost
    """
}

process CORE_PHISPY_UNASSESSED {
    tag "${meta.id}:missing_GenBank_annotations"
    cpus 1
    memory '1 GB'
    container 'sha256:1ad724df938c3e52641769fcc642e729cd728c22bbb4a23d8f4f7782e24debb1'
    publishDir "${params.outdir}/callers/${meta.id}/phispy", mode:'copy'
    input:
    tuple val(meta), path(prepared)
    path helper, stageAs:'core_callers.py'
    output:
    tuple val(meta), path('phispy'), emit: result
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 core_callers.py mark-unassessed --genome-id '${meta.id}' --outdir phispy
    """
}

process CORE_COMPARE_CALLERS {
    cpus 1
    memory '1 GB'
    container 'sha256:1ad724df938c3e52641769fcc642e729cd728c22bbb4a23d8f4f7782e24debb1'
    publishDir "${params.outdir}/caller_comparison", mode:'copy'
    input:
    path baseline, stageAs:'baseline.tsv'
    path pilot_dirs, stageAs:'pilot??/*'
    path helper, stageAs:'core_callers.py'
    output:
    path 'caller_comparison', emit: result
    script:
    def pilots = (pilot_dirs instanceof List ? pilot_dirs : [pilot_dirs]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 core_callers.py compare --baseline '${baseline}' \
        --pilot-dirs ${pilots} --outdir caller_comparison
    """
}
