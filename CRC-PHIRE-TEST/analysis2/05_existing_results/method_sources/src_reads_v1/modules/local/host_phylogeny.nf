process HOST_CORE_GFF {
    tag "${meta.id}"
    cpus 1
    memory '2 GB'
    input:
    tuple val(meta), path(genes), path(prepared)
    path host_config
    path helper
    path common
    output:
    tuple val(meta), path('host_input'), emit: adapted
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 ${helper} adapt --genome-id '${meta.id}' --genes-dir '${genes}' --prepared '${prepared}' \\
      --config '${host_config}' --outdir host_input
    """
}

process HOST_CORE_PANAROO {
    publishDir "${params.outdir}/host_phylogeny", mode: 'copy'
    cpus 6
    memory '46 GB'
    maxForks 1
    input:
    path(inputs), stageAs: 'inputs??/*'
    path qc
    path host_config
    path helper
    path common
    output:
    path 'host_clusters', emit: clusters
    script:
    def dirs = (inputs instanceof List ? inputs : [inputs]).collect { "'${it}'" }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 ${helper} cluster --inputs ${dirs} --qc '${qc}' --config '${host_config}' \\
      --cpus ${task.cpus} --outdir host_clusters
    """
}

process HOST_CORE_ALIGN {
    publishDir "${params.outdir}/host_phylogeny", mode: 'copy'
    cpus 6
    memory '12 GB'
    maxForks 1
    input:
    path clusters
    path host_config
    path helper
    path common
    output:
    path 'host_alignment', emit: alignment
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 ${helper} align --clusters '${clusters}' --config '${host_config}' \\
      --cpus ${task.cpus} --outdir host_alignment
    """
}

process HOST_TREE_PRELIMINARY {
    publishDir "${params.outdir}/host_phylogeny", mode: 'copy'
    cpus 6
    memory '46 GB'
    maxForks 1
    input:
    path alignment
    path host_config
    path helper
    path common
    output:
    path 'preliminary_tree', emit: tree
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 ${helper} tree --alignment '${alignment}' --config '${host_config}' \\
      --stage preliminary --cpus ${task.cpus} --outdir preliminary_tree
    """
}

process HOST_MOBILE_MASK {
    publishDir "${params.outdir}/host_phylogeny", mode: 'copy'
    cpus 1
    memory '8 GB'
    input:
    path alignment
    path preliminary_tree
    path master
    path mask_bundle
    path host_config
    path helper
    path common
    output:
    path 'host_mask', emit: masked
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 ${helper} mask --alignment '${alignment}' --preliminary-tree '${preliminary_tree}' \\
      --master '${master}' --mask-bundle '${mask_bundle}' --config '${host_config}' --outdir host_mask
    """
}

process HOST_TREE_FINAL {
    publishDir "${params.outdir}/host_phylogeny", mode: 'copy'
    cpus 6
    memory '46 GB'
    maxForks 1
    input:
    path alignment
    path host_config
    path helper
    path common
    output:
    path 'final_tree', emit: tree
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'unconfigured'}
    python3 ${helper} tree --alignment '${alignment}' --config '${host_config}' \\
      --stage final --cpus ${task.cpus} --outdir final_tree
    """
}
