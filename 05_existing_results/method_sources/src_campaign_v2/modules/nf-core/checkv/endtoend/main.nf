process CHECKV_ENDTOEND {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.3--pyhdfd78af_0':
        'quay.io/biocontainers/checkv:1.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path ("${prefix}/quality_summary.tsv") , emit: quality_summary
    tuple val(meta), path ("${prefix}/completeness.tsv")    , emit: completeness
    tuple val(meta), path ("${prefix}/contamination.tsv")   , emit: contamination
    tuple val(meta), path ("${prefix}/complete_genomes.tsv"), emit: complete_genomes
    tuple val(meta), path ("${prefix}/proviruses.fna")      , emit: proviruses
    tuple val(meta), path ("${prefix}/viruses.fna")         , emit: viruses
    tuple val("${task.process}"), val("checkv"), path("${prefix}/checkv.version.txt"), topic: versions, emit: versions_checkv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    checkv \\
        end_to_end \\
        $args \\
        -t $task.cpus \\
        -d $db \\
        $fasta \\
        $prefix
    # runtime_version_files_begin
    checkv -h 2>&1 | sed '1!d;s/^.*CheckV v//;s/:.*//' > ${prefix}/checkv.version.txt
    # runtime_version_files_end
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    mkdir -p ${prefix}
    touch ${prefix}/quality_summary.tsv
    touch ${prefix}/completeness.tsv
    touch ${prefix}/contamination.tsv
    touch ${prefix}/complete_genomes.tsv
    touch ${prefix}/proviruses.fna
    touch ${prefix}/viruses.fna
    echo 'stub_not_measured' > ${prefix}/checkv.version.txt
    """
}
