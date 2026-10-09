// Native pinned containers: no Docker daemon/socket exposed inside processes.
process CORE_GENES {
    publishDir "${params.outdir}/samples/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    cpus 1
    memory '2 GB'
    container 'community.wave.seqera.io/library/checkm2:1.1.0--60f287bc25d7a10d@sha256:fcd923f5d60786deede345f52e8ad3ec73954ca8ea2eb914b56e7d38c8e1159b'
    input:
    tuple val(meta), path(prepared)
    path(core_helper)
    path(phage_helper)
    output:
    tuple val(meta), path('genes'), emit: genes
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    mkdir -p genes/raw
    prodigal -i ${prepared}/genome.fna -p single -g 11 -f gff \\
      -o genes/raw/prodigal.gff -a genes/raw/prodigal.faa -d genes/raw/prodigal.ffn \\
      > genes/raw/prodigal.stdout.log 2> genes/raw/prodigal.stderr.log
    python3 ${core_helper} genes-normalize --genome-id '${meta.id}' \\
      --prepared ${prepared} --raw genes/raw --outdir genes
    """
}

process CORE_CHECKM2 {
    publishDir "${params.outdir}/samples/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    cpus 6
    memory '12 GB'
    maxForks 1
    container 'community.wave.seqera.io/library/checkm2:1.1.0--60f287bc25d7a10d@sha256:fcd923f5d60786deede345f52e8ad3ec73954ca8ea2eb914b56e7d38c8e1159b'
    input:
    tuple val(meta), path(prepared)
    path(database)
    path(metadata)
    val(policy)
    path(core_helper)
    path(phage_helper)
    output:
    tuple val(meta), path('qc'), emit: qc
    script:
    def metadataArg = metadata ? "--metadata ${metadata}" : ''
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    export TF_NUM_INTRAOP_THREADS=${task.cpus}
    export TF_NUM_INTEROP_THREADS=1
    export OMP_NUM_THREADS=${task.cpus}
    mkdir -p qc
    echo '1b86ef3eac0813c1853f53182c17657045e3763d66f384ec95747261a63ae46f  ${database}' > qc/database.sha256
    sha256sum -c qc/database.sha256 > qc/database_verification.log
    checkm2 predict --input ${prepared}/genome.fna --output-directory qc/raw \\
      --database_path ${database} --threads ${task.cpus} --lowmem \\
      > qc/checkm2.stdout.log 2> qc/checkm2.stderr.log
    python3 ${core_helper} qc-normalize --genome-id '${meta.id}' \\
      --prepared ${prepared} --checkm2-report qc/raw/quality_report.tsv --checkm2-name genome \\
      --completeness-min ${policy.completeness_min} --contamination-max ${policy.contamination_max} \\
      --high-completeness-min ${policy.high_completeness_min} --high-contamination-max ${policy.high_contamination_max} \\
      --fragmentation-n50 ${policy.fragmentation_n50} ${metadataArg} --outdir qc
    """
}

process CORE_FASTANI {
    tag "${meta.id}"
    cpus 6
    memory '2 GB'
    container 'quay.io/biocontainers/fastani:1.34--hb66fcc3_7@sha256:0723121d8154f7e4917304142824d12bff93494d2b2ed4ab072d51eccf6c8494'
    input:
    tuple val(meta), path(prepared)
    path(reference_dir)
    output:
    tuple val(meta), path(prepared), path('ani_raw'), emit: raw
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    mkdir -p ani_raw
    for reference in ${reference_dir}/*.fna; do
      [ -f "\$reference" ] || continue
      accession=\$(basename "\$reference" .fna)
      fastANI -q ${prepared}/genome.fna -r "\$reference" --threads ${task.cpus} \\
        --fragLen 3000 --minFraction 0.2 -o "ani_raw/\$accession.forward.tsv" \\
        > "ani_raw/\$accession.forward.stdout.log" 2> "ani_raw/\$accession.forward.stderr.log"
      fastANI -q "\$reference" -r ${prepared}/genome.fna --threads ${task.cpus} \\
        --fragLen 3000 --minFraction 0.2 -o "ani_raw/\$accession.reverse.tsv" \\
        > "ani_raw/\$accession.reverse.stdout.log" 2> "ani_raw/\$accession.reverse.stderr.log"
    done
    """
}

process CORE_TAXONOMY_NORMALIZE {
    publishDir "${params.outdir}/samples/${meta.id}", mode: 'copy'
    tag "${meta.id}"
    cpus 1
    memory '1 GB'
    container 'python:3.12.13-bookworm@sha256:3cd9086bdb30f7c9bc08a3fa621d9842e0d3f6f9291aeb4677e0547817c10b12'
    input:
    tuple val(meta), path(prepared), path(raw)
    path(reference_dir)
    path(metadata)
    val(policy)
    path(core_helper)
    path(phage_helper)
    output:
    tuple val(meta), path('taxonomy'), emit: taxonomy
    script:
    def metadataArg = metadata ? "--metadata ${metadata}" : ''
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    mkdir -p taxonomy/raw
    cp -a ${raw}/. taxonomy/raw/
    python3 ${core_helper} taxonomy-normalize --genome-id '${meta.id}' \\
      --prepared ${prepared} --raw taxonomy/raw --reference-dir ${reference_dir} \\
      --reference-manifest ${reference_dir}/reference_manifest.json ${metadataArg} \\
      --ani-min ${policy.ani_min} --fragment-fraction-min ${policy.fragment_fraction_min} --outdir taxonomy
    """
}
