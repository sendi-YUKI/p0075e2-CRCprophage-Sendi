// Production B1 adapters; metadata remain attached to every channel record.
def cohortQuote(value) { "'" + value.toString().replace("'", "'\\''") + "'" }

process COHORT_INDEX {
    publishDir "${params.outdir}/reference", mode:'symlink'
    cpus 3
    memory '4 GB'
    container params.cohort_mapping_container
    input:
    path(reference_dir)
    path(helpers)
    output:
    path('mapping_index'), emit:index
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 cohort_cli.py index --reference-dir ${cohortQuote(reference_dir)} --outdir mapping_index --threads ${task.cpus}
    """
}
process COHORT_MAPPABILITY {
    publishDir "${params.outdir}/reference", mode:'symlink'
    cpus 3
    memory '4 GB'
    container params.cohort_mapping_container
    input:
    path(index_dir)
    path(profile)
    path(helpers)
    output:
    path('mappability/mappability.json'), emit:summary
    path('mappability'), emit:evidence
    script:
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 cohort_cli.py mappability --index-dir ${cohortQuote(index_dir)} --profile ${cohortQuote(profile)} --outdir mappability --threads ${task.cpus}
    """
}
process COHORT_MEASURE {
    tag "${meta.id}"
    publishDir "${params.outdir}/units", mode:'symlink'
    cpus 3
    memory '4 GB'
    container params.cohort_mapping_container
    input:
    tuple val(meta), path(unit, stageAs:'unit.json'), path(read1, stageAs:'read1/*'), path(read2, stageAs:'read2/*')
    path(index_dir)
    path(profile)
    path(detections)
    path(mappability)
    path(helpers)
    output:
    tuple val(meta), path("unit_${meta.id}"), emit:results
    script:
    def second = read2 ? "--reads2 ${cohortQuote(read2 instanceof List ? read2[0] : read2)}" : ''
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 cohort_cli.py map --index-dir ${cohortQuote(index_dir)} --unit unit.json --profile ${cohortQuote(profile)} \
      --reads1 ${cohortQuote(read1)} ${second} --outdir unit_${meta.id}/raw --threads ${task.cpus}
    python3 cohort_cli.py measure --sam unit_${meta.id}/raw/alignments.qname.bam \
      --primary-sam unit_${meta.id}/raw/alignments.primary.qname.bam --reference-dir ${cohortQuote(index_dir)} --unit unit.json --profile ${cohortQuote(profile)} \
      --detections ${cohortQuote(detections)} --mappability ${cohortQuote(mappability)} --outdir unit_${meta.id}
    """
}
process COHORT_METAPHLAN {
    tag "${meta.id}"
    publishDir "${params.outdir}/host", mode:'copy'
    cpus 6
    memory { "${params.cohort_metaphlan_required_memory_gb} GB" }
    container params.cohort_metaphlan_container
    input:
    tuple val(meta), path(unit, stageAs:'unit.json'), path(read1, stageAs:'read1/*'), path(read2, stageAs:'read2/*')
    path(database, stageAs:'metaphlan_db')
    path(database_manifest)
    val(index_name)
    path(helpers)
    output:
    path("metaphlan_${meta.id}"), emit:results
    script:
    def second = read2 ? "--reads2 ${cohortQuote(read2 instanceof List ? read2[0] : read2)}" : ''
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 cohort_cli.py metaphlan --unit unit.json --reads1 ${cohortQuote(read1)} ${second} \
      --database ${cohortQuote(database)} --index ${cohortQuote(index_name)} \
      --database-manifest ${cohortQuote(database_manifest)} --threads ${task.cpus} --outdir metaphlan_${meta.id}
    """
}
process COHORT_REPORT {
    publishDir "${params.outdir}", mode:'copy'
    cpus 1
    memory '1 GB'
    container 'python:3.12.13-bookworm@sha256:3cd9086bdb30f7c9bc08a3fa621d9842e0d3f6f9291aeb4677e0547817c10b12'
    input:
    path(prepared)
    path(results)
    path(host_results)
    path(helpers)
    val(disabled_reason)
    output:
    path('summary'), emit:report
    script:
    def resultArgs = results.collect { cohortQuote(it) }.join(' ')
    def hostArgs = host_results.collect { cohortQuote(it) }.join(' ')
    """
    # native_runtime_identity=${task.ext.spark_native_identity ?: 'container-runtime'}
    python3 cohort_cli.py report --prepared ${cohortQuote(prepared)} --measurements ${resultArgs} \
      --metaphlan-results ${hostArgs} --metaphlan-disabled-reason ${cohortQuote(disabled_reason)} --outdir summary
    """
}
